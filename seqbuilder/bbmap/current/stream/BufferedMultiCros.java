package stream;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Allows output of reads to multiple different output streams.
 * Each output stream is controlled by a buffer,
 * which stores reads until there is a sufficient quantity to dump.
 * 
 * @author Brian Bushnell
 * @date May 14, 2019
 *
 */
public abstract class BufferedMultiCros extends Thread {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Details in primary constructor */
	public BufferedMultiCros(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, boolean threaded_){
		this(pattern1_, pattern2_, overwrite_, append_, allowSubprocess_, useSharedHeader_, defaultFormat_, threaded_, DEFAULT_MAX_STREAMS);
	}
	
	/**
	 * Primary constructor.
	 * @param pattern1_ Name pattern for file 1; must contain % (required)
	 * @param pattern2_ Name pattern for file 2; must contain % (optional)
	 * @param overwrite_ Permission to overwrite
	 * @param append_ Permission to append to existing files (this should generally be false)
	 * @param allowSubprocess_ Allow subprocesses such as pigz, bgzip, or samtools
	 * @param useSharedHeader_ Print the stored header (from an input sam file) in all output sam files 
	 * @param defaultFormat_ Assume files are in this format if they don't have a valid extension
	 * @param threaded_ Run this mcros in its own thread
	 * @param maxStreams_ Max allowed number of concurrent open streams
	 */
	public BufferedMultiCros(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, boolean threaded_, int maxStreams_){
		assert(pattern1_!=null && pattern1_.indexOf('%')>=0);
		assert(pattern2_==null || pattern1_.indexOf('%')>=0);
		
		//Perform # expansion for twin files
		if(pattern2_==null && pattern1_.indexOf('#')>=0){
			pattern1=pattern1_.replaceFirst("#", "1");
			pattern2=pattern1_.replaceFirst("#", "2");
		}else{
			pattern1=pattern1_;
			pattern2=pattern2_;
		}
		
		overwrite=overwrite_;
		append=append_;
		allowSubprocess=allowSubprocess_;
		useSharedHeader=useSharedHeader_;
		
		defaultFormat=defaultFormat_;
		
		threaded=threaded_;
		transferQueue=threaded ? new ArrayBlockingQueue<ArrayList<Read>>(8) : null;
		maxStreams=maxStreams_;
		
		memLimit=Tools.max(10000000, (long)(0.75*Shared.memAvailable()));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** True if no errors were encountered */
	public abstract boolean finishedSuccessfully();

	/** 
	 * Add a single read.  Should not be used in threaded mode.
	 * @param r Read to add.
	 * @param name Name of destination buffer.
	 */
	public abstract void add(Read r, String name);
	
	/** 
	 * Dump all buffered reads to disk, except when minReadsToDump forbids it.
	 * @return Number of reads dumped.
	 */
	abstract long dumpAll();
	
	/** 
	 * Dump all residual reads to this stream.
	 * @param rosu Destination stream.
	 * @return Number of residual reads dumped.
	 */
	public abstract long dumpResidual(ConcurrentReadOutputStream rosu);
	
	/** Dump everything and close any open streams. */
	abstract long closeInner();
	
	/** Generate a report on how many reads went to each file */
	public abstract ByteBuilder report();
	
	/*--------------------------------------------------------------*/
	/*----------------        Final Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Shut this down and perform any cleanup needed. */
	public final void close(){
		if(threaded){poisonAndWait();}
		else{closeInner();}
	}
	
	/** Primary file pattern */
	public final String fname(){return pattern1;}
	
	/** Return true if this stream has detected an error */
	public final boolean errorState(){
		return errorState;
	}
	
	/** 
	 * Send a list of reads to an output buffer.
	 * The reads must have a name attached to the object field in order to be written. 
	 */
	public final void add(ArrayList<Read> list) {
		if(threaded){//Send to the transfer queue
			try {
				transferQueue.put(list);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				KillSwitch.kill();
			}
		}else{//Add the reads from this thread
			addToBuffers(list);
		}
	}
	
	/** Send individual reads to their designated buffer */
	private final void addToBuffers(ArrayList<Read> list) {
		for(Read r : list){
			if(r.obj!=null){
				String name=(String)r.obj;
				add(r, name);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Threaded Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	/** For threaded mode */
	public final void run(){
		assert(threaded) : "This should only be called in threaded mode.";
		try {
			for(ArrayList<Read> list=transferQueue.take(); list!=poisonToken; list=transferQueue.take()){
				if(verbose){System.err.println("Got list; size=\"+transferQueue.size())");}
				addToBuffers(list);
				if(verbose){System.err.println("Added list; size="+transferQueue.size());}
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
			//Terminate JVM if something goes wrong
			KillSwitch.kill();
		}
		closeInner();
	}
	
	/** Indicate that no more reads will be sent, for threaded mode */
	public final void poison(){
		assert(threaded) : "This should only be called in threaded mode.";
		transferQueue.add(poisonToken);
	}
	
	/** Indicate that no more reads will be sent, for threaded mode */
	public final void poisonAndWait(){
		assert(threaded) : "This should only be called in threaded mode.";
		poison();
		waitForFinish();
	}
	
	/** Wait for this object's thread to terminate */
	public final void waitForFinish(){
		assert(threaded);
		if(verbose){System.err.println("Waiting for finish.");}
		while(this.getState()!=Thread.State.TERMINATED){
			if(verbose){System.err.println("Attempting join.");}
			try {
				this.join(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output file patterns containing a % symbol */
	public final String pattern1, pattern2;
	
	/** True if an error was encountered */
	boolean errorState=false;
	
	/** File overwrite permission */
	final boolean overwrite;
	
	/** File append permission */
	final boolean append;
	
	/** Subprocess spawning permission (e.g., for pigz) */
	final boolean allowSubprocess;
	
	/** Output file format, if unclear from file extension */
	final int defaultFormat;
	
	/** Buffers for each ReadStreamWriter */
	int rswBuffers=1;
	
	/** Print the shared header (for sam files) */
	final boolean useSharedHeader;
	
	/** Dump everything if this limit is reached from buffered reads */
	long memLimit;

	/** Allow this many active streams, for MCros3 */
	public final int maxStreams;
	
	/** Dump a buffer once it holds this many reads */
	public int readsPerBuffer=2000;
	
	/** Dump a buffer once it holds this many bytes (estimated) */
	public int bytesPerBuffer=3000000;
	
	/** Never write files with fewer than this many reads */
	public long minReadsToDump=0;

	/** Number of reads encountered that were not written */
	public long residualReads=0, residualBases=0;
	
	/** Current number of buffered reads */
	long readsInFlight=0;
	
	/** Current number of buffered bytes (estimated) */
	long bytesInFlight=0;
	
	/** Used when MultiCros is run in threaded mode */
	private final ArrayBlockingQueue<ArrayList<Read>> transferQueue;
	
	/** Signal to terminate when in threaded mode */
	private final ArrayList<Read> poisonToken=new ArrayList<Read>(0);
	
	/** True if this object is intended to run in a separate thread */
	public final boolean threaded;
	
	/** Use a LogLog to track cardinality for each output file */
	public boolean trackCardinality=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int DEFAULT_MAX_STREAMS=4;
	public static boolean verbose=false;

}
