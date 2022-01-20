package stream;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import cardinality.CardinalityTracker;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * This implementation allows multiple concurrent open streams.
 * 
 * @author Brian Bushnell
 * @date May 14, 2019
 *
 */
public class MultiCros3 extends BufferedMultiCros {
	
	/** 
	 * For testing.<br>
	 * args should be:
	 * {input file, output pattern, names...}
	 */
	public static void main(String[] args){
		
		//Do some basic parsing
		String in=args[0];
		String pattern=args[1];
		ArrayList<String> names=new ArrayList<String>();
		for(int i=2; i<args.length; i++){names.add(args[i]);}
		
		//Create an mcros for this pattern
		MultiCros3 mcros=new MultiCros3(pattern, null, false, false, false, false, FileFormat.FASTQ, false, 4);
		
		//Create the input stream
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, false, in);
		cris.start();
		
		//Fetch the first list
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		//Process all remaining lists
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			
			//Add the reads by barcode
			for(Read r1 : reads){
				mcros.add(r1, r1.barcode(true));
			}
			
			//Return the old list and get a new one
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		
		//Close the streams
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
		mcros.close();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** @See Details in superclass constructor */
	public MultiCros3(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, boolean threaded_, int maxStreams_){
		super(pattern1_, pattern2_, overwrite_, append_, allowSubprocess_, useSharedHeader_, defaultFormat_, threaded_, maxStreams_);
		
		bufferMap=new LinkedHashMap<String, Buffer>();
		streamQueue=new ArrayDeque<String>(maxStreams);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public boolean finishedSuccessfully(){
		return !errorState;
	}
	
	@Override
	public void add(Read r, String name){
		Buffer b=bufferMap.get(name);
		if(b==null){
			b=new Buffer(name);
			bufferMap.put(name, b);
			//Note: I could adjust bytesPerBuffer threshold here in response to the number of buffers.
		}
		b.add(r);
	}
	
	@Override
	public long dumpResidual(ConcurrentReadOutputStream rosu){
		//For each Buffer, check if it contains residual reads
		//If so, dump it into the stream
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			Buffer b=e.getValue();
			assert((b.readsIn<minReadsToDump) == (b.list!=null && !b.list.isEmpty()));
			if(b.readsIn>0 && b.readsIn<minReadsToDump){
				assert(b.list!=null && !b.list.isEmpty());
				residualReads+=b.readsIn;
				residualBases+=b.basesIn;
				if(rosu!=null){rosu.add(b.list, 0);}
			}
			b.list=null;
		}
		return residualReads;
	}
	
	@Override
	public ByteBuilder report(){
		ByteBuilder bb=new ByteBuilder(1024);
		
		//Add a line for residual reads dumped 
		if(minReadsToDump>0){
			bb.append("Residual").tab().append(residualReads).tab().append(residualBases).nl();
		}
		
		//Add a line for each Buffer
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			Buffer buffer=e.getValue();
			if(buffer.numDumps>0){//Only add a line if this Buffer actually created a file
				buffer.appendTo(bb);
			}
		}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	long closeInner() {
		//First dump everything
		final long x=dumpAll();
		//Then, retire any active streams
		while(!streamQueue.isEmpty()){retire();}
		return x;
	}
	
	@Override
	long dumpAll(){
		long dumped=0;
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			dumped+=e.getValue().dump();
		}
		return dumped;
	}
	
	/** Close the least-recently-used stream */
	private void retire(){
		if(verbose){System.err.println("Enter retire(); streamQueue="+streamQueue);}
		
		//Select the first name in the queue, which is the least-recently-used.
		String name=streamQueue.removeFirst();
		Buffer b=bufferMap.get(name);
		
		if(verbose){System.err.println("retire("+name+"); ros="+(b.currentRos!=null)+", streamQueue="+streamQueue);}
		assert(b!=null);
		assert(b.currentRos!=null) : name+"\n"+streamQueue+"\n";
		
		//Purge any remaining reads first
		b.dump(b.currentRos);
		
		//Then close the stream
		if(closeFast){
			b.currentRos.close();//Faster, but the error state is not caught
			//Also the stream could be re-opened before writing is done, which is very bad, so this is unsafe as implemented.
		}else{
			errorState=ReadWrite.closeStream(b.currentRos) | errorState;//Traditional synchronous close-and-wait
		}
		
		//Delete the pointer to output stream
		b.currentRos=null;
		if(verbose){System.err.println("Exit retire("+name+"); ros="+(b.currentRos!=null)+", streamQueue="+streamQueue);}
//		assert(!streamQueue.contains(name)); //Slow
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * A Buffer holds reads destined for to a specific file.
	 * When sufficient reads are present, it opens a stream and writes them.
	 * If too many streams are open, it closes another stream first.
	 */
	private class Buffer {
		
		Buffer(String name_){
			name=name_;
			String s1=pattern1.replaceFirst("%", name);
			String s2=pattern2==null ? null : pattern2.replaceFirst("%", name);
			
			//These are created with overwrite=false append=true because 
			//the files will be appended to if the stream gets prematurely retired.
			//Therefore, files must be explicitly deleted first.
			//Alternative would be to create a new FileFormat each time.
			ff1=FileFormat.testOutput(s1, defaultFormat, null, allowSubprocess, false, true, false);
			ff2=FileFormat.testOutput(s2, defaultFormat, null, allowSubprocess, false, true, false);
			
			list=new ArrayList<Read>(readsPerBuffer);
			if(trackCardinality){loglog=CardinalityTracker.makeTracker();}
			if(verbose){System.err.println("Made buffer for "+name);}
		}
		
		/** 
		 * Add a read to this buffer, and update all the tracking variables.
		 * This may trigger a dump.
		 */
		void add(Read r){
			//Add the read
			list.add(r);
			
			//Gather statistics
			long size=r.countPairBytes();
			int count=r.pairCount();
			currentBytes+=size;
			bytesInFlight+=size;
			basesIn+=size;
			readsInFlight+=count;
			readsIn+=count;
			if(trackCardinality){loglog.hash(r);}
			
			//Decide whether to dump
			handleLoad();
		}
		
		/**
		 * Determine whether to dump this buffer based on its current size.
		 * Then, determine whether to dump all buffers based on their combined size.
		 */
		private void handleLoad(){
			//3rd term allows preemptive dumping
			//More generally, this triggers a dump if the reads in this buffer exceed the maximum allowed reads or bytes 
			if(list.size()>=readsPerBuffer || currentBytes>=bytesPerBuffer || (currentRos!=null && list.size()>=400)){
				if(verbose){
					System.err.println("list.size="+list.size()+"/"+readsPerBuffer+
							", bytes="+currentBytes+"/"+bytesPerBuffer+", bytesInFlight="+bytesInFlight+"/"+memLimit);
				}
				dump();
			}
			
			//Too much buffered data in ALL buffers; dump everything.
			if(bytesInFlight>=memLimit){
				long dumped=dumpAll();
				if(dumped<1 && Shared.EA()){//Dump failed; exit
					KillSwitch.kill("\nThis program ran out of memory."
							+ "\nTry increasing the -Xmx flag or get rid of the minreads flag,"
							+ "\nor disable assertions to skip this message and try anyway.");
				}
			}
		}
		
		/** Dump buffered reads, creating a stream if needed */
		long dump(){
			if(list.isEmpty() || readsIn<minReadsToDump){return 0;}
			ConcurrentReadOutputStream ros=getStream();
			return dump(ros);
		}
		
		/** 
		 * Dump buffered reads to the stream.
		 * If the buffer is empty, nothing happens. */
		long dump(final ConcurrentReadOutputStream ros){
			if(verbose){System.err.println("Dumping "+name);}
			if(list.isEmpty()){return 0;}
			final long size0=list.size();
			
			//Send the list to the output stream
			ros.add(list, numDumps);
			
			//Create a new list, since the old one is busy
			list=new ArrayList<Read>(readsPerBuffer);
			
			//Manage statistics
			bytesInFlight-=currentBytes;
			readsInFlight-=size0;
			readsWritten+=size0;
			currentBytes=0;
			numDumps++;
			return size0;
		}
		
		/** Fetch the stream for this buffer, creating a new one if needed */
		private ConcurrentReadOutputStream getStream(){
			if(verbose){System.err.println("Enter getStream("+name+"); ros="+(currentRos!=null)+", +streamQueue="+streamQueue);}
			
			if(currentRos!=null){//The stream already exists
//				assert(streamQueue.contains(name));//slow
				if(streamQueue.peekLast()!=name){
					//Move to end to prevent early retirement
					boolean b=streamQueue.remove(name);
					assert(b) : "streamQueue did not contain "+name+", but the ros was open.";
					streamQueue.addLast(name);
				}
			}else{//The stream does not exist, so create it
				createStream();
			}
			
			assert(currentRos!=null) : "The stream for "+name+" was not created.";
			assert(streamQueue.peekLast()==name) : "The stream for "+name+" was not placed in the queue.";
			if(verbose){System.err.println("Exit  getStream("+name+"); ros="+(currentRos!=null)+", +streamQueue="+streamQueue);}
			return currentRos;
		}
		
		/** Create a stream for this buffer, and stick it in the queue */
		private ConcurrentReadOutputStream createStream(){
			assert(currentRos==null) : "This should never be called if there is an existing stream.";
			if(numDumps==0 && overwrite){
				//First time, an existing file must be deleted first, because the ff is set to append mode
				if(verbose){System.err.println("Deleting "+name+" ; exists? "+ff1.exists());}
				delete(ff1);
				delete(ff2);
			}
			
//			assert(!streamQueue.contains(name));//slow
			
			//Active streams should never exceed maxStreams
			assert(streamQueue.size()<=maxStreams) : "Too many streams: "+streamQueue+", "+maxStreams;
			if(streamQueue.size()>=maxStreams){
				//Too many open streams; retire one.
				retire();
			}
			//After retirement, open streams must be less than maxStreams
			assert(streamQueue.size()<maxStreams) : "Too many streams: "+streamQueue+", "+maxStreams;
			
			//Create a stream
			currentRos=ConcurrentReadOutputStream.getStream(ff1, ff2, rswBuffers, null, useSharedHeader && numDumps==0);
			currentRos.start();
			if(verbose){System.err.println("Created ros "+name+"; ow="+ff1.overwrite()+", append="+ff1.append());}
			
			//Add it to the queue
			streamQueue.addLast(name);
			return currentRos;
		}
		
		/** Delete this file if it exists */
		private void delete(FileFormat ff){
			if(ff==null){return;}
			assert(overwrite || !ff.exists()) : "Trying to delete file "+ff.name()+", but overwrite=f.  Please add the flag overwrite=t.";
			ff.deleteIfPresent();
		}
		
		/** 
		 * Format this buffer's summary as a line of text.
		 * @param bb ByteBuilder to append the text
		 * @return The modified ByteBuilder
		 */
		ByteBuilder appendTo(ByteBuilder bb) {
			bb.append(name).tab().append(readsIn).tab().append(basesIn);
			if(trackCardinality){bb.tab().append(loglog.cardinality());}
			return bb.nl();
		}
		
		@Override
		public String toString(){
			return appendTo(new ByteBuilder()).toString();
		}
		
		/** Stream name, which is the variable part of the file pattern */
		private final String name;
		/** Output file 1 */
		private final FileFormat ff1;
		/** Output file 2 */
		private final FileFormat ff2;
		
		/** Once created, the stream sticks around to be re-used unless it is retired. */
		private ConcurrentReadOutputStream currentRos;
		
		/** Current list of buffered reads */
		private ArrayList<Read> list;
		
		/** Number of reads entering the buffer */
		private long readsIn=0;
		/** Number of bases entering the buffer */
		private long basesIn=0;
		/** Number of reads written to disk */
		@SuppressWarnings("unused")
		private long readsWritten=0;//This does not count read2!
		/** Number of bytes currently in this buffer (estimated) */
		private long currentBytes=0;
		/** Number of dumps executed */
		private long numDumps=0;
		/** Optional, for tracking cardinality */
		private CardinalityTracker loglog;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Open stream names */
	private final ArrayDeque<String> streamQueue;
	
	/** Map of names to buffers */
	public final LinkedHashMap<String, Buffer> bufferMap;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Trigger stream close, but don't wait for it to finish.
	 * Prevents error state from being captured.
	 * 
	 * THIS IS UNSAFE AND CAUSES ERRORS,
	 * because the stream might get reopened again before writing is finished.
	 */
	private static final boolean closeFast=false;

}
