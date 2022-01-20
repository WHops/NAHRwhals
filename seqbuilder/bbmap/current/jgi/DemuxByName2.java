package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.BufferedMultiCros;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.MultiCros2;
import stream.MultiCros3;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;
import structures.ListNum;


/**
 * This class is designed to handle very large numbers of output files
 * with a fixed number of file handles.  An example use case is demultiplexing
 * Illumina Novaseq runs.
 * 
 * @author Brian Bushnell
 * @date May 1, 2019
 *
 */
public class DemuxByName2 {

	/** Code entrance from the command line */
	public static void main(String[] args){
		
		//Capture values of static variables that might be modified in case this is called by another class.
		final int oldCap=Shared.numBuffers(), oldZipThreads=ReadWrite.MAX_ZIP_THREADS, oldZl=ReadWrite.ZIPLEVEL;
		
		//External compressor selection
		final boolean oldPigz=ReadWrite.USE_PIGZ, oldUnpigz=ReadWrite.USE_UNPIGZ;
		final boolean oldBgzip=ReadWrite.USE_BGZIP, oldPreferBgzip=ReadWrite.PREFER_BGZIP;
		
		//Preserve quality scores even for reads with with incorrect quality scores
		final boolean oldCQ=Read.CHANGE_QUALITY;
		
		//Create and start a timer, for statistics
		Timer t=new Timer();
		
		//Create the demultiplexer object
		DemuxByName2 demultiplexer=new DemuxByName2(args);
		
		//Process all data
		demultiplexer.process(t);
		
		//Restore values of static variables.
		Shared.setBuffers(oldCap);
		ReadWrite.ZIPLEVEL=oldZl;
		ReadWrite.USE_PIGZ=oldPigz;
		ReadWrite.USE_BGZIP=oldBgzip;
		ReadWrite.PREFER_BGZIP=oldPreferBgzip;
		ReadWrite.USE_UNPIGZ=oldUnpigz;
		ReadWrite.MAX_ZIP_THREADS=oldZipThreads;
		Read.CHANGE_QUALITY=oldCQ;
		
		//Close the print stream if it was redirected
		Shared.closeStream(demultiplexer.outstream);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public DemuxByName2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set some static variables
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=true;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		SamLine.SET_FROM_OK=true;
		
		//Reduce default compression level to increase speed
		//This automatically gets bumped up to 4 if pigz/bgzip is detected
		ReadWrite.ZIPLEVEL=2;
		Read.CHANGE_QUALITY=false;
		
		
		{//Argument-parsing block
			
			//Parse arguments
			final Parser parser=parse(args);
			
			//Process parser fields
			//This is in the constructor rather than a function to allow final fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			trackCardinality=parser.loglog;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		//Ensure the input and output files are specified correctly
		validate();

		//Create FileFormat opjects for input files
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

		//If the input and output are in sam/bam format, print the sam header
		if(ffin1!=null && out1!=null && ffin1.samOrBam()){
			String ext=ReadWrite.rawExtension(out1);
			useSharedHeader=FileFormat.isSamOrBam(ext);
		}
		
		//Set up data structures
		handleModeAndNames();
	}
	
	/** 
	 * Parse command-line arguments
	 * @param args Command-line arguments
	 * @return A Parser object with fields filled from the parameters
	 */
	private Parser parse(final String[] args){
		
		//Handles parsing of common flags
		Parser parser=new Parser();
		parser.overwrite=overwrite;
		
		//Handles parsing of flags specific to this program
		for(int i=0; i<args.length; i++){//Parsing loop
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				BufferedMultiCros.verbose=verbose;
			}else if(a.equals("names") || a.equals("name") || a.equals("affixes")){
				if(b!=null){
					String[] x=b.split(",");
					for(String s : x){
						names.put(s, s);
					}
				}
			}else if(a.equals("threaded")){
				threaded=Parse.parseBoolean(b);
			}else if(a.equals("chrom") || a.equals("chrommode") || a.equals("scaffold")){
				if(Parse.parseBoolean(b)){mode=CHROM_MODE;}
			}else if(a.equals("substringmode") || a.equals("substring")){
				if(Parse.parseBoolean(b)){mode=SUBSTRING_MODE;}
			}else if(a.equals("barcode") || a.equals("barcodemode") || a.equals("index")){
				if(Parse.parseBoolean(b)){mode=BARCODE_MODE;}
			}else if(a.equals("perheader") || a.equals("persequence") || a.equals("header") || a.equals("headermode")){
				if(Parse.parseBoolean(b)){mode=HEADER_MODE;}
			}else if(a.equals("affix") || a.equals("affixmode")){
				if(Parse.parseBoolean(b)){mode=AFFIX_MODE;}
			}else if(a.equals("delimiter")){
				delimiter=Parse.parseSymbol(b);
				if(delimiter!=null){mode=DELIMITER_MODE;}
			}else if(a.equals("prefixmode") || a.equals("prefix") || a.equals("pm")){
				prefixMode=Parse.parseBoolean(b);
			}else if(a.equals("suffixmode") || a.equals("suffix") || a.equals("sm")){
				prefixMode=!Parse.parseBoolean(b);
			}else if(a.equals("column")){
				column=Integer.parseInt(b);
				assert(column>0 || column==-1) : "Column is 1-based; must be 1+ or else -1 to disable.";
				column--;
			}else if(a.equalsIgnoreCase("length") || a.equalsIgnoreCase("len") || a.equalsIgnoreCase("affixlength") || a.equalsIgnoreCase("affixlen")){
				fixedLength=Integer.parseInt(b); //also set affix mode, perhaps
			}else if(a.equalsIgnoreCase("stats") || a.equalsIgnoreCase("report") || a.equalsIgnoreCase("results")){
				stats=b;
			}else if(a.equalsIgnoreCase("minreadstodump") || a.equalsIgnoreCase("minreads")){
				minReadsToDump=Parse.parseKMG(b);
			}else if(a.equalsIgnoreCase("mcrostype")){
				mcrosType=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("streams") || a.equalsIgnoreCase("maxstreams")){
				maxStreams=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("readsperbuffer") || a.equalsIgnoreCase("rpb")){
				readsPerBuffer=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("bytesPerBuffer") || a.equalsIgnoreCase("bpb")){
				bytesPerBuffer=Parse.parseIntKMG(b);
			}else if(a.equals("hdist") || a.equals("hamming") || a.equals("hammingdistance")){
				hdist=Integer.parseInt(b);
			}else if(a.equals("outu") || a.equals("outu1")){
				outu1=b;
			}else if(a.equals("outu2")){
				outu2=b;
			}else if(a.equals("pattern")){
				parser.out1=b;
			}else if(a.equals("remap") || a.equals("replace")){
				remap=Parse.parseRemap(b);
			}
			//Old version, allowed more complicated replacements
//			else if(a.equals("replace")){
//				if(b==null){
//					replaceA=replaceB=null;
//				}else{
//					if(b.length()==2){
//						replaceA=b.substring(0, 1);
//						replaceB=b.substring(1);
//					}else{
//						String[] pair=b.split(",");
//						replaceA=Parse.parseSymbol(pair[0]);
//						replaceB=Parse.parseSymbol(pair[1]);
//					}
//				}
//			}
			else if(parser.parse(arg, a, b)){
				//Flag was captured by the parser; do nothing
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/** 
	 * Ensure input and output files are specified correctly,
	 * and that interleaved mode is consistent with the number of input files.
	 */
	private void validate(){
		
		//Validate patterns
		assert(out1==null || out1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(out2==null || out2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout1==null || qfout1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout2==null || qfout2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";

		//Perform # replacement for twin files
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(outu1!=null && outu2==null && outu1.indexOf('#')>-1){
			outu2=outu1.replace("#", "2");
			outu1=outu1.replace("#", "1");
		}

		//Disable interleaving if in2 is specified
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		{//Perform various file validation
			assert(FastaReadInputStream.settingsOK());
			if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
			if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
			
			//ByteFile2 is multithreaded; only use it if there are plenty of threads
			if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
				ByteFile.FORCE_MODE_BF2=true;
			}
			
			//Do some interleaving-detection logic if unspecified
			if(!setInterleaved){
				assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
				if(in2!=null){ //If there are 2 input streams.
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}else{ //There is one input stream.
					if(out2!=null){
						FASTQ.FORCE_INTERLEAVED=true;
						FASTQ.TEST_INTERLEAVED=false;
						outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
					}
				}
			}
			
			//If the user explicitly used the word "null", change it to null.
			if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
			if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
			
			//Ensure output files have no duplicates and can be written
			//Unfortunately, this does not prevent patterns expanding into a matching string
			if(!Tools.testOutputFiles(overwrite, append, false, outu1, outu2, stats)){
				outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2+", "+outu1+", "+outu2+", "+stats);
				throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+", "+outu1+", "+outu2+", "+stats+"\n");
			}

			//This does not prevent expansion into a matching string
			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
		}
	}
	
	/**
	 * Perform necessary setup of variables and data structures
	 * related to the processing mode and the names table. 
	 */
	private void handleModeAndNames(){
		
		//Handle names
		if(names.isEmpty()){
			names=null;
		}else{
			{//Process names, because they can either be literals or filenames at this point
				String[] keys=names.keySet().toArray(new String[names.size()]);
				names.clear();
				for(String key : keys){//Process either literal names or files of names
					File f=new File(key);
					if(f.exists() && f.isFile()){//Read the file
						TextFile tf=new TextFile(key);
						String[] lines=tf.toStringLines();
						for(String line : lines){names.put(line, line);}//Add each line as a name
					}else{names.put(key, key);}//Re-add each key as a literal name
				}
			}
		}
		
		//Handle affix lengths
		if(mode==AFFIX_MODE && names!=null){
			
			{//Find the lengths of all names to fill the lengthArray array
				BitSet bs=new BitSet();//A set of observed affix lengths
				
				//Add the fixed affix length, if present
				if(fixedLength>0){bs.set(fixedLength);}
				
				//Add the lengths of everything in the names set
				for(String s : names.keySet()){bs.set(s.length());}
				
				//Put the set members into an array for faster access
				lengthArray=new int[bs.cardinality()];
				for(int i=0, bit=-1; i<lengthArray.length; i++){
					bit=bs.nextSetBit(bit+1);
					lengthArray[i]=bit;
				}
				
				//Sort and reverse so that the lengths are descending
				Arrays.sort(lengthArray);
				Tools.reverseInPlace(lengthArray);
			}
			
			assert((lengthArray.length>0 && lengthArray[0]>0)) : 
				"Must include at least one name, an affix length, or a delimiter.";
			
			//If there are more than one length, disable the fixedLength mode
			if(lengthArray!=null && lengthArray.length>1){fixedLength=-1;}
			else if(lengthArray.length==1){
				assert(fixedLength<1 || fixedLength==lengthArray[0]) : 
					"\nLength flag ("+fixedLength+") does not match detected name length of "+lengthArray[0]+
					"\nPlease omit the length flag or set it correctly.\n";
				fixedLength=lengthArray[0];
				lengthArray=null;
			}
		}
		
		//For length-1 delimiters, use a char instead of a String
		assert(column<0 || delimiter!=null) : "Column may not be set if there is no delimiter.";
		if(delimiter!=null){
			if(delimiter.length()==1){
				delimiterChar=delimiter.charAt(0);
			}else if(delimiter.length()==2 && delimiter.charAt(0)=='\\'){
				delimiterChar=delimiter.charAt(1);
			}else{
				delimiterChar=0;
			}
			delimiterPattern=Pattern.compile(delimiter);
		}
		
		//Populate names table with mutants for hamming distance
		if(hdist>0 && names.size()>0){
			@SuppressWarnings("unused")
			int numMutants=mutate(names, hdist, outstream);
		}
		
		//Populate name list with keys from the name map
		if(names==null){
			nameList=null;
		}else{
			for(Entry<String, String> e : names.entrySet()){
				nameList.add(e.getKey());
			}
		}
		
		//Ensure there are names in substring mode
		if(mode==SUBSTRING_MODE){
			assert(nameList!=null && !nameList.isEmpty()) : "Empty names list is not allowed in substring mode.";
		}
		
		//Perform symbol replacements in output filenames
		if(remap!=null && names!=null && !names.isEmpty()){
			for(String key : nameList){
				String old=names.get(key);
				String value=Tools.remap(remap, old);
				if(!value.equals(old)){names.put(key, value);}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Primary method.
	 * Starts and stops I/O streams, processes the data, and prints results.
	 * @param t A timer that was already started.
	 */
	void process(Timer t){
		
		//Create stream for input reads
		final ConcurrentReadInputStream cris=makeInputStream();
		
		//Create streams for output reads other than unmatched
		final BufferedMultiCros mcros=makeMatchedOutputStream(cris.paired());
		
		//Create stream for unmatched output reads
		final ConcurrentReadOutputStream rosu=makeUnmatchedOutputStream();
		
		//Streams are set up, so process the reads
		processInner(cris, mcros, rosu);
		//At this point processing has finished.
		
		//Close streams
		cleanup(cris, mcros, rosu);
		
		//Report statistics to file
		if(stats!=null){printReport(mcros);}
		
		//Stop the timer
		t.stop();
		
		//Print results
		printResultsToScreen(t);
		
		if(errorState){//Exit with a message and error code
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Create and start the input stream (cris). */
	private ConcurrentReadInputStream makeInputStream(){
		//Stream for input reads
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
		if(verbose){outstream.println("Started cris");}
		cris.start();
		
		//Report whether the input is being processed as paired.
		final boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}
	
	/** Create and start the normal output streams (mcros). */
	private BufferedMultiCros makeMatchedOutputStream(boolean pairedInput){
		//Streams for output reads other than unmatched
		final BufferedMultiCros mcros;
		if(out1!=null){
			
			if(mcrosType==2){//Slow, synchronous mcros type
				mcros=new MultiCros2(out1, out2, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, threaded);
			}else if(mcrosType==3){//Faster, asynchronous type
				mcros=new MultiCros3(out1, out2, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, threaded, maxStreams);
			}else{
				throw new RuntimeException("Bad mcrosType: "+mcrosType);
			}
			
			//Set mcros fields not available in constructor
			mcros.readsPerBuffer=readsPerBuffer;
			mcros.bytesPerBuffer=bytesPerBuffer;
			mcros.minReadsToDump=minReadsToDump;
			mcros.trackCardinality=trackCardinality;
			
			//Report output format
			if(pairedInput && out2==null && (in1==null || !ffin1.samOrBam())){
				outstream.println("Writing interleaved.");
			}
			
			//Start the mcros thread if necessary
			if(threaded){mcros.start();}
		}else{
			mcros=null;
		}
		return mcros;
	}
		
	/** Create and start the unmatched output stream (rosu). */
	private ConcurrentReadOutputStream makeUnmatchedOutputStream(){
		//Stream for unmatched output reads
		final ConcurrentReadOutputStream rosu;
		if(outu1!=null){
			//Number of output buffers; does not need to be high since access is single-threaded
			final int buff=4;
			
			FileFormat ffout1=FileFormat.testOutput(outu1, FileFormat.FASTQ, extout, true, overwrite, append, false);
			FileFormat ffout2=(outu2==null ? null : FileFormat.testOutput(outu2, FileFormat.FASTQ, extout, true, overwrite, append, false));
			rosu=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, true);
			rosu.start();
		}else{
			rosu=null;
		}
		
		return rosu;
	}
	
	/** 
	 * Process all reads.
	 * @param cris Input stream (required).
	 * @param mcros Matched read output stream (optional).
	 * @param rosu Unmatched read output stream (optional).
	 */
	private void processInner(final ConcurrentReadInputStream cris, final BufferedMultiCros mcros, final ConcurrentReadOutputStream rosu) {
		
		//Fetch the first list
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		//Ensure pairing seems correct
		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
		}
		
		//While there are more reads to process...
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			//List to accumulate unmatched reads
			ArrayList<Read> unmatched=new ArrayList<Read>();
			
			//For each read...
			for(Read r : reads){
				
				//Get the target file identifier
				String name=getValue(r);
				final int pairCount=r.pairCount(), pairLen=r.pairLength();
				if(name!=null){
					//Set the name so that the mcros will send it to the right place
					r.obj=name;
					readsOut+=pairCount;
					basesOut+=pairLen;
				}else{
					//Send it to unmatched
					unmatched.add(r);
					readsUnmatched+=pairCount;
					basesUnmatched+=pairLen;
				}

				readsProcessed+=pairCount;
				basesProcessed+=pairLen;
			}
			if(rosu!=null){rosu.add(unmatched, ln.id);}//Send unmatched reads to rosu
			if(mcros!=null){mcros.add(reads);}//Send matched reads to mcros
			
			//Notify the input stream that the list has been processed
			cris.returnList(ln);
			
			//Fetch a new list
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		
		//Notify the input stream that this thread is done processing
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	/** 
	 * Finish writing and close streams.
	 * @param cris Input stream (required).
	 * @param mcros Matched read output stream (optional).
	 * @param rosu Unmatched read output stream (optional).
	 */
	private void cleanup(final ConcurrentReadInputStream cris, final BufferedMultiCros mcros, final ConcurrentReadOutputStream rosu){
		
		//Shut down mcros
		if(mcros!=null){
			mcros.close();
			errorState|=mcros.errorState();
		}

		if(minReadsToDump>0 && mcros!=null){
			//Dump the residual reads into the unmatched file
			mcros.dumpResidual(rosu);

			//Adjust statistics to reflect that residuals did not get demultiplexed
			readsOut-=mcros.residualReads;
			basesOut-=mcros.residualBases;
		}

		//errorState|=ReadStats.writeAll(); //Currently unused; allows tracking statistics.
		
		//Shut down normal streams
		errorState|=ReadWrite.closeStreams(cris, rosu);
	}
	
	/**
	 * @param t A Timer that has already been stopped
	 */
	private void printResultsToScreen(Timer t) {
		//Report statistics to screen
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		outstream.println("Time:               "+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format(Locale.ROOT, "%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+String.format(Locale.ROOT, "%.2fm bases/sec", bpnano*1000));
		outstream.println("Reads Out:          "+readsOut);
		outstream.println("Bases Out:          "+basesOut);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Generate mutants for all names according to the hamming distance
	 * @param names Map of names to mutate
	 * @param hdist Hamming distance
	 * @param outstream Report results here
	 * @return Number of mutants added (net)
	 */
	private static int mutate(HashMap<String, String> names, int hdist, PrintStream outstream){
		
		//Put the names in a list for iteration, to prevent concurrent modification of the table
		ArrayList<String> list=new ArrayList<String>(names.size());
		list.addAll(names.values());
		
		//Number of mutants generated
		int mutants=0;
		
		//Set of colliding keys to be removed
		final HashSet<String> collisions=new HashSet<String>();
		
		//Generate the mutants
		for(String name : list){
			mutants+=addMutants(name.getBytes(), name, hdist, names, collisions);
		}
		outstream.println("Added "+mutants+" mutants to the initial "+list.size()+" names.");
		
		//Remove original keys from collisions so that it only contains mutants
		collisions.removeAll(list);
		
		int removed=0;
		if(!collisions.isEmpty()){//Remove colliding mutants from the name table
			for(String key : collisions) {
				names.remove(key);
				removed++;
			}
			outstream.println("Removed "+removed+" collisions due to ambiguity.");
		}
		
		final int netChange=mutants-removed;
		return netChange;
	}
	
	/** 
	 * Generate mutants for one name, recursively
	 * @param keyArray Current key (name) as a byte array, to be modified
	 * @param value Name to which this key should map
	 * @param hdist Remaining Hamming distance
	 * @param names Names table to insert the mutants
	 * @param collisions Set of keys mapping to multiple values
	 * @return Number of mutants added
	 */
	private static int addMutants(final byte[] keyArray, final String value, final int hdist, final HashMap<String, String> names, final HashSet<String> collisions){
		assert(hdist>0);
		int added=0;
		
		//For each position in the array, change it to all possible values
		for(int i=0; i<keyArray.length; i++){
			final byte old=keyArray[i];
			if(AminoAcid.isACGTN(old)){//Only mutate bases
				for(byte b : symbols){
					if(b!=old){
						keyArray[i]=b;
						{
							String key=new String(keyArray);
							String oldValue=names.get(key);
							if(oldValue==null){
								names.put(key, value);
								added++;
//								outstream.println("Added "+key+"->"+value+"; oldValue="+oldValue);
							}else if(!oldValue.equals(value)){
//								assert(value.equals(names.get(s))) : "Collision between "+value+" and "+names.get(s)+" for mutant "+s;
								collisions.add(key);
//								outstream.println("Collision for "+key+"->"+value+"; oldValue="+oldValue);
							}
						}
						if(hdist>1){
							//Recur on the mutant, with a reduced Hamming distance
							added+=addMutants(keyArray, value, hdist-1, names, collisions);
						}
					}
				}
				//Restore the original value before moving to the next position
				keyArray[i]=old;
			}
		}
		return added;
	}
	
	/** Returns the value this read maps to; meaning, the variable part of the output filename */
	private String getValue(Read r){
		final String key=getKey(r);
		if(verbose){System.err.println("Got key "+key);}
		if(key==null){return null;}
		
		//Return the value this key maps to
		if(names!=null){
			if(verbose){System.err.println("Got value "+names.get(key));}
			return names.get(key);
		}
		
		//Return the value generated by remapping the key's symbols, if needed
		final String value=(remap==null ? key : Tools.remap(remap,  key));
		if(verbose){System.err.println("Returning value "+value);}
		return value;
	}
	
	/** Generates a key from the read header */
	private String getKey(Read r){
		final String id=r.id;
		final int idlen=id.length();
		final String key;
		
		//Simplest case; every header is its own key
		if(mode==HEADER_MODE){return id;}
		
		//Parse the barcode
		if(mode==BARCODE_MODE){return r.barcode(true);}
		
		//Use mapping information
		if(mode==CHROM_MODE){return r.rnameS();}
		
		//Return a prefix or suffix
		if(mode==AFFIX_MODE){
			if(fixedLength>0){//If there is a single fixed affix length, return the affix
				key=(id.length()<=fixedLength ? id : prefixMode ? id.substring(0, fixedLength) : id.substring(idlen-fixedLength));
				return key;
			}else{//Make a substring of each possible length and look it up in the hashmap
				for(int affixLen : lengthArray){
					final String sub=idlen>=affixLen ? prefixMode ? id.substring(0, affixLen) : id.substring(idlen-affixLen) : id;
					if(names.containsKey(sub)){return sub;}
				}
				return null;
			}
		}
		
		//Return a delimited substring
		if(mode==DELIMITER_MODE){
			
			if(column>-1){//Use a split operation, which is slow
				
//				String[] split=id.split(delimiter);
				String[] split=delimiterPattern.split(id);//Faster; uses a precompiled pattern
				assert(split.length>1) : "Delimiter '"+delimiter+"' was not found in name '"+id+"'";
				
				int col=Tools.min(column, split.length-1);
				key=split[col];
				if(col!=column && !warned){
					outstream.println("*** WARNING! ***\n"
							+ "Only "+(col+1)+" columns for record "+id+"\n"
							+ "Further warnings will be suppressed.\n");
					warned=true;
					assert(errorState=true); //Change error state to true if assertions are enabled.
				}
			}else if(prefixMode){//Use a delimited prefix, which is fast
				//Find the first index of the delimiter, using the char when possible
				int idx=(delimiterChar>0 ? id.indexOf(delimiterChar) : id.indexOf(delimiter));
				assert(idx>=0) : "Delimiter\r\n" + 
						"		System.err.println(\"a\"); '"+delimiter+"' was not found in name '"+id+"'";
				key=id.substring(0, idx);
			}else{//Use a delimited suffix, which is fast
				//Find the last index of the delimiter, using the char when possible
				int idx=(delimiterChar>0 ? id.lastIndexOf(delimiterChar) : id.lastIndexOf(delimiter));
				assert(idx>=0) : "Delimiter '"+delimiter+"' was not found in name '"+id+"'";
				key=id.substring(idx+delimiter.length());
			}
			return key;
		}
		
		//Substring mode requires brute-force matching and is very slow
		if(mode==SUBSTRING_MODE){
			if(nameList.size()>0){
				for(String s : nameList){
					if(id.contains(s)){return s;}
				}
				return null;
			}
		}
		
		throw new RuntimeException("No mode is being used: "+mode); //Should be unreachable
	}
	
	/** 
	 * Print statistics about demultiplexing to a file.
	 * @param mcros Output stream, after processing is completely finished.
	 */
	void printReport(BufferedMultiCros mcros){
		if(stats==null){return;}
		
		//Make a writer for the stats file
		ByteStreamWriter bsw=new ByteStreamWriter(stats, overwrite, append, true);
		bsw.start();
		
		{
			ByteBuilder bb=new ByteBuilder();

			//Print the header
			bb.append("#ReadsIn\t").append(readsProcessed).nl();
			bb.append("#BasesIn\t").append(basesProcessed).nl();
			bb.append("#ReadsOut\t").append(readsOut).nl();
			bb.append("#BasesOut\t").append(basesOut).nl();
			bb.append("#Name\tReads\tBases"+(trackCardinality ? "\tCardinality" : "")+"\n");
			
			//Print results of unmatched reads, which is not captured by the mcros
			if(names!=null && !names.isEmpty()){
				bb.append("Unmatched\t").append(readsUnmatched).tab().append(basesUnmatched).nl();
			}
			bsw.print(bb);
		}
		
		if(mcros!=null){//Print results from the mcros
			ByteBuilder bb=mcros.report();
			bsw.print(bb);
		}
		
		//Finish writing
		bsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file */
	private String in1=null;
	/** Optional input file for read 2 of paired reads in twin files */
	private String in2=null;
	
	/** Input qual file 1 */
	private String qfin1=null;
	/** Input qual file 2 */
	private String qfin2=null;

	/** Output file pattern 1 */
	private String out1=null;
	/** Output file pattern 2 */
	private String out2=null;

	/** Output qual file 1 */
	private String qfout1=null;
	/** Output qual file 2 */
	private String qfout2=null;

	/** Unmatched read output file 1 */
	private String outu1=null;
	/** Unmatched read output file 2 */
	private String outu2=null;
	
	/** File extension override for input files */
	private String extin=null;
	/** File extension override for output files */
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	/** Primary input file */
	private final FileFormat ffin1;
	
	/** Read 2 input file */
	private final FileFormat ffin2;
	
	/*--------------------------------------------------------------*/
	
	/** Input reads */
	long readsProcessed=0;
	long basesProcessed=0;
	
	/** Demultiplexed output reads */
	long readsOut=0;
	long basesOut=0;
	
	/** Output reads that did not get demultiplexed */
	long readsUnmatched=0;
	long basesUnmatched=0;

	/** Stop after this many input reads */
	private long maxReads=-1;

	/** File to print number of reads sent to each output file */
	private String stats=null;

	/** For splitting headers on a symbol */
	private String delimiter=null;
	
	/** For splitting headers faster if the delimiter is just one character */
	private char delimiterChar=0;
	
	/** Precompiled pattern matching the delimiter */
	private Pattern delimiterPattern=null;
	
	/** For replacing symbols in the filename */
	private byte[] remap=null;
	
	/** If there is a delimiter, use this column after splitting.
	 * Column is 1-based. */
	private int column=-1;
	
	/** Use the prefix of a header.  If false, use the suffix. 
	 * This is ignored if column is set. */
	private boolean prefixMode=true;
	
	/** How to select a read's key; see static modes array */
	private int mode=AFFIX_MODE;
	
//	/** Assume Illumina header format, and used the barcode. */
//	private boolean mode==BARCODE_MODE=false;
//	
//	/** Allow names to be any substring of a header */
//	private boolean substringMode=false;
//	
//	/** Demultiplex every sequence into its own file */
//	private boolean mode==HEADER_MODE=false;
//
//	/** Demultiplex mapped reads into one file per reference sequence */
//	private boolean mode==CHROM_MODE=false;

	/** Track per-output-file cardinality; may be slow */
	private boolean trackCardinality=false;
	
	/** Prevents issuing warnings multiple times */
	private boolean warned=false;

	/** Mcros buffers will be dumped after accumulating this many reads */
	private int readsPerBuffer=1600;

	/** Mcros buffers will be dumped after accumulating this many reads */
	private int bytesPerBuffer=3000000;
	
	/** Hamming distance for read indexes */
	private int hdist=0;
	
	/** 
	 * Do not create files with under this many reads.
	 * Increases memory usage since reads must be retained until processing is finished.
	 * Upon running out of memory the program may go very slowly.
	 */
	private long minReadsToDump=0;
	
	/** Affix length if fixed */
	private int fixedLength=-1;
	
	/** All possible affix lengths if there are multiple.
	 * Basically, the lengths of everything in names. */
	private int[] lengthArray;
	
	/** 
	 * A set of recognized names and the output names they map to.
	 * The keys and values will be identical unless a Hamming distance is used.
	 * This is only filled if explicit names are provided.
	 */
	private HashMap<String, String> names=new HashMap<String, String>();
	
	/** Names in list form, for substring matching */
	private ArrayList<String> nameList=new ArrayList<String>();
	
	/** 
	 * Select MultiCros2 or MultiCros3.
	 * MultiCros2 only opens 1 stream at a time and does not do concurrent writes.
	 * MultiCros3 maintains multiple streams and does concurrent writes even with only 1 stream.
	 */
	private int mcrosType=3;
	
	/** Max concurrent streams for MultiCros3.  Does not seem to speed things up. */
	private int maxStreams=BufferedMultiCros.DEFAULT_MAX_STREAMS;
	
	/** 
	 * This puts the mcros handling in a second thread.
	 * Does not seem to speed things up.
	 */
	private boolean threaded=false;
	
	/** Whether input file interleaving was explicitly set */
	private boolean setInterleaved=false;
	
	/*--------------------------------------------------------------*/
	
	/** Print messages here */
	private PrintStream outstream=System.err;
	
	/** True if errors were encountered */
	public boolean errorState=false;
	
	/** Permission to overwrite existing files. */
	private boolean overwrite=true;
	
	/** 
	 * Append to existing files rather than overwriting.
	 * It is not advisable to set this flag for this class. 
	 * */
	private boolean append=false;
	
	/** Retain header of sam/bam files. */
	private boolean useSharedHeader=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Operation mode for selection of key */
	public static final int AFFIX_MODE=1, DELIMITER_MODE=2, BARCODE_MODE=3, SUBSTRING_MODE=4, HEADER_MODE=5, CHROM_MODE=6;
	
	/** 
	 * Symbols allowed to substitute for Hamming distance.
	 * Hamming distance is intended for barcodes only.
	 */
	private static final byte[] symbols={'A', 'C', 'G', 'T', 'N'};
	
	/** Verbose messages for debugging */
	public static boolean verbose=false;
	
}
