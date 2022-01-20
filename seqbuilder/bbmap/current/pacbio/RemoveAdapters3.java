package pacbio;

import java.util.ArrayList;
import java.util.Locale;

import align2.MultiStateAligner9PacBio;
import aligner.MultiStateAligner9PacBioAdapter;
import dna.AminoAcid;
import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Increased sensitivity to nearby adapters.
 * Does reverse-complementation on low-scoring suspects.
 * @author Brian Bushnell
 * @date Nov 5, 2012
 *
 */
public class RemoveAdapters3 {

	public static void main(String[] args){
		
		boolean verbose=false;
		int ziplevel=-1;

		String in1=null;
		String in2=null;
		long maxReads=-1;
		
		String outname1=null;
		String outname2=null;
		
		String query=pacbioAdapter;
		
		boolean splitReads=false;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(a.equals("path") || a.equals("root") || a.equals("tempdir")){
				Data.setPath(b);
			}else if(a.equals("fasta") || a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
				if(b.indexOf('#')>-1){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("query") || a.equals("adapter")){
				query=b;
			}else if(a.equals("split")){
				splitReads=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
				System.out.println("Set overwrite to "+overwrite);
			}else if(a.equals("threads") || a.equals("t")){
				if(b.equalsIgnoreCase("auto")){THREADS=Shared.LOGICAL_PROCESSORS;}
				else{THREADS=Integer.parseInt(b);}
				System.out.println("Set threads to "+THREADS);
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.startsWith("outname") || a.startsWith("outfile") || a.equals("out")){
				if(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none") || split.length==1){
					System.out.println("No output file.");
					outname1=null;
					OUTPUT_READS=false;
				}else{
					OUTPUT_READS=true;
					if(b.indexOf('#')>-1){
						outname1=b.replace('#', '1');
						outname2=b.replace('#', '2');
					}else{
						outname1=b;
					}
				}
			}else if(a.equals("perfectmode")){
				PERFECTMODE=Parse.parseBoolean(b);
				if(ziplevel==-1){ziplevel=2;}
			}else if(a.equals("minratio")){
				MINIMUM_ALIGNMENT_SCORE_RATIO=Float.parseFloat(b);
				System.out.println("Set MINIMUM_ALIGNMENT_SCORE_RATIO to "+MINIMUM_ALIGNMENT_SCORE_RATIO);
			}else if(a.startsWith("verbose")){
				verbose=Parse.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter: "+args[i]);
			}
		}
		
		Parser.processQuality();
		
		assert(FastaReadInputStream.settingsOK());
		if(in1==null){throw new RuntimeException("Please specify input file.");}
		
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
//			if(verbose){System.err.println("Started cris");}
//			cris.start(); //4567
//			th.start();
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		ConcurrentReadOutputStream ros=null;
		if(OUTPUT_READS){
			final int buff=(!ordered ? THREADS : Tools.max(24, 2*THREADS));
			
			FileFormat ff1=FileFormat.testOutput(outname1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(outname2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			ros=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
		}
		process(cris, ros, query, splitReads);
	}
	
	public static void process(ConcurrentReadInputStream cris, ConcurrentReadOutputStream ros, String query, boolean split){

		Timer t=new Timer();
		
		cris.start(); //4567
		
		System.out.println("Started read stream.");
		

		if(ros!=null){
			ros.start();
			System.out.println("Started output threads.");
		}
		ProcessThread[] pts=new ProcessThread[THREADS];
		for(int i=0; i<pts.length; i++){
			pts[i]=new ProcessThread(cris, ros, MINIMUM_ALIGNMENT_SCORE_RATIO, query, split);
			pts[i].start();
		}
		System.out.println("Started "+pts.length+" processing thread"+(pts.length==1 ? "" : "s")+".");
		
		for(int i=0; i<pts.length; i++){
			ProcessThread pt=pts[i];
			synchronized(pt){
				while(pt.getState()!=Thread.State.TERMINATED){
					try {
						pt.join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			if(i==0){
				System.out.print("Detecting finished threads: 0");
			}else{
				System.out.print(", "+i);
			}
		}
		System.out.println('\n');
		ReadWrite.closeStreams(cris, ros);
		
		printStatistics(pts);
		
		t.stop();
		System.out.println("Time: \t"+t);
		
	}
	
	public static void printStatistics(ProcessThread[] pts){

		long plusAdaptersFound=0;
		long minusAdaptersFound=0;
		long goodReadsFound=0;
		long badReadsFound=0;
		
		long truepositive=0;
		long truenegative=0;
		long falsepositive=0;
		long falsenegative=0;
		long expected=0;
		long unexpected=0;
		
		for(ProcessThread pt : pts){
			plusAdaptersFound+=pt.plusAdaptersFound;
			minusAdaptersFound+=pt.minusAdaptersFound;
			goodReadsFound+=pt.goodReadsFound;
			badReadsFound+=pt.badReadsFound;

			truepositive+=pt.truepositive;
			truenegative+=pt.truenegative;
			falsepositive+=pt.falsepositive;
			falsenegative+=pt.falsenegative;
			expected+=pt.expected;
			unexpected+=pt.unexpected;
		}
		
		long totalReads=goodReadsFound+badReadsFound;
		if(expected<1){expected=1;}
		if(unexpected<1){unexpected=1;}

		System.out.println("Good reads:              \t"+goodReadsFound);
		System.out.println("Bad reads:               \t"+badReadsFound);
		System.out.println("Plus adapters:           \t"+plusAdaptersFound);
		System.out.println("Minus adapters:          \t"+minusAdaptersFound);
		System.out.println();
		if(truepositive>0 || truenegative>0 || falsepositive>0 || falsenegative>0){
			System.out.println("Adapters Expected:       \t"+expected);
			System.out.println("True Positive:           \t"+truepositive+" \t"+String.format(Locale.ROOT, "%.3f%%", truepositive*100f/expected));
			System.out.println("True Negative:           \t"+truenegative+" \t"+String.format(Locale.ROOT, "%.3f%%", truenegative*100f/unexpected));
			System.out.println("False Positive:          \t"+falsepositive+" \t"+String.format(Locale.ROOT, "%.3f%%", falsepositive*100f/unexpected));
			System.out.println("False Negative:          \t"+falsenegative+" \t"+String.format(Locale.ROOT, "%.3f%%", falsenegative*100f/expected));
		}
		
	}
	
	private static class ProcessThread extends Thread{
		
		public ProcessThread(ConcurrentReadInputStream cris_,
				ConcurrentReadOutputStream ros_, float minRatio_, String query_, boolean split_) {
			cris=cris_;
			ros=ros_;
			minRatio=minRatio_;
			query1=query_.getBytes();
			query2=AminoAcid.reverseComplementBases(query1);
			
			stride=(int)(query1.length*0.95f);
			window=(int)(query1.length*2.5f+10);
			
			ALIGN_ROWS=Tools.max(query1.length, rcompDistance)+1;
			ALIGN_COLUMNS=Tools.max(window, rcompDistance)+5;
			SPLIT=split_;
			
			assert(window<ALIGN_COLUMNS);

			msa=new MultiStateAligner9PacBioAdapter(ALIGN_ROWS, ALIGN_COLUMNS);
			msaR=new MultiStateAligner9PacBio(ALIGN_ROWS, ALIGN_COLUMNS); //TODO: Why is this not 'adapter' version?

			maxSwScore=msa.maxQuality(query1.length);
			minSwScore=(int)(maxSwScore*MINIMUM_ALIGNMENT_SCORE_RATIO);
			minSwScoreRcomp=(int)(msa.maxQuality(rcompDistance)*MINIMUM_ALIGNMENT_SCORE_RATIO);
			minSwScoreSuspect=(int)(maxSwScore*Tools.min(MINIMUM_ALIGNMENT_SCORE_RATIO*0.7f, MINIMUM_ALIGNMENT_SCORE_RATIO-0.05f));
			maxImperfectSwScore=msa.maxImperfectScore(query1.length);

			rcomp1=new byte[rcompDistance];
			rcomp2=new byte[rcompDistance];
		}
		
		@Override
		public void run(){
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> readlist=ln.list;
			
			while(!readlist.isEmpty()){
				
				//System.err.println("Got a list of size "+readlist.size());
				for(int i=0; i<readlist.size(); i++){
					Read r=readlist.get(i);

//					System.out.println("Got read: "+r.toText());
//					System.out.println("Synthetic: "+r.synthetic());
					
					if(r.synthetic()){syntheticReads++;}
					
					processRead(r);
					if(r.mate!=null){processRead(r.mate);}
					
				}
				
//				System.err.println("outputStream = "+outputStream==null ? "null" : "real");
				if(ros!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
					if(DONT_OUTPUT_BROKEN_READS){removeDiscarded(readlist);}
					for(Read r : readlist){
						if(r!=null){
							r.obj=null;//Not sure what r.obj is here
							assert(r.bases!=null);
							if(r.sites!=null && r.sites.isEmpty()){r.sites=null;}
						}
					}
//					System.err.println("Adding list of length "+readlist.size());
					if(SPLIT){
						ros.add(split(readlist), ln.id);
					}else{
						ros.add(readlist, ln.id);
					}
				}
				
				cris.returnList(ln.id, readlist.isEmpty());
				
				//System.err.println("Waiting on a list...");
				ln=cris.nextList();
				readlist=ln.list;
			}
			
			//System.err.println("Returning a list... (final)");
			assert(readlist.isEmpty());
			cris.returnList(ln.id, readlist.isEmpty());
		}

		/**
		 * @param readlist
		 * @return
		 */
		private static ArrayList<Read> split(ArrayList<Read> in) {
			ArrayList<Read> out=new ArrayList<Read>(in.size());
			for(Read r : in){
				if(r!=null){
					assert(r.mate==null);
					if(!r.hasAdapter()){out.add(r);}
					else{out.addAll(split(r));}
				}
			}
			return out;
		}

		/**
		 * @param r
		 * @return
		 */
		private static ArrayList<Read> split(Read r) {
			ArrayList<Read> sections=new ArrayList<Read>();
			
			int lastX=-1;
			for(int i=0; i<r.length(); i++){
				if(r.bases[i]=='X'){
					if(i-lastX>minContig){
						byte[] b=KillSwitch.copyOfRange(r.bases, lastX+1, i);
						byte[] q=r.quality==null ? null : KillSwitch.copyOfRange(r.quality, lastX+1, i);
						Read r2=new Read(b, q, r.id+"_"+(lastX+1), r.numericID);
						sections.add(r2);
					}
					lastX=i;
				}
			}
			int i=r.length();
			if(i-lastX>minContig){
				byte[] b=KillSwitch.copyOfRange(r.bases, lastX+1, i);
				byte[] q=r.quality==null ? null : KillSwitch.copyOfRange(r.quality, lastX+1, i);
				Read r2=new Read(b, q, r.id+"_"+(lastX+1), r.numericID);
				sections.add(r2);
			}
			return sections;
		}

		/**
		 * @param r
		 */
		private int processRead(Read r) {
			
			int begin=0;
			while(begin<r.length() && r.bases[begin]=='N'){begin++;} //Skip reads made of 'N'
			if(begin>=r.length()){return 0;}
			
			final byte[] array=npad(r.bases, npad);
			
			int lim=array.length-npad-stride;
			
			int plusFound=0;
			int minusFound=0;
			
			int lastSuspect=-1;
			
			for(int i=begin; i<lim; i+=stride){
				int j=Tools.min(i+window, array.length-1);
				if(j-i<window){
					lim=0; //Last loop cycle
//					i=Tools.max(0, array.length-2*query1.length);
				}
				
				{
					int[] rvec=msa.fillAndScoreLimited(query2, array, i, j, minSwScoreSuspect);
					if(rvec!=null && rvec[0]>=minSwScoreSuspect){
						int score=rvec[0];
						int start=rvec[1];
						int stop=rvec[2];
						assert(score>=minSwScoreSuspect);
						if((i==0 || start>i) && (j==array.length-1 || stop<j)){
							boolean kill=(score>=minSwScore || (lastSuspect>0 && start>=lastSuspect && start-lastSuspect<suspectDistance));
							
							if(!kill && array.length-stop>window){//Look ahead
								rvec=msa.fillAndScoreLimited(query2, array, stop, stop+window, minSwScoreSuspect);
								if(rvec!=null && rvec[0]>=minSwScoreSuspect && rvec[1]-stop<suspectDistance){kill=true;}
							}
							
							if(RCOMP && !kill && start>rcompDistance && array.length-stop>rcompDistance+1){
								kill=testRcomp(array, start, stop);
//								System.out.print(kill ? "#" : ".");
							}
							
							if(kill){
//								System.out.println("-:"+score+", "+minSwScore+", "+minSwScoreSuspect+"\t"+lastSuspect+", "+start+", "+stop);
								minusFound++;
								for(int x=Tools.max(0, start); x<=stop; x++){array[x]='X';}
							}
						}
//						System.out.println("Set lastSuspect="+stop+" on score "+score);
						lastSuspect=stop;
					}
				}
				
				{
					int[] rvec=msa.fillAndScoreLimited(query1, array, i, j, minSwScoreSuspect);
					if(rvec!=null && rvec[0]>=minSwScoreSuspect){
						int score=rvec[0];
						int start=rvec[1];
						int stop=rvec[2];
						if((i==0 || start>i) && (j==array.length-1 || stop<j)){
							boolean kill=(score>=minSwScore || (lastSuspect>0 && start>=lastSuspect && start-lastSuspect<suspectDistance));
							
							if(!kill && array.length-stop>window){//Look ahead
								rvec=msa.fillAndScoreLimited(query1, array, stop, stop+window, minSwScoreSuspect);
								if(rvec!=null && rvec[0]>=minSwScoreSuspect && rvec[1]-stop<suspectDistance){kill=true;}
							}
							
							if(RCOMP && !kill && start>rcompDistance && array.length-stop>rcompDistance+1){
								kill=testRcomp(array, start, stop);
//								System.out.print(kill ? "#" : ".");
							}
							
							if(kill){
//								System.out.println("+:"+score+", "+minSwScore+", "+minSwScoreSuspect+"\t"+lastSuspect+", "+start+", "+stop);
								plusFound++;
								for(int x=Tools.max(0, start); x<=stop; x++){array[x]='X';}
							}
						}
//						System.out.println("Set lastSuspect="+stop+" on score "+score);
						lastSuspect=stop;
					}
				}
			}
			
			if(r.synthetic()){
				if(/*r.hasadapter() && */(r.numericID&3)==0){
					if(plusFound>0){truepositive++;}else{falsenegative++;}
					if(plusFound>1){falsepositive+=(plusFound-1);}
					falsepositive+=minusFound;
					expected++;
				}else if(/*r.hasadapter() && */(r.numericID&3)==1){
					if(minusFound>0){truepositive++;}else{falsenegative++;}
					if(minusFound>1){falsepositive+=(minusFound-1);}
					falsepositive+=plusFound;
					expected++;
				}else{
					falsepositive=falsepositive+plusFound+minusFound;
					if(plusFound+minusFound==0){truenegative++;}
					unexpected++;
				}
			}
			
			plusAdaptersFound+=plusFound;
			minusAdaptersFound+=minusFound;
			int found=plusFound+minusFound;
			if(found>0){
				for(int i=npad, j=0; j<r.length(); i++, j++){r.bases[j]=array[i];}
				if(DONT_OUTPUT_BROKEN_READS){r.setDiscarded(true);}
				badReadsFound++;
			}else{
				goodReadsFound++;
			}
			
			r.setHasAdapter(found>0);
			
			return found;
			
		}
		
		/**
		 * @param array
		 * @param start
		 * @param stop
		 * @return
		 */
		private boolean testRcomp(byte[] array, int start, int stop) {
			
			for(int i=0, j=start-rcompDistance, k=stop+1; i<rcompDistance; i++, j++, k++){
				rcomp1[i]=array[j];
				if(rcomp1[i]=='X'){rcomp1[i]='N';}
				rcomp2[i]=array[k];
				if(rcomp2[i]=='X'){rcomp2[i]='N';}
			}
			AminoAcid.reverseComplementBasesInPlace(rcomp2);
			
//			System.out.println(new String(array).substring(start-rcompDistance, stop+rcompDistance));
//			System.out.println(new String(rcomp1));
//			System.out.println(new String(rcomp2));

//			int[] rvec=msa.fillAndScoreLimited(rcomp1, rcomp2, 0, rcompDistance-1, minSwScoreSuspect);
			int[] rvec=msaR.fillAndScoreLimited(rcomp1, rcomp2, 0, rcompDistance-1, 0, null);
			int score=rvec==null ? -99999 : rvec[0];
//			System.out.println(minSwScoreSuspect+", "+score+"\n");
			return score>=1000;
//			return score>=minSwScoreRcomp;
		}

		private byte[] npad(final byte[] array, final int pad){
			final int len=array.length+2*pad;
			if(padbuffer==null || padbuffer.length!=len){padbuffer=new byte[len];}
			byte[] r=padbuffer;
			for(int i=0; i<pad; i++){r[i]='N';}
			for(int i=pad, j=0; j<array.length; i++, j++){r[i]=array[j];}
			for(int i=array.length+pad, limit=Tools.min(r.length, len+50); i<limit; i++){r[i]='N';}
			padbuffer=null; //Kills the buffer.  Causes more memory allocation, but better cache/NUMA locality if threads move around.
			return r;
		}
		
		private byte[] padbuffer=null;
		private final byte[] rcomp1, rcomp2;
		private final byte[] query1, query2;
		private final ConcurrentReadInputStream cris;
		private final ConcurrentReadOutputStream ros;
		private final float minRatio;
		private final MultiStateAligner9PacBioAdapter msa;
		private final MultiStateAligner9PacBio msaR;
		private final int ALIGN_ROWS;
		private final int ALIGN_COLUMNS;
		private final int stride;
		private final int window;
		private final boolean SPLIT;

		long plusAdaptersFound=0;
		long minusAdaptersFound=0;
		long goodReadsFound=0;
		long badReadsFound=0;
		long truepositive=0;
		long truenegative=0;
		long falsepositive=0;
		long falsenegative=0;
		long expected=0;
		long unexpected=0;
		
		private final int maxSwScore;
		private final int minSwScore;
		private final int minSwScoreSuspect;
		private final int minSwScoreRcomp;
		private final int maxImperfectSwScore;
		
		long syntheticReads=0;
		
	}
	
	private static int removeDiscarded(ArrayList<Read> list){
		int removed=0;
		for(int i=0; i<list.size(); i++){
			Read r=list.get(i);
			if(r.discarded()){
				if(r.mate==null || r.mate.discarded()){
					list.set(i, null);
					removed++;
				}
			}
		}
		return removed;
	}

	public static boolean DONT_OUTPUT_BROKEN_READS;
	/** Permission to overwrite existing files */
	private static boolean overwrite=true;
	/** Permission to append to existing files */
	private static boolean append=false;
	private static int THREADS=Shared.LOGICAL_PROCESSORS;
	private static boolean OUTPUT_READS=false;
	private static boolean ordered=false;
	private static boolean PERFECTMODE=false;
	private static float MINIMUM_ALIGNMENT_SCORE_RATIO=0.31f; //0.31f: At 250bp reads, approx 0.01% false-positive and 94% true-positive.
	public static boolean RCOMP=true;
	private static int npad=35;
	public static int minContig=20;
	public static int suspectDistance=100;
	public static int rcompDistance=80;
	
	public static final String pacbioAdapter="ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT";
	public static final String pacbioStandard_v1="TCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAGAAGGCTGGGCAGGCTATGCACCCTGGTCCAGGTCAAA" +
			"AGCTGCGGAACCCGCTAGCGGCCATCTTGGCCACTAGGGGTCCCGCAGATTCATATTGTCGTCTAGCATGCACAATGCTGCAAACCCAGCTTGCAATGCCCACAGCA" +
			"AGCGGCCAATCTTTACGCCACGTTGAATTGTTTATTACCTGTGACTGGCTATGGCTTGCAACGCCACTCGTAAAACTAGTACTTTGCGGTTAGGGGAAGTAGACAAA" +
			"CCCATTACTCCACTTCCCGGAAGTTCAACTCATTCCAACACGAAATAAAAGTAAACTCAACACCCCAAGCAGGCTATGTGGGGGGGTGATAGGGGTGGATTCTATTT" +
			"CCTATCCCATCCCCTAGGATCTCAATTAAGTTACTAGCGAGTTAAATGTCTGTAGCGATCCCGTCAGTCCTATCGCGCGCATCAAGACCTGGTTGGTTGAGCGTGCA" +
			"GTAGATCATCGATAAGCTGCGAGTTAGGTCATCCCAGACCGCATCTGGCGCCTAAACGTTCAGTGGTAGCTAAGGCGTCACCTTCGACTGTCTAAAGGCAATATGTC" +
			"GTCCTTAGCTCCAAGTCCCTAGCAAGCGTGTCGGGTCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGACCCGACACGCTTGCTAGGGACTTGGAGCT" +
			"AAGGACGACATATTGCCTTTAGACAGTCGAAGGTGACGCCTTAGCTACCACTGAACGTTTAGGCGCCAGATGCGGTCTGGGATGACCTAACTCGCAGCTTATCGATG" +
			"ATCTACTGCACGCTCAACCAACCAGGTCTTGATGCGCGCGATAGGACTGACGGGATCGCTACAGACATTTAACTCGCTAGTAACTTAATTGAGATCCTAGGGGATGG" +
			"GATAGGAAATAGAATCCACCCCTATCACCCCCCCACATAGCCTGCTTGGGGTGTTGAGTTTACTTTTATTTCGTGTTGGAATGAGTTGAACTTCCGGGAAGTGGAGT" +
			"AATGGGTTTGTCTACTTCCCCTAACCGCAAAGTACTAGTTTTACGAGTGGCGTTGCAAGCCATAGCCAGTCACAGGTAATAAACAATTCAACGTGGCGTAAAGATTG" +
			"GCCGCTTGCTGTGGGCATTGCAAGCTGGGTTTGCAGCATTGTGCATGCTAGACGACAATATGAATCTGCGGGACCCCTAGTGGCCAAGATGGCCGCTAGCGGGTTCC" +
			"GCAGCTTTTGACCTGGACCAGGGTGCATAGCCTGCCCAGCCTTCTCTCTCTCTTT";

	
}
