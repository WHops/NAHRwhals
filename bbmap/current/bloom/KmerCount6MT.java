package bloom;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.FileFormat;
import shared.Parse;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Jul 5, 2012
 *
 */
public class KmerCount6MT extends KmerCountAbstract {
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		
		String fname1=args[0];
		String fname2=(args.length>1 ? args[1] : null);
		int k=14;
		int cbits=16;
		int gap=0;
		int matrixbits=-1;
		int hashes=1;
		
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.startsWith("gap")){
				gap=Integer.parseInt(b);
			}else if(a.startsWith("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.startsWith("matrixbits")){
				matrixbits=Integer.parseInt(b);
			}else if(a.startsWith("hashes")){
				hashes=Integer.parseInt(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		int kbits=2*k;
		if(matrixbits<0){
			matrixbits=kbits;
		}
		matrixbits=Tools.min(kbits, matrixbits);
		
		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			assert(!FastaReadInputStream.SPLIT_READS);
			FastaReadInputStream.MIN_READ_LEN=k;
		}
		
		KCountArray count=KCountArray.makeNew(1L<<kbits, 1L<<matrixbits, cbits, gap, hashes);
		count=count(fname1, fname2, k, cbits, gap, true, count);
		count.shutdown();
		
//		verbose=true;
		
		t.stop();
		System.out.println("Finished counting; time = "+t);
		
		printStatistics(count);
		
	}

	public static void printStatistics(KCountArray count){
		long[] freq=count.transformToFrequency();

//		System.out.println(count+"\n");
//		System.out.println(Arrays.toString(freq)+"\n");
		
		long sum=sum(freq);
		System.out.println("Kmer fraction:");
		int lim1=8, lim2=16;
		for(int i=0; i<lim1; i++){
			String prefix=i+"";
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+String.format(Locale.ROOT, "%.3f%%   ",(100l*freq[i]/(double)sum))+"\t"+freq[i]);
		}
		while(lim1<=freq.length){
			int x=0;
			for(int i=lim1; i<lim2; i++){
				x+=freq[i];
			}
			String prefix=lim1+"-"+(lim2-1);
			if(lim2>=freq.length){prefix=lim1+"+";}
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+String.format(Locale.ROOT, "%.3f%%   ",(100l*x/(double)sum))+"\t"+x);
			lim1*=2;
			lim2=min(lim2*2, freq.length);
		}
		
		long sum2=sum-freq[0];
		long x=freq[1];
		System.out.println();
		System.out.println("Keys Counted:  \t         \t"+keysCounted);
		System.out.println("Unique:        \t         \t"+sum2);
		System.out.println("Avg Sites/Key: \t         \t"+String.format(Locale.ROOT, "%.3f    ",(keysCounted*1d/sum2)));
		System.out.println();
		System.out.println("Singleton:     \t"+String.format(Locale.ROOT, "%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
		x=sum2-x;
		System.out.println("Useful:        \t"+String.format(Locale.ROOT, "%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles, int k, int cbits){
		return makeKca(fname1, fname2, extraFiles, k, cbits, 0, Tools.min(2*k, 35), 1, minQuality, true, maxReads, 1, 1, 1, 2);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int k, int cbits, int gap, int matrixbits, int hashes, int minqual, boolean rcomp, long maxreads){
		return makeKca(fname1, fname2, extraFiles, k, cbits, gap, matrixbits, hashes, minqual, rcomp, maxreads, 1, 1, 1, 2);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int k, int cbits, int gap, int matrixbits, int hashes, int minqual, boolean rcomp, long maxreads, int passes, int stepsize, int thresh1, int thresh2){
		final int kbits=2*k;
//		verbose=true;
		if(verbose){System.err.println("Making kca from ("+fname1+", "+fname2+")\nk="+k+", gap="+gap+", matrixbits="+matrixbits+", cbits="+cbits);}
		
		boolean oldsplit=FastaReadInputStream.SPLIT_READS;
		long oldmax=maxReads;
		byte oldq=minQuality;
		maxReads=maxreads;
		minQuality=(byte)minqual;
		
		//	System.out.println("kbits="+(kbits)+" -> "+(1L<<kbits)+", matrixbits="+(matrixbits)+" -> "+(1L<<matrixbits)+", cbits="+cbits+", gap="+gap+", hashes="+hashes);
		KCountArray kca=KCountArray.makeNew(1L<<kbits, 1L<<matrixbits, cbits, gap, hashes);
		
		if(extraFiles!=null){
			for(String s : extraFiles){
				if(fileIO.FileFormat.hasFastaExtension(s)){
					assert(!FastaReadInputStream.SPLIT_READS);
				}
			}
		}
		
		if(passes==1){

			count(fname1, fname2, k, cbits, gap, rcomp, kca);
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					count(s, null, k, cbits, gap, rcomp, kca);
				}
			}
			kca.shutdown();

		}else{
			assert(passes>1);
			KCountArray trusted=null;
			for(int i=1; i<passes; i++){
				boolean conservative=i>2;// /*or, alternately, (trusted==null || trusted.capacity()>0.3)
				int step=(stepsize==1 ? 1 : stepsize+i%2);
				//			if(!conservative){step=(step+3)/4;}
				if(!conservative){step=Tools.min(3, (step+3)/4);}

				count(fname1, fname2, k, cbits, true, kca, trusted, maxreads, thresh1, step, conservative);
				if(extraFiles!=null){
					maxReads=-1;
					for(String s : extraFiles){
						count(s, null, k, cbits, true, kca, trusted, maxreads, thresh1, step, conservative);
					}
				}
				kca.shutdown();
				
				System.out.println("Trusted:   \t"+kca.toShortString());
				trusted=kca;
				kca=KCountArray.makeNew(1L<<kbits, 1L<<matrixbits, cbits, gap, hashes);

			}

			count(fname1, fname2, k, cbits, true, kca, trusted, maxreads, thresh2, stepsize, true);
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					count(s, null, k, cbits, true, kca, trusted, maxreads, thresh2, stepsize, true);
				}
			}
			kca.shutdown();
		}
		
		minQuality=oldq;
		maxReads=oldmax;
		FastaReadInputStream.SPLIT_READS=oldsplit;
		
		
		return kca;
	}
	
	public static KCountArray count(String reads1, String reads2, int k, int cbits, int gap, boolean rcomp, KCountArray count){
		assert(k<32 && k>=1 && (count!=null || k<20));
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
//		System.err.println("countFastq...  making a new cris");
		if(count==null){
			final long cells=1L<<kbits;
			if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
			count=KCountArray.makeNew(cells, cbits, gap);
		}
		assert(gap==count.gap);
		
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.start(); //4567
		}
		
		assert(cris!=null) : reads1;
		System.err.println("Started cris");
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
//		count(cris, k, rcomp, count);
//		assert(false) : THREADS;
		CountThread[] cta=new CountThread[THREADS];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(cris, k, rcomp, count);
			cta[i].start();
		}
		
		for(int i=0; i<cta.length; i++){
			CountThread ct=cta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		cris.close();
		if(verbose){System.err.println("Closed stream");}
		if(verbose){System.err.println("Processed "+readsProcessed+" reads.");}

		
		return count;
	}
	

	
	

	
	public static KCountArray count(final String reads1, final String reads2, final int k, final int cbits, final boolean rcomp,
			KCountArray count, final KCountArray trusted, final long maxReads, final int thresh, final int detectStepsize, final boolean conservative){
		
		assert(k<32 && k>=1 && (count!=null || k<20));
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		
//		System.out.println("k="+k+", kbits="+kbits+", mask="+Long.toHexString(mask)+", thresh="+thresh);
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);
		
		if(count==null){
			final long cells=1L<<kbits;
			if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
			count=KCountArray.makeNew(cells, cbits, 0);
		}
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.start(); //4567
		}
		
		assert(cris!=null) : reads1;
		System.err.println("Started cris");
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		

//		count(cris, k, rcomp, count, trusted, thresh, detectStepsize, conservative);

//		assert(false) : THREADS;
		CountThread[] cta=new CountThread[THREADS];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(cris, k, rcomp, count, trusted, thresh, detectStepsize, conservative);
			cta[i].start();
		}
		
		for(int i=0; i<cta.length; i++){
			CountThread ct=cta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		cris.close();
		if(verbose){System.err.println("Closed stream");}
		
//		System.out.println("*** after ***");
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);
		
		return count;
	}
	
	private static class CountThread extends Thread{
		
		CountThread(final ConcurrentReadInputStream cris_, final int k_, final boolean rcomp_, final KCountArray count_){
			this(cris_, k_, rcomp_, count_, null, 2, 1, true);
		}
		
		CountThread(final ConcurrentReadInputStream cris_, final int k_, final boolean rcomp_,
		final KCountArray count_, final KCountArray trusted_, final int thresh_, final int detectStepsize_, final boolean conservative_){
			cris=cris_;
			k=k_;
			rcomp=rcomp_;
			counts=count_;
			trusted=trusted_;
			thresh=thresh_;
			detectStepsize=detectStepsize_;
			conservative=conservative_;
			MAKE_NEW_ARRAY=(counts.getClass()!=KCountArray4MT.class);
		}
		
		@Override
		public void run(){
			buffer=new long[BUFFERLEN];
			
			if(trusted==null){
				count(cris, k, rcomp, counts);
			}else{
				count(cris, k, rcomp, counts, trusted, thresh, detectStepsize,  conservative);
			}
			
			synchronized(getClass()){
				keysCounted+=keysCountedLocal;
				readsProcessed+=readsProcessedLocal;
				
				if(bufflen>0){
					if(bufflen<BUFFERLEN){
						buffer=Arrays.copyOf(buffer, bufflen);
					}
					counts.increment(buffer);
				}
				buffer=null;
				bufflen=0;
			}
		}
		

		
		
		private void count(ConcurrentReadInputStream cris, int k, boolean rcomp, KCountArray count){
			assert(k<32 && k>=1 && count!=null);
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			
			if(count.gap==0){
				final int kbits=2*k;
				final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
				
				
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					//System.err.println("reads.size()="+reads.size());
					for(Read r : reads){
						readsProcessedLocal++;
						
						addRead(r, count, k, mask, rcomp);
						if(r.mate!=null){
							addRead(r.mate, count, k, mask, rcomp);
						}

					}
					//System.err.println("returning list");
					cris.returnList(ln);
					//System.err.println("fetching list");
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
			}else{
				final int k1=(k+1)/2;
				final int k2=k/2;
				final int kbits1=2*k1;
				final int kbits2=2*k2;
				final int gap=count.gap;
				final long mask1=~((-1L)<<(kbits1));
				final long mask2=~((-1L)<<(kbits2));
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					//System.err.println("reads.size()="+reads.size());
					for(Read r : reads){
						readsProcessedLocal++;

						addReadSplit(r, count, k1, k2, mask1, mask2, gap, rcomp);
						if(r.mate!=null){
							addReadSplit(r.mate, count, k1, k2, mask1, mask2, gap, rcomp);
						}

					}
					//System.err.println("returning list");
					cris.returnList(ln);
					//System.err.println("fetching list");
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
			}
			
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){System.err.println("Returned list");}
		}
		
		


		private void count(final ConcurrentReadInputStream cris, final int k, final boolean rcomp,
				final KCountArray count, final KCountArray trusted, final int thresh, final int detectStepsize, final boolean conservative){
			if(count.gap>0){countFastqSplit(cris, k, rcomp, count, trusted, thresh, detectStepsize, conservative);}
			assert(k<32 && k>=1 && (count!=null || k<20));
			final int kbits=2*k;
			final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					
					Read r2=r.mate;
					{
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r, trusted, k, thresh, detectStepsize) :
								ErrorCorrect.detectTrusted(r, trusted, k, thresh, detectStepsize));
//							System.out.println("\n"+toString(bs, r.length()));
//							System.out.println(new String(r.bases));
							for(int i=bs.nextClearBit(0); i<r.length(); i=bs.nextClearBit(i+1)){
								r.bases[i]='N';
								if(r.quality!=null){r.quality[i]=0;}
							}
//							System.out.println(new String(r.bases));
//							System.out.println("used = "+String.format(Locale.ROOT, "%.3f%%",count.usedFraction()*100));
//							System.out.println("used = "+((KCountArray4)count).cellsUsed());
//							if(bs.length()<r.length()){r=null;}
						}
//						if(r!=null){addRead(r, count, k, mask, rcomp);}
						addRead(r, count, k, mask, rcomp);
					}
					if(r2!=null){
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r2, trusted, k, thresh, detectStepsize) :
								ErrorCorrect.detectTrusted(r2, trusted, k, thresh, detectStepsize));
							for(int i=bs.nextClearBit(0); i<r2.length(); i=bs.nextClearBit(i+1)){
								r2.bases[i]='N';
								if(r2.quality!=null){r2.quality[i]=0;}
							}
						}
						addRead(r2, count, k, mask, rcomp);
					}

				}
				//System.err.println("returning list");
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){System.err.println("Returned list");}
		}
		
		
		private void countFastqSplit(final ConcurrentReadInputStream cris, final int k, final boolean rcomp,
				final KCountArray count, final KCountArray trusted, final int thresh, final int detectStepsize, final boolean conservative){
			assert(false) : cris.paired();
			assert(count.gap>0);
			assert(k<32 && k>=1 && (count!=null || k<20));
			final int kbits=2*k;
			final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
			

			final int k1=(k+1)/2;
			final int k2=k/2;
			final int kbits1=2*k1;
			final int kbits2=2*k2;
			final int gap=count.gap;
			final long mask1=~((-1L)<<(kbits1));
			final long mask2=~((-1L)<<(kbits2));
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					
					Read r2=r.mate;
					{
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r, trusted, k, thresh, detectStepsize) :
								ErrorCorrect.detectTrusted(r, trusted, k, thresh, detectStepsize));
//							System.out.println("\n"+toString(bs, r.length()));
//							System.out.println(new String(r.bases));
							for(int i=bs.nextClearBit(0); i<r.length(); i=bs.nextClearBit(i+1)){
								r.bases[i]='N';
								r.quality[i]=0;
							}
//							System.out.println(new String(r.bases));
//							System.out.println("used = "+String.format(Locale.ROOT, "%.3f%%",count.usedFraction()*100));
//							System.out.println("used = "+((KCountArray4)count).cellsUsed());
//							if(bs.length()<r.length()){r=null;}
						}
//						if(r!=null){addRead(r, count, k, mask, rcomp);}
						
						addReadSplit(r, count, k1, k2, mask1, mask2, gap, rcomp);
					}
					if(r2!=null){
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r2, trusted, k, thresh, detectStepsize) :
								ErrorCorrect.detectTrusted(r2, trusted, k, thresh, detectStepsize));
							for(int i=bs.nextClearBit(0); i<r2.length(); i=bs.nextClearBit(i+1)){
								r2.bases[i]='N';
								r2.quality[i]=0;
							}
						}
						addReadSplit(r2, count, k1, k2, mask1, mask2, gap, rcomp);
					}

				}
				//System.err.println("returning list");
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){System.err.println("Returned list");}
		}
		
		
		
		private void addRead(Read r, final KCountArray count, final int k, final long mask, boolean rcomp){
			
			if(PREJOIN && r.mate!=null && r.insert()>0){
				r.mate.reverseComplement();
				r=r.joinRead();
			}
			
			int len=0;
			long kmer=0;
			byte[] bases=r.bases;
			byte[] quals=r.quality;
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				int x=AminoAcid.baseToNumber[b];
				if(x<0 || (quals!=null && quals[i]<minQuality)){
					len=0;
					kmer=0;
				}else{
					kmer=((kmer<<2)|x)&mask;
					len++;
					if(len>=k){
						keysCountedLocal++;
//						System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));

//						System.out.println("Arrays.toString(buffer));
						buffer[bufflen]=kmer;
						bufflen++;
						if(bufflen>=buffer.length){
//							assert(false) : "Submitting "+Arrays.toString(buffer);
							count.increment(buffer);
							bufflen=0;
							if(MAKE_NEW_ARRAY){buffer=new long[BUFFERLEN];}
						}
//						count.increment(kmer);
						
//						System.out.println(" -> "+count.read(kmer));
//						System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//						array[(int)kmer]++;
//						System.out.println(" -> "+array[(int)kmer]+"\n");
//						assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
					}
				}
			}
			if(rcomp){
				r.reverseComplement();
				addRead(r, count, k, mask, false);
			}
		}
		
		private void addReadSplit(Read r, final KCountArray count, final int k1, final int k2, final long mask1, final long mask2, final int gap, boolean rcomp){
			
			if(PREJOIN && r.mate!=null && r.insert()>0){
				if(verbose){System.err.println("Prejoining "+r.numericID+" at "+r.insert());}
				r.mate.reverseComplement();
				r=r.joinRead();
			}
			
			int len=0;
			int shift=k2*2;
			long kmer1=0;
			long kmer2=0;
			byte[] bases=r.bases;
			byte[] quals=r.quality;
			
			assert(kmer1>=kmer2);
			
//			assert(false) : k1+", "+k2+", "+mask1+", "+mask2+", "+gap;
			
			if(verbose){System.err.println("Hashing read "+r.numericID+"; loop limits "+(k1+gap)+"-"+(bases.length));}
			for(int i=0, j=i+k1+gap; j<bases.length; i++, j++){
				int x1=AminoAcid.baseToNumber[bases[i]];
				int x2=AminoAcid.baseToNumber[bases[j]];
				
				if(x1<0 || x2<0 || (quals!=null && (quals[i]<minQuality || quals[j]<minQuality))){
					len=0;
					kmer1=0;
					kmer2=0;
				}else{
					kmer1=((kmer1<<2)|x1)&mask1;
					kmer2=((kmer2<<2)|x2)&mask2;
					len++;
					if(len>=k1){
						
						keysCountedLocal++;
//						System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
						
						long key=(kmer1<<shift)|kmer2;
//						System.err.println(Long.toHexString(key));
						
						if(verbose){System.err.println("Hashing key "+Long.toHexString(key)+" at length "+len);}
						
						buffer[bufflen]=key;
						bufflen++;
						if(bufflen>=buffer.length){
							count.increment(buffer);
							bufflen=0;
							if(MAKE_NEW_ARRAY){buffer=new long[BUFFERLEN];}
						}
//						count.increment(kmer);
						
						
//						System.out.println(" -> "+count.read(kmer));
//						System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//						array[(int)kmer]++;
//						System.out.println(" -> "+array[(int)kmer]+"\n");
//						assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
					}
				}
			}
			if(rcomp){
				r.reverseComplement();
				addReadSplit(r, count, k1, k2, mask1, mask2, gap, false);
			}
		}

		private final ConcurrentReadInputStream cris;
		private final int k;
		private final boolean rcomp;
		private final KCountArray counts;
		private final KCountArray trusted;
		private final int thresh;
		private final int detectStepsize;
		private final boolean conservative;
		private long keysCountedLocal=0;
		private long readsProcessedLocal=0;
		private long[] buffer;
		private int bufflen=0;
		private final boolean MAKE_NEW_ARRAY;
	}
	
}
