package bloom;

import java.util.ArrayList;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.FileFormat;
import shared.Timer;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Jul 5, 2012
 *
 */
public class KmerCount4 extends KmerCountAbstract {
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		
		String fname1=args[0];
		String fname2=(args.length>3 || args[1].contains(".") ? args[1] : null);
		int k=14;
		int cbits=16;
		int gap=0;
		
		for(int i=(fname2==null ? 1 : 2); i<args.length; i++){
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
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		KCountArray2 count=null;
		
		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			assert(!FastaReadInputStream.SPLIT_READS);
			FastaReadInputStream.MIN_READ_LEN=k;
		}
		
		if(gap==0){
			count=count(fname1, fname2, k, cbits, true);
		}else{
			count=countFastqSplit(fname1, fname2, (k+1)/2, k/2, gap, cbits, true, null);
		}
		
		
		t.stop();
		System.out.println("Finished counting; time = "+t);
		
		printStatistics(count);
		
	}

	public static void printStatistics(KCountArray2 count){
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
	
	public static KCountArray2 count(String reads1, String reads2, int k, int cbits, boolean rcomp){
		return count(reads1, reads2, k, cbits, rcomp, null);
	}
	
	public static KCountArray2 count(String reads1, String reads2, int k, int cbits, boolean rcomp, KCountArray2 count){
		assert(k>=1 && k<20);
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		
		if(count==null){
			final long cells=1L<<kbits;
			if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
			count=new KCountArray2(cells, cbits);
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
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert(paired==(r.mate!=null));
		}

		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			//System.err.println("reads.size()="+reads.size());
			for(Read r : reads){
				readsProcessed++;

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
		if(verbose){System.err.println("Finished reading");}
		cris.returnList(ln);
		if(verbose){System.err.println("Returned list");}
		cris.close();
		if(verbose){System.err.println("Closed stream");}
		if(verbose){System.err.println("Processed "+readsProcessed+" reads.");}

		
		return count;
	}
	
	public static KCountArray2 countFastqSplit(String reads1, String reads2, int k1, int k2, int gap, int cbits, boolean rcomp, KCountArray2 count){
		assert(k1+k2>=1 && k1+k2<20);
		assert(gap>=0);
		final int kbits1=2*k1;
		final int kbits2=2*k2;
		final long mask1=~((-1L)<<(kbits1));
		final long mask2=~((-1L)<<(kbits2));
		
		if(count==null){
			final long cells=1L<<(kbits1+kbits2);
			if(verbose){System.err.println("k1="+k1+", k2="+k2+", kbits1="+kbits1+", kbits2="+kbits2+", cells="+cells+
					", mask1="+Long.toHexString(mask1)+", mask2="+Long.toHexString(mask2));}
			count=new KCountArray2(cells, cbits, gap);
		}
		assert(count.gap==gap);
		
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
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert(paired==(r.mate!=null));
		}

		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			//System.err.println("reads.size()="+reads.size());
			for(Read r : reads){
				readsProcessed++;

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
		if(verbose){System.err.println("Finished reading");}
		cris.returnList(ln);
		if(verbose){System.err.println("Returned list");}
		cris.close();
		if(verbose){System.err.println("Closed stream");}
		if(verbose){System.err.println("Processed "+readsProcessed+" reads.");}

		
		return count;
	}
	
	public static void addRead(final Read r, final KCountArray2 count, final int k, final long mask, boolean rcomp){
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
					keysCounted++;
//					System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
					count.increment(kmer, 1);
//					System.out.println(" -> "+count.read(kmer));
//					System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//					array[(int)kmer]++;
//					System.out.println(" -> "+array[(int)kmer]+"\n");
//					assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
				}
			}
		}
		if(rcomp){
			r.reverseComplement();
			addRead(r, count, k, mask, false);
		}
	}
	
	public static void addReadSplit(final Read r, final KCountArray2 count, final int k1, final int k2, final long mask1, final long mask2, final int gap, boolean rcomp){
		int len=0;
		int shift=k2*2;
		long kmer1=0;
		long kmer2=0;
		byte[] bases=r.bases;
		byte[] quals=r.quality;
		
		assert(kmer1>=kmer2);
		
//		assert(false) : k1+", "+k2+", "+mask1+", "+mask2+", "+gap;
		
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
					keysCounted++;
//					System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
					
					long key=(kmer1<<shift)|kmer2;
//					System.err.println(Long.toHexString(key));
					count.increment(key, 1);
//					System.out.println(" -> "+count.read(kmer));
//					System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//					array[(int)kmer]++;
//					System.out.println(" -> "+array[(int)kmer]+"\n");
//					assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
				}
			}
		}
		if(rcomp){
			r.reverseComplement();
			addReadSplit(r, count, k1, k2, mask1, mask2, gap, false);
		}
	}
	
	public static void addReadSplit(final byte[] bases, final KCountArray2 count, final int k1, final int k2, final long mask1, final long mask2, final int gap, boolean rcomp){
		int len=0;
		int shift=k2*2;
		long kmer1=0;
		long kmer2=0;
		byte[] quals=null;
		
		assert(kmer1>=kmer2);
		
//		assert(false) : k1+", "+k2+", "+mask1+", "+mask2+", "+gap;
		
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
					keysCounted++;
//					System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
					
					long key=(kmer1<<shift)|kmer2;
					System.out.println(Long.toHexString(kmer1));
					System.out.println(Long.toHexString(kmer2));
					System.out.println(Long.toHexString(key));
					count.increment(key, 1);
//					System.out.println(" -> "+count.read(kmer));
//					System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//					array[(int)kmer]++;
//					System.out.println(" -> "+array[(int)kmer]+"\n");
//					assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
				}
			}
		}
		if(rcomp){
			AminoAcid.reverseComplementBasesInPlace(bases);
			addReadSplit(bases, count, k1, k2, mask1, mask2, gap, false);
		}
	}
	
}
