 package var2;

import shared.Parse;
import shared.Shared;
import shared.Tools;
import stream.SamLine;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.CoverageArray3;
import structures.CoverageArray3A;

public class Scaffold {
	
	/** Assumes SAM format.
	 * e.g.<br> @SQ	SN:scaffold_0	LN:1785514	AS:build 9 */
	public Scaffold(byte[] line, int scafnum){
		assert(Tools.startsWith(line, "@SQ\t")) : new String(line);
		number=scafnum;
		int a=0, b=0;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		assert(Tools.startsWith(line, "SN:", a));
		name=new String(line, a+3, b-a-3);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		assert(Tools.startsWith(line, "LN:", a));
		length=Parse.parseInt(line, a+3, b);
		b++;
		a=b;
	}
	
	public Scaffold(String name_, int scafnum_, int len_){
		name=name_;
		number=scafnum_;
		length=len_;
	}
	
	public void add(SamLine sl){
		int start=sl.pos-1;
		int stop=sl.stop(start, false, false);
		increment(start, stop, sl.strand());
	}
	
	public void increment(int from, int to, int strand){
//		assert(trackStrand);
		if(!initialized()){
			synchronized(this){
				if(!initialized()){
//					assert(ca==null);
//					assert(caMinus==null);
//					assert(trackStrand);
					ca=useCA3A ? new CoverageArray3A(number, length) : useCA3 ? new CoverageArray3(number, length) : new CoverageArray2(number, length);
					if(trackStrand){
						caMinus=useCA3A ? new CoverageArray3A(number, length) : useCA3 ? new CoverageArray3(number, length) : new CoverageArray2(number, length);
					}
				}
				initialized=true;
			}
//			assert(ca!=null);
//			assert(!trackStrand || caMinus!=null) : trackStrand;
		}
		assert(initialized());
		assert(ca!=null);
//		assert(!trackStrand || caMinus!=null) : trackStrand;
//		assert(trackStrand);
		ca.incrementRangeSynchronized(from, to, 1);
		if(trackStrand && strand==Shared.MINUS){caMinus.incrementRangeSynchronized(from, to, 1);}
	}
	
	public synchronized void incrementOld(int from, int to, int strand){
		if(ca==null){
			ca=useCA3 ? new CoverageArray3(number, length) : new CoverageArray2(number, length);
		}
		ca.incrementRange(from, to);
		if(trackStrand && strand==Shared.MINUS){
			if(caMinus==null){
				caMinus=useCA3 ? new CoverageArray3(number, length) : new CoverageArray2(number, length);
			}
			caMinus.incrementRange(from, to);
		}
	}
	
	public String getSequence(SamLine sl) {
		int start=sl.start(false, false);
		int stop=sl.stop(start, false, false);
		return getSequence(start, stop);
	}
	
	public String getSequence(int start, int stop) {
		assert(bases!=null) : this;
		start=Tools.max(0, start);
		stop=Tools.min(bases.length-1, stop);
		String s=new String(bases, start, stop-start+1);
		return s;
	}
	
	public int calcCoverage(Var v){
		return calcCoverage(v, ca);
	}
	
	public int minusCoverage(Var v){
		assert(trackStrand);
		return calcCoverage(v, caMinus);
	}
	
	public int calcCoverage(Var v, CoverageArray ca){
		final int a=v.start;
		final int b=v.stop;
//		assert(false) : ca.maxIndex+", "+a;
		if(ca==null || ca.maxIndex<a){return 0;}
		final int type=v.type();
		final int avg;
		final int rlen=v.reflen();
		long sum=0;
		if(type==Var.SUB || type==Var.NOCALL || type==Var.DEL){
			for(int i=a; i<b; i++){
				sum+=ca.get(i);
			}
			avg=(int)Math.round(sum/(double)rlen);
		}else if(type==Var.INS){
			assert(rlen==0 && a==b);
//			if(a<=0){sum=2*ca.get(0);}
//			else if(b>ca.maxIndex)
			if(b>=ca.maxIndex){
				sum=2*ca.get(ca.maxIndex);
				avg=(int)(sum/2);
			}else{
				sum=ca.get(a)+ca.get(b);
				avg=(int)Math.ceil(sum/2);
			}
		}else if(type==Var.LJUNCT){
			//Take coverage from right of junction, unless that's off the end.
			//but it should be impossible for that to be off the end.
			avg=ca.get(Tools.min(ca.maxIndex, a+1));
		}else if(type==Var.RJUNCT){
			//Take coverage from left of junction, unless that's off the end.
			//but it should be impossible for that to be off the end.
			avg=ca.get(Tools.max(0, a-1));
		}else{
			throw new RuntimeException("Unknown type "+type+"\n"+v);
		}
		return avg;
	}
	
	@Override
	public String toString(){
		return "@SQ\tSN:"+name+"\tLN:"+length+"\tID:"+number;
	}
	
	public synchronized void clearCoverage(){
		ca=null;
		caMinus=null;
		initialized=false;
	}
	
	public final String name;
	public final int number;
	public final int length;
	private CoverageArray ca;
	private CoverageArray caMinus;
	public byte[] bases;
	private boolean initialized(){return initialized;};
	private boolean initialized;

	public static void setCA3(boolean b){useCA3=b;}
	public static void setCA3A(boolean b){useCA3A=b;}
	public static void setTrackStrand(boolean b){trackStrand=b;}
	public static boolean trackStrand(){return trackStrand;}

	private static boolean useCA3=false;
	private static boolean useCA3A=true;
	private static boolean trackStrand=false;
	
}
