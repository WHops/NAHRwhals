package structures;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Tools;


public abstract class CoverageArray implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -7175422489330746676L;
	
	
	public static final CoverageArray read(String fname){
		
		if(!fname.contains(".ca")){
			throw new RuntimeException();
//			ca=new CoverageArray2();
//			ca.load(new TsvCoverageFile(fname));
//			return ca;
		}
		
		fname=ReadWrite.findFileExtension(fname);
//		System.err.println("Found "+fname);
		
		return ReadWrite.read(CoverageArray.class, fname, true);

//		if(fname.endsWith(".ca2") || fname.contains(".ca2.")){return ReadWrite.read(CoverageArray2.class, fname);}
//		else if(fname.endsWith(".ca") || fname.contains(".ca.")){return ReadWrite.read(CoverageArray1.class, fname);}
//		else{return ReadWrite.read(CoverageArray.class, fname);}
	}
	
	public CoverageArray(int chrom){chromosome=chrom;}
	
	/**
	 * @param loc
	 * @param amt
	 */
	public abstract void increment(int loc, int amt);
	
	/**
	 * @param loc
	 */
	public abstract void increment(int loc);

	public final void incrementRange(int min, int max){incrementRange(min, max, 1);}
	public abstract void incrementRange(int min, int max, int amt);
	public abstract void incrementRangeSynchronized(int min, int max, int amt);
	
	public void incrementRanges(IntList ranges, int amt){
		for(int i=0; i<ranges.size; i+=2){
			int a=ranges.get(i), b=ranges.get(i+1);
			incrementRange(a, b-1, 1);
		}
	}
	
	public abstract void set(int loc, int val);
	
	public abstract int get(int loc);
	
	public abstract void resize(int newlen);
	
	
	public final double[][] toGraph(int blocksize, int min, int max){
		
		min=max(min, minIndex);
		max=min(max, maxIndex);
		int length=max-min;
		
		ArrayList<double[]> list=new ArrayList<double[]>();
		
		int block;
		
		if(blocksize<=0){
//			block=((array.length+62999)/63000);//For Excel
//			block=((length+62999)/63000);//For Excel
			block=((length+31499)/31500);//For Excel
		}else{
			block=blocksize;
		}
		block=max(block, 1);
		
		int current=0;
		double[] sum=new double[2];
		for(int loc=min; loc<=max; loc++){
			if(current==block){
				for(int i=0; i<sum.length; i++){
					sum[i]=sum[i]/current;
				}
				sum[0]=Math.round(sum[0]);
				list.add(sum);
				sum=new double[2];
				current=0;
			}
			
			sum[0]+=loc;
			sum[1]+=get(loc);
			
			current++;
		}
		
		return list.toArray(new double[0][]);
		
	}
	
	
	public static final void print(double[][] data){
		
//		data=stats.Smoother.weightedAveragePlank(data, 24);
		assert(false) : "Smoother disabled in this code purely to reduce dependancies.";
		StringBuilder sb=new StringBuilder(data.length*20);
		for(double[] d : data){
			sb.append(String.format(Locale.ROOT, "%d\t%.2f\n",(int)d[0],d[1]));
		}
		System.out.print(sb);
	}
	
	public static CoverageArray makeArray(int num, int size, Class<? extends CoverageArray> c){
		if(c==CoverageArray2.class){
			return new CoverageArray2(num, size);
		}else if(c==CoverageArray3.class){
			return new CoverageArray3(num, size);
		}else if(c==CoverageArray3A.class){
			return new CoverageArray3(num, size);
		}
		throw new RuntimeException("Unhandled class: "+c);
	}
	
	//TODO: Extremely slow due to string processing
	public static HashMap<String, CoverageArray> loadDepth(FileFormat ffdepth, Class<? extends CoverageArray> c) {
		ByteFile bf=ByteFile.makeByteFile(ffdepth);
		HashMap<String, CoverageArray> map=new HashMap<String, CoverageArray>();
		
		String prevName=null;
		CoverageArray prevArray=null;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line[0]!='#'){
//				ArrayList<byte[]> split=Tools.split(line, 0, (byte)'\t');
				String[] split=Tools.tabPattern.split(new String(line));
				String name=split[0];
				int pos=Integer.parseInt(split[1]);
				int depth=Integer.parseInt(split[2]);
				CoverageArray current;
				if(name.equals(prevName)){
					current=prevArray;
				}else{
					assert(!map.containsKey(name)) : name; //Could do a lookup but should not be needed
					current=makeArray(map.size()+1, 64, c);
					map.put(name, current);
					prevName=name;
					prevArray=current;
				}
				if(depth>0){current.set(pos, depth);}
			}
		}
		return map;
	}
	
	@Override
	public abstract String toString();
	
	static final long min(long x, long y){return x<y ? x : y;}
	static final long max(long x, long y){return x>y ? x : y;}
	static final int min(int x, int y){return x<y ? x : y;}
	static final int max(int x, int y){return x>y ? x : y;}
	
	public int chromosome;
	
	public int maxIndex=-1;
	public int minIndex=Integer.MAX_VALUE;
	public int length(){return maxIndex-minIndex+1;}
	public abstract int arrayLength();
	
	private static boolean OVERFLOWED=false;
	
}
