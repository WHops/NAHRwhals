package dna;

/**
 * @author Brian Bushnell
 * @date Jan 4, 2013
 *
 */
public class Scaffold implements Comparable<Scaffold> {
	
	public Scaffold(String name_, String assembly_, int length_){
		name=name_;
		assembly=assembly_;
		length=length_;
	}
	
	/** Assumes SAM format.
	 * e.g.<br> @SQ	SN:scaffold_0	LN:1785514	AS:build 9 */
	public Scaffold(byte[] s){
		this(new String(s).split("\t"));
	}
	
	/** Assumes SAM format.
	 * e.g.<br> @SQ	SN:scaffold_0	LN:1785514	AS:build 9 */
	public Scaffold(String s){
		this(s.split("\t"));
	}
	
	/** Assumes SAM format */
	public Scaffold(String[] split) {
		assert(split.length>2 && split[0].equals("@SQ"));
		for(String s : split){
			if(s.equals("@SQ")){
				//Do nothing
			}else if(s.startsWith("SN:")){
				assert(name==null);
				name=new String(s.substring(3)); //Data.forceIntern(s.substring(3));
			}else if(s.startsWith("LN:")){
				length=Integer.parseInt(s.substring(3));
			}else if(s.startsWith("AS:")){
				assembly=Data.forceIntern(s.substring(3));
			}
		}
		assert(length>-1);
		assert(name!=null);
	}
	
	public Scaffold(String name_, int length_) {
		name=name_;
		length=length_;
	}
	
	@Override
	public int hashCode(){
		return name.hashCode();
	}
	
	@Override
	public int compareTo(Scaffold other){
		return name.compareTo(other.name);
	}
	
	@Override
	public String toString(){
		return "@SQ\tSN:"+name+"\tLN:"+length+(assembly==null ? "" : "\tAS:"+assembly);
	}
	
	public String name;
	public String assembly;
	public int length=-1;
	public long basehits=0;
	public long readhits=0;
	/** For calculating FPKM */
	public long fraghits=0;
	public long readhitsMinus=0;
	
	/** {A,C,G,T,N} */
	public long[] basecount;
	public float gc;
	
	/** For attaching things */
	public Object obj1;
	
	/** For attaching more things */
	public Object obj2;
	
}
