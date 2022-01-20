package var2;

//This class is not currently thought to be used
/** Allows simpler hashing, without complete alleles. */
public class VarKey implements Comparable<VarKey> {
	
	public static VarKey toVarKey(Var v){
		if(v.type==Var.INS){
			return new VarKey(v.scafnum, v.start, v.allele.length, v.type, v.allele[0]);
		}else if(v.type==Var.DEL){
			return new VarKey(v.scafnum, v.start, v.reflen(), v.type, 0);
		}
		return new VarKey(v.scafnum, v.start, v.reflen(), v.type, v.allele[0]);
	}
	
	public VarKey(int scafNum_, int start_, int length_, int type_, int allele_){
		scafNum=scafNum_;
		start=start_;
		length=length_;
		type=type_;
		allele=allele_;
	}
	
	@Override
	public int hashCode(){
		return scafNum^Integer.rotateLeft(start, 4)^Integer.rotateRight(start, 18)^Integer.rotateLeft(type, 8)^Integer.rotateLeft(allele, 12);
	}
	
	@Override
	public boolean equals(Object b){
		return equals((VarKey)b);
	}
	
	public boolean equals(VarKey b){
		if(b==null){return false;}
		return scafNum==b.scafNum && start==b.start && length==b.length && type==b.type && allele==b.allele;
	}
	
	@Override
	public int compareTo(VarKey b){
		if(b==null){return -1;}
		if(scafNum!=b.scafNum){return scafNum-b.scafNum;}
		if(start!=b.start){return start-b.start;}
		if(length!=b.length){return length-b.length;}
		if(type!=b.type){return type-b.type;}
		if(allele!=b.allele){return allele-b.allele;}
		return 0;
	}
	
	int scafNum;
	int start;
	int length;
	int type;
	//0 for DEL, otherwise first letter of allele 
	int allele;
}
