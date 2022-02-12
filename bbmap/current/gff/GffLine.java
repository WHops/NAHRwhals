package gff;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Locale;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import prok.ProkObject;
import shared.Parse;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import var2.ScafMap;
import var2.VCFLine;
import var2.Var;

/**
 * Used by both the var2 and prok packages for processing gff files.
 * @author Brian Bushnell
 * @date Sep 12, 2018
 *
 */
public class GffLine {
	
	//#seqid	source	type	start	end	score	strand	phase	attributes
	public GffLine(byte[] line){
		int a=0, b=0;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		seqid=parseSeqid ? intern(new String(line, a, b-a)) : null;
//		assert(seqid==null || seqid.equals(new String(line, a, b-a)));
//		assert(seqid!=null) : new String(line, a, b-a)+", "+a+", "+b+"\n"+line;
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		if(b==a+1 && line[a]=='.'){source=DOTS;}
		else{source=paseSource ? intern(new String(line, a, b-a)) : null;}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		if(b==a+1 && line[a]=='.'){type=DOTS;}
		else{
			try {//This was to catch a probably intermittent hardware error; can't replicate.
				type=(parseType ? intern(new String(line, a, b-a)) : null);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.err.println("\n"+new String(line)+"\n"+a+", "+b+", "+(b-a));
				assert(false);
			}
		}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		start=Parse.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		stop=Parse.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		if(b<=a){
			//Badly formatted line; common in IMG
			return;
		}
		assert(b>a) : "Missing field 5: "+new String(line);
		if(b==a+1 && line[a]=='.'){score=-1;}
		else{score=Parse.parseFloat(line, a, b);}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 6: "+new String(line);
		assert(b==a+1);
		strand=find(line[a], STRANDS);
//		assert(strand>0) : line[a]+", "+Arrays.toString(STRANDS)+", "+(char)line[b];
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 7: "+new String(line);
		assert(b==a+1);
		if(line[a]=='.'){phase=-1;}
		else{phase=Parse.parseInt(line, a, b);}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 8: "+new String(line);
		if(b==a+1 && line[a]=='.'){attributes=DOTS;}
		else{attributes=parseAttributes ? new String(line, a, b-a) : null;}
		b++;
		a=b;
		
//		assert(strand>=0) : "\n"+this.toString()+"\n"+new String(line);
	}
	
	public GffLine(VCFLine vcf){
		seqid=vcf.scaf;
		source=DOTS;
		type="sequence_variant_obs";
		start=vcf.start()+1;
		stop=vcf.stop()+1;
		score=(float)vcf.qual;
		strand=PLUS;
		phase=-1;
		final int vtype=vcf.type();
		ByteBuilder bb=new ByteBuilder(16);
		bb.append("ID=").append(Var.typeArray[vtype]).append(' ');
		if(vtype==Var.SUB){
			bb.append(vcf.ref).append('>').append(vcf.alt);
		}else if(vtype==Var.DEL){
			bb.append("length ").append(vcf.reflen()-vcf.readlen());
		}else if(vtype==Var.INS){
			int offset=vcf.reflen();
			int length=vcf.readlen()-offset;
			bb.append(vcf.alt, offset, length);
		}else if(vtype==Var.NOCALL){
			bb.append("length ").append(vcf.reflen());
		}
		attributes=bb.toString();
		bb.clear();
	}
	
	public GffLine(Var v, double properPairRate, double totalQualityAvg, double totalMapqAvg, double readLengthAvg, double rarity, int ploidy, ScafMap map){
		seqid=v.scafName();
		source=DOTS;
		type="sequence_variant_obs";
		start=v.start+1;
		stop=Tools.max(v.start+1, v.stop);
		score=(float)v.score(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map);
		strand=PLUS;
		phase=-1;
		final int vtype=v.type();
		ByteBuilder bb=new ByteBuilder(16);
		bb.append("ID=").append(Var.typeArray[vtype]);
		if(vtype==Var.SUB || vtype==Var.INS){
			bb.append(' ').append(v.allele);
		}else if(vtype==Var.DEL || vtype==Var.NOCALL){
			bb.append(" length ").append(v.reflen());
		}else{assert(false) : vtype+"\n"+v;}
		attributes=bb.toString();
		bb.clear();
	}
	
	public GffLine(Var v){
		seqid=v.scafName();
		source="BBTools";
		type="sequence_variant_obs";
		start=v.start+1;
		stop=Tools.max(v.start+1, v.stop);
		score=-1;
		strand=PLUS;
		phase=-1;
		final int vtype=v.type();
		ByteBuilder bb=new ByteBuilder(16);
		bb.append("ID=").append(Var.typeArray[vtype]);
		if(vtype==Var.SUB || vtype==Var.INS){
			bb.append(' ').append(v.allele);
		}else if(vtype==Var.DEL || vtype==Var.NOCALL){
			bb.append(" length ").append(v.reflen());
		}else{assert(false) : vtype+"\n"+v;}
		attributes=bb.toString();
		bb.clear();
	}
	
	public static ArrayList<GffLine> loadGffFile(String fname, String types, boolean banUnprocessed){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.GFF, null, false, false);
		return loadGffFile(ff, types, banUnprocessed);
	}
	
	public static ArrayList<GffLine>[] loadGffFileByType(FileFormat ff, String types, boolean banUnprocessed){
		ArrayList<GffLine> list=loadGffFile(ff, types, banUnprocessed);
		String[] typeArray=types.split(",");
		ArrayList<GffLine>[] lists=new ArrayList[typeArray.length];
		for(int i=0; i<typeArray.length; i++){
			String type=typeArray[i];
			lists[i]=new ArrayList<GffLine>();
			for(GffLine gline : list){
				if(gline.type.equals(type)){
					lists[i].add(gline);
				}
			}
		}
		return lists;
	}
	
	public static ArrayList<GffLine> loadGffFile(FileFormat ff, String types, boolean banUnprocessed){
		HashSet<String> set=null;
		if(types!=null){
			String[] split=types.split(",");
			set=new HashSet<String>(split.length*2);
			for(String s : split){
				set.add(s);
			}
		}
		
		ArrayList<GffLine> list=new ArrayList<GffLine>();
		ByteFile bf=ByteFile.makeByteFile(ff);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line[0]=='#'){
				//skip
			}else{
				GffLine gline=new GffLine(line);
				assert(gline.strand>=0) : "\n"+gline.toString()+"\n"+new String(line)+"\n";
				if(set==null || (gline.type!=null && set.contains(gline.type))){
					if(!banUnprocessed || ProkObject.processType(gline.prokType())){
						list.add(gline);
					}
				}
			}
		}
		
		boolean error=bf.close();
		assert(!error) : "Problem with file "+ff.name();
		return list;
	}
	
	public static void toText(ByteBuilder bb, Var v, double properPairRate, double totalQualityAvg, 
			double totalMapqAvg, double readLengthAvg, double rarity, int ploidy, ScafMap map){
//		assert(false);
		bb.append(v.scafName(map)).append('\t');
		bb.append('.').append('\t');
		bb.append("sequence_variant_obs").append('\t');
		bb.append(v.start+1).append('\t');
		bb.append(Tools.max(v.start+1, v.stop)).append('\t');
		bb.append(v.score(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map), 2).append('\t');
		bb.append('+').append('\t');
		bb.append('.').append('\t');
//		System.err.println(v.typeString()+", "+v.start+", "+v.stop);
		final int vtype=v.type();
		bb.append("ID=").append(Var.typeArray[vtype]);
		if(vtype==Var.SUB || vtype==Var.INS){
			bb.append(' ').append(v.allele);
		}else if(vtype==Var.DEL || vtype==Var.NOCALL){
			bb.append(" length ").append(v.reflen());
		}else{assert(false) : vtype+"\n"+v;}
	}
	
	public static String toHeader(double properPairRate, double totalQualityAvg, double mapqAvg, double rarity, double minAlleleFraction, int ploidy, 
			long reads, long pairs, long properPairs, long bases, String ref){
		StringBuilder sb=new StringBuilder();
		
		final double readLengthAvg=bases/Tools.max(1.0, reads);
		sb.append("##gff-version 3\n");
		sb.append("#BBMapVersion\t"+Shared.BBMAP_VERSION_STRING+"\n");
		sb.append("#ploidy\t"+ploidy+"\n");
		sb.append(String.format(Locale.ROOT, "#rarity\t%.5f\n", rarity));
		sb.append(String.format(Locale.ROOT, "#minAlleleFraction\t%.4f\n", minAlleleFraction));
		sb.append("#reads\t"+reads+"\n");
		sb.append("#pairedReads\t"+pairs+"\n");
		sb.append("#properlyPairedReads\t"+properPairs+"\n");
		sb.append(String.format(Locale.ROOT, "#readLengthAvg\t%.2f\n", readLengthAvg));
		sb.append(String.format(Locale.ROOT, "#properPairRate\t%.4f\n", properPairRate));
		sb.append(String.format(Locale.ROOT, "#totalQualityAvg\t%.4f\n", totalQualityAvg));
		sb.append(String.format(Locale.ROOT, "#mapqAvg\t%.2f\n", mapqAvg));
		if(ref!=null){sb.append("#reference\t"+ref+"\n");}
		
		sb.append("#seqid	source	type	start	end	score	strand	phase	attributes");
		return sb.toString();
	}
	
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		appendTo(bb);
		return bb.toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append(seqid==null ? "." : seqid).append('\t');
		bb.append(source==null ? "." : source).append('\t');
		bb.append(type==null ? "." : type).append('\t');
		bb.append(start).append('\t');
		bb.append(stop).append('\t');
		if(score<0){bb.append('.').append('\t');}
		else{bb.append(score, 2).append('\t');}
		
		bb.append((strand>=0 ? STRANDS[strand] : (byte)'.')).append('\t');
		
		if(phase<0){bb.append('.').append('\t');}
		else{bb.append(phase).append('\t');}
		
		bb.append(attributes==null ? "." : attributes);
		return bb;
	}
	
	public int length() {
		return stop-start+1;
	}
	
	private static int find(byte a, byte[] array){
		for(int i=0; i<array.length; i++){
			if(array[i]==a){return i;}
		}
		return -1;
	}
	
	private static String intern(String s){
		return Data.forceIntern(s);
	}
	
	@Override
	public int hashCode(){
		return trueStop()^seqid.hashCode();
	}
	
	@Override
	public boolean equals(Object o){
		GffLine b=(GffLine)o;
		if(start!=b.start){return false;}
		if(stop!=b.stop){return false;}
		if(strand!=b.strand){return false;}
		if(!seqid.equals(b.seqid)){return false;}
		if(!type.equals(b.type)){return false;}
		return true;
	}
	
	public int trueStart(){
		return strand==0 ? start : stop;
	}
	
	public int trueStop(){
		return strand==0 ? stop : start;
	}
	
	public final int prokType(){
		if(type.equals("CDS")){
			return ProkObject.CDS;
		}else if(type.equals("tRNA")){
			return ProkObject.tRNA;
		}else if(type.equals("rRNA")){
			if(attributes.contains("16S")){
				return ProkObject.r16S;
			}else if(attributes.contains("23S")){
				return ProkObject.r23S;
			}else if(attributes.contains("18S")){
				return ProkObject.r18S;
			}else if(attributes.contains("5S") && length()<300){
				return ProkObject.r5S;
			}
		}
		return -1;
	}
	
	public final boolean partial(){return attributes!=null && attributes.contains("partial=true");}
	
	public final boolean inbounds(int scaflen){return start>=0 && stop<scaflen;}
	
	public String seqid;
	public String source;
	public String type;
	public int start;
	public int stop;
	public float score;
	public int strand;
	public int phase;
	public String attributes;
	
	private static final byte[] STRANDS=new byte[] {'+', '-', '?', '.'};
	public static final int PLUS=0, MINUS=1, QMARK=2, DOT=3;
	public static final String DOTS=".";
	
	public static boolean parseSeqid=true;
	public static boolean paseSource=false;
	public static boolean parseType=true;
	public static boolean parseScore=false;
	public static boolean parseAttributes=true;
	
//	public static boolean parseSeqid=true;
//	public static boolean paseSource=true;
//	public static boolean parseType=true;
//	public static boolean parseScore=true;
//	public static boolean parseAttributes=true;
	
}
