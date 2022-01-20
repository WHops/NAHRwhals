package tax;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Shared;
import shared.Tools;
import structures.IntList;

/**
 * @author Brian Bushnell
 * @date Mar 10, 2015
 *
 */
public class GiToTaxid {
	
	public static void main(String[] args){
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		ReadWrite.ZIPLEVEL=9;
		ReadWrite.PIGZ_BLOCKSIZE=256;
//		ReadWrite.PIGZ_ITERATIONS=30;
		
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			shared.Parser.parseZip(arg, a, b);
		}
//		if(args.length>2 && false){//Run a test
//			test(args);
//		}else 
		if(args.length>=2){//Write array
			initialize(args[0]);
			ReadWrite.write(array, args[1], true);
		}
	}
	
	public static void test(String[] args){
		System.err.println(getID(1000));
		System.err.println(getID(10000));
		System.err.println(getID(10001));
		System.err.println(getID(10002));
		System.err.println(getID(10003));
		System.err.println(getID(10004));
		System.err.println(getID(10005));
		System.err.println(getID(100000));
		System.err.println(getID(1000000));
		System.err.println(getID(10000000));
		
		TaxTree tree=null;
		if(args.length>1){
			tree=TaxTree.loadTaxTree(args[0], System.err, true, true);
		}
		
		System.err.println("Strings:");
		int x;
		x=getID("gi|18104025|emb|AJ427095.1| Ceratitis capitata centromeric or pericentromeric satellite DNA, clone 44");
		System.err.println(x);
		if(tree!=null){
			System.err.println(tree.getNode(x));
			tree.incrementRaw(x, 30);
		}
		x=getID("gi|15982920|gb|AY057568.1| Arabidopsis thaliana AT5g43500/MWF20_22 mRNA, complete cds");
		System.err.println(x);
		if(tree!=null){
			System.err.println(tree.getNode(x));
			tree.incrementRaw(x, 40);
		}
		x=getID("gi|481043749|gb|KC494054.1| Plesiochorus cymbiformis isolate ST05-58 internal transcribed spacer 2, partial sequence");
		System.err.println(x);
		if(tree!=null){
			System.err.println(tree.getNode(x));
			tree.incrementRaw(x, 20);
		}
		
		if(tree!=null){
			tree.percolateUp();
			ArrayList<TaxNode> nodes=tree.gatherNodesAtLeastLimit(35);
			for(TaxNode n : nodes){
				System.err.println(n);
			}
		}
	}
	
	public static int parseGiToTaxid(String s){return parseGiToTaxid(s, '|');}
	public static int parseGiToTaxid(String s, char delimiter){
		int x=parseGiNumber(s, delimiter);
		assert(x>=0) : s;
		assert(array!=null) : "To use gi numbers, you must load a gi table.";
//		if(x>=array.length || array[x]<0){x=(int)(Math.random()*array.length);} //Test to make sure array is nonempty.
		if(x>=0 && x<array.length){return array[x];}
		assert(x<array.length) : "The GI number "+x+" is too big.\n"
				+ "Please update the gi table with the latest version from NCBI as per the instructions in gitable.sh.\n"
				+ "To ignore this problem, please run with the -da flag.\n";
		return -1;
	}
	

	public static int parseGiToTaxid(byte[] s){return parseGiToTaxid(s, '|');}
	public static int parseGiToTaxid(byte[] s, char delimiter){
		int x=parseGiNumber(s, delimiter);
		if(x>=0){return array[x];}
		return -1;
	}
	
	/** Parse a gi number, or return -1 if formatted incorrectly. */
	static int parseGiNumber(String s, char delimiter){
		if(s==null || s.length()<4){return -1;}
//		System.err.println("a");
		if(s.charAt(0)=='>'){return getID(s.substring(1), delimiter);}
//		System.err.println("b");
		if(!s.startsWith("gi")){return -1;}
//		System.err.println("c");
//		System.err.println("d");
		int initial=s.indexOf(delimiter);
//		System.err.println("e");
		if(initial<0){
			if(delimiter!='~'){
				delimiter='~';
				initial=s.indexOf(delimiter);
			}
			if(initial<0){
				delimiter='_';
				initial=s.indexOf(delimiter);
			}
			if(initial<0){return -1;}
//			System.err.println("f");
//			System.err.println("g");
		}
//		System.err.println("h");
		if(!Tools.isDigit(s.charAt(initial+1))){return -1;}
//		System.err.println("i");
		
		int number=0;
		for(int i=initial+1; i<s.length(); i++){
			char c=s.charAt(i);
			if(c==delimiter){break;}
			assert(Tools.isDigit(c));
			number=(number*10)+(c-'0');
		}
//		System.err.println("j: "+number);
		return number;
	}
	
	/** Parse a ncbi number, or return -1 if formatted incorrectly. */
	public static int parseTaxidNumber(String s, char delimiter){
		if(s==null || s.length()<5){return -1;}
		if(s.charAt(0)=='>'){return parseTaxidNumber(s.substring(1), delimiter);}
		if(!s.startsWith("ncbi") && !s.startsWith("tid")){return -1;}
		int initial=s.indexOf(delimiter);
		if(initial<0){
			delimiter='_';
			initial=s.indexOf(delimiter);
			if(initial<0){return -1;}
		}
		if(!Tools.isDigit(s.charAt(initial+1))){return -1;}
		
		int number=0;
		for(int i=initial+1; i<s.length(); i++){
			char c=s.charAt(i);
			if(c==delimiter || c==' '){break;}
			assert(Tools.isDigit(c)) : c+"\n"+s;
			number=(number*10)+(c-'0');
		}
		return number;
	}
	

	public static int getID(String s){return getID(s, '|');}
	/** Get the taxID from a header starting with a taxID or gi number */
	public static int getID(String s, char delimiter){
		int x=parseTaxidNumber(s, delimiter);
		if(x>=0){return x;}
		x=parseGiNumber(s, delimiter);
		if(x>=0){return array[x];}
		return -1;
	}
	
	/** Parse a gi number, or return -1 if formatted incorrectly. */
	static int parseGiNumber(byte[] s, char delimiter){
		if(s==null || s.length<4){return -1;}
		if(!Tools.startsWith(s, "gi") && !Tools.startsWith(s, ">gi")){return -1;}
		int initial=Tools.indexOf(s, (byte)delimiter);
		if(initial<0){
			delimiter='_';
			initial=Tools.indexOf(s, (byte)delimiter);
			if(initial<0){return -1;}
		}
		if(!Tools.isDigit(s[initial+1])){return -1;}
		
		int number=0;
		for(int i=initial+1; i<s.length; i++){
			byte c=s[i];
			if(c==delimiter){break;}
			assert(Tools.isDigit(c));
			number=(number*10)+(c-'0');
		}
		return number;
	}
	
	/** Parse a gi number, or return -1 if formatted incorrectly. */
	static int parseNcbiNumber(byte[] s, char delimiter){
		if(s==null || s.length<3){return -1;}
		if(!Tools.startsWith(s, "ncbi") && !Tools.startsWith(s, ">ncbi") && !Tools.startsWith(s, "tid") && !Tools.startsWith(s, ">tid")){return -1;}
		int initial=Tools.indexOf(s, (byte)delimiter);
		if(initial<0){
			delimiter='_';
			initial=Tools.indexOf(s, (byte)delimiter);
			if(initial<0){return -1;}
		}
		if(!Tools.isDigit(s[initial+1])){return -1;}
		
		int number=0;
		for(int i=initial+1; i<s.length; i++){
			byte c=s[i];
			if(c==delimiter){break;}
			assert(Tools.isDigit(c));
			number=(number*10)+(c-'0');
		}
		return number;
	}

	public static int getID(byte[] s){return getID(s, '|');}
	/** Get the taxID from a header starting with a taxID or gi number */
	public static int getID(byte[] s, char delimiter){
		int x=parseGiNumber(s, delimiter);
		if(x>=0){return array[x];}
		return parseNcbiNumber(s, delimiter);
	}
	
	/** Get the taxID from a gi number */
	public static int getID(int gi){
		assert(gi>=0) : gi;
		assert(gi<array.length) : gi+", "+array.length;
		return array[gi];
	}
	
	public static void initialize(String fname){
		assert(fname!=null);
		if(fileString==null || !fileString.equals(fname)){
			synchronized(GiToTaxid.class){
				if(!initialized || fileString==null || !fileString.equals(fname)){
					fileString=fname;
					if(fname.contains(".int1d")){
						array=ReadWrite.read(int[].class, fname, true);
					}else{
						array=makeArray(fname);
					}
				}
				initialized=true;
			}
		}
	}
	
	public static boolean isInitialized(){return initialized;}
	
	public static synchronized void unload(){
		array=null;
		fileString=null;
		initialized=false;
	}
	
	private static int[] makeArray(String fnames){
		String[] split;
		if(new File(fnames).exists()){split=new String[] {fnames};}
		else if(fnames.indexOf(',')>=0){split=fnames.split(",");}
		else if(fnames.indexOf('#')>=0){
			assert(fnames.indexOf("/")<0) : "Note: Wildcard # only works for relative paths in present working directory.";
			File dir=new File(System.getProperty("user.dir"));
			String prefix=fnames.substring(0, fnames.indexOf('#'));
			String suffix=fnames.substring(fnames.indexOf('#')+1);
			
			File[] array=dir.listFiles();
			StringBuilder sb=new StringBuilder();
			String comma="";
			for(File f : array){
				String s=f.getName();
				if(s.startsWith(prefix) && s.startsWith(suffix)){
					sb.append(comma);
					sb.append(s);
					comma=",";
				}
			}
			split=sb.toString().split(",");
		}else{
			throw new RuntimeException("Invalid file: "+fnames);
		}
		
		IntList list=new IntList();
//		assert(max<Integer.MAX_VALUE) : "Overflow.";
//		int[] x=new int[(int)max+1];
//		Arrays.fill(x, -1);
		
		long total=0;
		for(String s : split){
			long count=addToList(s, list);
			total+=count;
		}
		return list.shrink().array;
	}
	
	private static long addToList(String fname, IntList list){
		boolean warned=false;
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		long count=0, invalid=0;
		byte[] line=bf.nextLine();
		while(line!=null){
			if(line.length>0 && Tools.isDigit(line[line.length-1])){//Invalid lines will end with tab or na
				count++;
				int tab2=Tools.indexOfNth(line, '\t', 2);
				int tab3=Tools.indexOfNth(line, '\t', 1, tab2+1);
				assert(tab2>0 && (tab2<tab3) && tab3<line.length) : tab2+", "+tab3+", "+line.length;
				assert(tab2<line.length && line[tab2]=='\t') : tab2+", "+tab3+", '"+new String(line)+"'";
				assert(tab3<line.length && line[tab3]=='\t') : tab2+", "+tab3+", '"+new String(line)+"'";
//				assert(false) : tab2+", "+tab3+", '"+new String(line)+"'";
				int tid=Parse.parseInt(line, tab2+1, tab3);
				int gi=Parse.parseInt(line, tab3+1, line.length);
				if(gi>=Shared.MAX_ARRAY_LEN || gi<0){//A gi over 2.5b was observed May 3, 2021.
					invalid++;
				}else{
				assert(gi>=0) : "tid="+tid+", gi="+gi+", line=\n'"+new String(line)+"'";
				int old=list.get(gi);
				assert(old==0 || old==tid) : "Contradictory entries for gi "+gi+": "+old+" -> "+tid+"\n'"+new String(line)+"'\ntab2="+tab2+", tab3="+tab3;
				
				list.set(gi, tid);
				
				//assert(x[gi]==-1 || x[gi]==ncbi) : "Contradictory entries for gi "+gi+": "+x[gi]+" -> "+ncbi;
//				if(x[gi]!=-1 && x[gi]!=ncbi){
//					if(!warned){
//						System.err.println("***WARNING*** For file "+fname+":\n"+
//								("Contradictory entries for gi "+gi+": mapped to both taxID "+x[gi]+" and taxID "+ncbi)+
//								"\nThis may be an error from NCBI and you may wish to report it, but it is\n"
//								+ "being suppressed because NCBI data is known to contain multi-mapped gi numbers,\n"
//								+ "at least between nucleotide and protein, and gi numbers are deprecated anyway.");
//						warned=true;
//					}
//				}else{
//					x[gi]=ncbi;
//				}
				}
			}else{
				if(line.length==0){System.err.println(fname+", "+count);}//debug
				invalid++;
			}
			line=bf.nextLine();
		}
		if(verbose){System.err.println("Count: "+count+"; \tInvalid: "+invalid);}
		bf.close();
		return count;
	}
	
//	private static int[] makeArrayOld(String fnames){
//		String[] split;
//		if(new File(fnames).exists()){split=new String[] {fnames};}
//		else{split=fnames.split(",");}
//		
//		long max=0;
//		for(String s : split){
//			max=Tools.max(max, findMaxID(s));
//		}
//		
//		assert(max<Integer.MAX_VALUE) : "Overflow.";
//		int[] x=new int[(int)max+1];
//		Arrays.fill(x, -1);
//		
//		long total=0;
//		for(String s : split){
//			long count=fillArray(s, x);
//			total+=count;
//		}
//		return x;
//	}
//	
//	private static long findMaxID(String fname){
//		ByteFile bf=ByteFile.makeByteFile(fname, true);
//		long count=0, max=0;
//		byte[] line=bf.nextLine();
//		while(line!=null){
//			count++;
//			int tab=Tools.indexOf(line, (byte)'\t');
//			long gi=Parse.parseLong(line, 0, tab);
//			max=Tools.max(max, gi);
//			line=bf.nextLine();
//		}
//		bf.close();
//		return max;
//	}
//	
//	private static long fillArray(String fname, int[] x){
//		boolean warned=false;
//		ByteFile bf=ByteFile.makeByteFile(fname, true);
//		long count=0;
//		byte[] line=bf.nextLine();
//		while(line!=null){
//			count++;
//			int tab=Tools.indexOf(line, (byte)'\t');
//			int gi=Parse.parseInt(line, 0, tab);
//			int ncbi=Parse.parseInt(line, tab+1, line.length);
//			//assert(x[gi]==-1 || x[gi]==ncbi) : "Contradictory entries for gi "+gi+": "+x[gi]+" -> "+ncbi;
//			if(x[gi]!=-1 && x[gi]!=ncbi){
//				if(!warned){
//					System.err.println("***WARNING*** For file "+fname+":\n"+
//							("Contradictory entries for gi "+gi+": mapped to both taxID "+x[gi]+" and taxID "+ncbi)+
//							"\nThis may be an error from NCBI and you may wish to report it, but it is\n"
//							+ "being suppressed because NCBI data is known to contain multi-mapped gi numbers,\n"
//							+ "at least between nucleotide and protein, and gi numbers are deprecated anyway.");
//					warned=true;
//				}
//			}else{
//				x[gi]=ncbi;
//			}
//			line=bf.nextLine();
//		}
//		if(verbose){System.err.println("Count: "+count);}
//		bf.close();
//		return count;
//	}
	
	private static int[] array;
	private static String fileString;
	
	public static boolean verbose=false;
	private static boolean initialized=false;
}
