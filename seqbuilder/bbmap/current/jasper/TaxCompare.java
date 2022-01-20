package jasper;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import shared.Timer;


public class TaxCompare {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws FileNotFoundException, IOException {
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		TaxCompare x=new TaxCompare(args);
		
		x.process(t);
		
		//Organism cat = new Organism(9, "snake");
		//cat.printOrg();
		
		}

	/**
	 * Constructor. Handle inputs.
	 * @param args Command line arguments
	 */
	public TaxCompare(String[] args){
		
		String arg=args[0];
		String[] split=arg.split("=");
		String a=split[0].toLowerCase();
		String b=split.length>1 ? split[1] : null;
		if(b!=null && b.equalsIgnoreCase("null")){b=null;}
		in = b;
	}
	
	
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	/**
	 * Process input sketch comparison file, pass names of organisms (eventually taxon IDs) to matrix.
	 * 
	 * @param t
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	void process(Timer t) throws FileNotFoundException, IOException{
		
		//Set for taxon names.
		Set<String> nameSet = new HashSet<String>();
		
		//ArrayList of all names in sketch comparison file, not unique values.
		ArrayList<String> names = new ArrayList<String>();
		
		//ArrayList of all lines. need these later to fill in values for the matrix
		//This will get ugly with large comparisons.
		ArrayList<String> lines = new ArrayList<String>();
		
		//Read in file, add header line and add to header variable
	    try (BufferedReader br = new BufferedReader(new FileReader(in))) {
	        String line;
	        
	        while ((line = br.readLine()) != null) {
	        	
	        	//if line is the header line, split and assign to variable.
	        	if(line.startsWith("#")) {header=line.split("\t");
	        	} else {
	        		String[] data = line.split("\t");
	        		
	        		//make sure the data in column 1 isnt in the header line
		        	//column 1 should be query names
	        		if(!Arrays.asList(header).contains(data[0])) {nameSet.add(data[0]);}
	        		lines.add(line);
	        	
	        	}
	        	
	        }
	    
	    //create matrix and add values
	    double[][] matrix = new double[nameSet.size() + 1 ][nameSet.size() + 1];
	    
	    //loop over lines and fill in matrix
	    for(int i=0; i<lines.size(); i++) {
	    	
	    	fillMatrix(matrix, nameSet, lines.toArray()[i]);
	    }
	    
	    System.out.println(nameSet);
	    for (int i = 0; i < matrix.length; i++) {
		    for (int j = 0; j < matrix[i].length; j++) {
		        System.out.print(matrix[i][j] + " ");
		    }
		    System.out.println();
		}
	    
	    
	    }
	    t.stop();
	    System.out.println(t);
	    
	}
	
	/**
	 * Fill matrix with relationship information of organisms output from sketch comparison.
	 * 
	 * @param matrix Matrix of comparison percentage values.
	 * @param setNames Set of names of included organisms.
	 * @param object Line of sketch comparison output file.
	 */
	void fillMatrix(double[][] matrix, Set<String> setNames, Object object) {
		//System.out.println(object);
		//cast line as string
		String stringLine = (String) object;
		
		//split line
		String[] lineData = stringLine.split("\t");
		
		//place both organism names in variables
		//qName is the query, column 1
		String qName = lineData[0];
		String mName = lineData[1];
		double similarity = Double.parseDouble(lineData[2]);
		
		//set positions variables
		int qPos = -1;
		int mPos = -1;
		
		//convert setNames to array that can be iterated over
		String[] nameArray = setNames.toArray(new String[setNames.size()]);
		
		//loop over setNames and get each organisms positions within the matrix
		for(int i = 0; i<nameArray.length; i++) {
			if(nameArray[i].contentEquals(qName)) {qPos = i;}
			else if(nameArray[i].contentEquals(mName)) {mPos = i;}
			
			
			//after finding both name positions, add similarity value to matrix
			if(qPos!=-1 && mPos!=-1) {matrix[qPos][mPos] = similarity;}
		}
		
		//displays matrix
		/*
		for (int i = 0; i < matrix.length; i++) {
		    for (int j = 0; j < matrix[i].length; j++) {
		        System.out.print(matrix[i][j] + " ");
		    }
		    System.out.println();
		}*/
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in=null;
	private String out=null;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private String[] header;
	//private final FileFormat ffin;
	//private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	//private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
}
