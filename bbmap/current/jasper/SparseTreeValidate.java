package jasper;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import shared.Parser;
import shared.PreParser;
import shared.Timer;



public class SparseTreeValidate {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws FileNotFoundException, IOException {
		
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		SparseTreeValidate x=new SparseTreeValidate(args);
		
		//Run the object
		x.process(t);
		
	
	
	}
	
	/**
	 * Handles pre-parsing and parsing of user flags.
	 * Reads in the sketch similarity file (sim) and organism-parent relationship file (tree).
	 * 
	 * @param args string of the arguments input at the commandline.
	 */
	public SparseTreeValidate(String[] args) {
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Primary parsing of standard arguments found in all bbmap programs (maxReads, parseSam, parseZip, etc).
		Parser parser=new Parser();
		
		//Loop through arguments up to the maximum number of arguments input.
		//process all remaining arguments. 
		for(int i=0; i<args.length; i++){
			
			//Grab argument string at index.
			String arg=args[i];
			
			//Split argument string on "=".
			String[] split=arg.split("=");
			
			//Convert the left side to lowercase.
			String a=split[0].toLowerCase();
			
			//Ternary conditional statement: is the length of the split greater than 1 (thus, an actual input)?
			//if so, the right side of the split is the b variable, if not, b is null.
			String b=split.length>1 ? split[1] : null;
			
			//If b isn't null but a string "null" was input, convert b to null.
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			//Unused example statement. does nothing currently. start here for adding new flag parsing.
			if(a.equals("parse_flag_goes_here")){
				
			//Handle similarities file variable assignment.
			}else if(a.equals("sim")){
				sim=b;
			
			//Handle tree file variable assignment.
			}else if(a.equals("tree")){
				tree=b;
				
			//Tells program to write tree graphs in .dot format.
			}else if(a.equals("writetrees")) {
				writeTrees = true;
	
			//If writetrees is true and the path to an output dir is included
			//set variable to the path to the desired directory
			}else if(a.equals("outpath") && writeTrees == true) {
				outpath=b;
				if(outpath.endsWith("/")) {System.out.println("output path correct");}
				else {outpath = outpath + "/";}
					
			//Parses in and out flags, handles all flags not recognized earlier in class.
			}else if(parser.parse(arg, a, b)){
				
			//If not one of the known parameters, let the user know they made a mistake.
			}else{
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
			
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Creates the similarity matrix and the relationship tree.
	 * Also passes these objects to the inner processing method for analysis.
	 * 
	 * @param t
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	void process(Timer t) throws FileNotFoundException, IOException{
		
		//Pass input file to Tree class to create tree
		SparseTree relationshipTree=new SparseTree(tree);
		
		//Pass similarity file to create similarity matrix object
		SparseSimilarityMatrix matrix=new SparseSimilarityMatrix(sim, relationshipTree);
		
		//Add parent node similarity percentages to each node in the tree.
		//addRelationSims(relationshipTree, matrix);
		
		System.out.println(relationshipTree.getNode(5));
		
		//Traverse the tree and add levels to all nodes.
		//Hardcoded to start at node "0" or "life" node.
		relationshipTree.beginTraverse("0");
		
		//Sets the identities by beginning at a particular node and working backwards.
		//relationshipTree.setIdentity(relationshipTree.getNode(10), matrix);
		
		//TODO: remove 10 since this is just part of testing/percolate
		//relationshipTree.root.percolateIdentityUp(10);
		
		//Check similarities.
		checkSimilarities(relationshipTree, matrix);
		
		//System.out.println(relationshipTree);
		
		
		t.stop();
		
		outstream.println("Time:                         \t"+t);
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
//	/**
//	 * Iterate over nodes in the Tree and compare present relationships
//	 * with values found in the similarity matrix.
//	 * 
//	 * @param tree Tree object containing TreeNode objects detailing the parent and children of each node.
//	 * @param matrix SimilarityMatrix2 object containing percentage similarity of sketches.
//	 */
//	void addRelationSims(SparseTree tree, SparseSimilarityMatrix matrix){
//
//		//Iterate over organisms/nodes in the tree.
//		for ( String keyOrg : tree.keySet() ) {
//
//			//add the node of current focus to a variable in the loop.
//			TreeNode keyNode = tree.getNode(keyOrg);
//
//			//Get Parent node of keyNode.
//			TreeNode parentNode = keyNode.parentNode;
//
//			//If statement to ignore the root/"life" node when obtaining similarities.
//			//Then get the parent nodes similarity to the key node.
//			//add the similarity value to the key node.
//			if(!parentNode.orgName.equals("0")) {
//
//				keyNode.addParSim(parentNode.averageIdentity());
//				System.out.println(parentNode.averageIdentity());
//			}
//		}
//	}
	
	
	/**
	 * Method to check tree for surprising similarity values for each node in the tree.
	 * Ignores the root/"life" node, any node without any sequence.
	 * Any node with a higher similarity to another node than to its parent that
	 * isn't an descendant is flagged as possibly erroneous.
	 * 
	 * @param tree The tree object storing taxon nodes (TreeNode objects).
	 * @param matrix The sparse similarity matrix containing all possible pairwise similarity values.
	 */
	void checkSimilarities(SparseTree tree, SparseSimilarityMatrix matrix) {

		//Iterate over organisms/nodes in the tree.
		for ( String keyOrg : tree.keySet() ) {

			//If the organism isn't the life/0 node.
			if(!keyOrg.equals("0")) {

				//Get the node from the tree
				TreeNode keyNode = tree.getNode(keyOrg);

				//Reset the identity values of the TreeNodes
				//This is done so all identity values are relative to the current keyNode.
				tree.root.resetIdentity();
				
				//Set the identities of the tree based on the current keyNode.
				tree.setIdentity(keyNode, matrix);
				
				//Percolate the average identities of each node upwards through the tree
				//relative to the current keyNode.
				tree.root.percolateIdentityUp(keyNode.nodeId);
				
				//System.out.println(tree.root.toDot());
				
				//Handle writing tree relationships in .dot format.
				//If selected by user.
				if(writeTrees == true) {
					
					//Creates the tree output file name along with the path to the desired directory.
					String currentTreeName = outpath + "SimilarityTree_" + keyNode.orgName + "_.dot";
					
					//Create the output file.
					createFile(currentTreeName);

					//Write the .dot formatted tree to the file.
					writeToFile(currentTreeName, tree.root.toDot());

					//assert false;
				}
				//Prevents analysis of "empty" nodes that don't contain sequences (genus/phylum/etc).
				//TODO: if there are no sibling nodes, the parent sim could be 0.
				if(keyNode.parentNode.averageIdentity() != 0.0) {

					//Identify parent node.
					String parentName = keyNode.getParentName();
					
					//Get parent node ID
					//int parentID = tree.getNode(parentName).getNodeId();

					//Get the row of similarity values associated with
					//the key node and each other node.
					ArrayList<Comparison> keyOrgRow = matrix.getOrgRow(keyOrg);

					//Iterate over the node organism names.
					for(Comparison rowOrgComparison : keyOrgRow) {
						
						//Get the node being iterated over in the tree.
						TreeNode matrixOrgNode = tree.getNode(rowOrgComparison.refID);

						//if we aren't comparing similarities of the node to itself and
						//if we aren't examining a child node and
						//if we aren't examining a parent node
						if(!matrixOrgNode.isDescendantOf(keyNode.parentNode) && 
								!matrixOrgNode.isAncestorOf(keyNode)) {

							//Get similarity between the key node and any pairwise compared node
							double matrixOrgSim = rowOrgComparison.identity;

							//If the similarity value is higher than the similarity
							//between the key node and its parent.
							if(matrixOrgSim > keyNode.parentSimilarity()) {

								//Currently prints out a bunch of node/similarity info
								//TODO: fix flagged node handling
								System.out.println();
								System.out.println("problem");
								System.out.println("key org " + keyOrg);

								System.out.println("par name " + parentName);
								System.out.println("other org " + rowOrgComparison);
								System.out.println("par sim " + keyNode.parentSimilarity());

								System.out.println("matrix sim " + matrixOrgSim);

								//keyNode.flagRelation(rowOrgComparison, matrixOrgSim);
								System.out.println(keyNode.getFlaggedRelations());
							}

						}
					}
				}
			}
		}
	}

	/**
	 * Method for creating a file with an input file name.
	 * @param fileName String for the file name.
	 */
	void createFile(String fileName) {
		
		try {
			
			//Creates a file object using the input file name.
			File myObj = new File(fileName);
			
			//If the file is successfully created, print a message saying so.
			if (myObj.createNewFile()) {
				System.out.println("File created: " + myObj.getName());
			
			//If the file cannot be created, the file already exists
			} else {
				System.out.println("File already exists.");
			}
			
		//Catch a message if their are exceptions.
		} catch (IOException e) {
			System.out.println("An error occurred.");
			e.printStackTrace();
		}
	}

	/**
	 * Method to write information to a file.
	 * @param fileName Name of the file created using the createFile method.
	 * @param sb StringBuilder containing information to be written to the file.
	 */
	void writeToFile(String fileName, StringBuilder sb) {
		
		try {
			
			//Create a FileWriter object using the name of the file.
			FileWriter myWriter = new FileWriter(fileName);
			
			//Write the StringBuilder contents to the file
			myWriter.write(sb.toString());
			
			//Close the FileWriter object.
			myWriter.close();
			
			//Message on success.
			System.out.println("Successfully wrote to the file.");
		
		//Catch exceptions and write a message to alert the user.
		} catch (IOException e) {
			System.out.println("An error occurred.");
			e.printStackTrace();
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Pairwise similarity file generated by BBSketch.
	 */
	private String sim=null;
	
	/**
	 * Tree file containing taxonomic relationships.
	 */
	private String tree=null;
	
	/**
	 * Output file name.
	 */
	private String out=null;
	
	/**
	 * Output path for any output files.
	 */
	private String outpath=null;
	
	/**
	 * Flag for indicating whether to write tree .dot files.
	 * If present in commandline args, set to true.
	 */
	private boolean writeTrees = false;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	private long taxa=0;
	
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
	/**Output stream that output statistics are piped through to the output file. */
	private java.io.PrintStream outstream=System.err;
	
	
}
