package jasper;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;

import shared.Parser;
import shared.PreParser;
import shared.Timer;



public class DenseTreeValidate {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws FileNotFoundException, IOException {
		
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		DenseTreeValidate x=new DenseTreeValidate(args);
		
		//Run the object
		x.process(t);
		
	
	
	}
	
	/**
	 * Handles pre-parsing and parsing of user flags.
	 * Reads in the sketch similarity file (sim) and organism-parent relationship file (tree).
	 * 
	 * @param args string of the arguments input at the commandline.
	 */
	public DenseTreeValidate(String[] args) {
		
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
				
			//Handle reference variable assignment.
			}else if(a.equals("sim")){
				sim=b;
			
			//Handle kmer variable assignment.
			}else if(a.equals("tree")){
				tree=b;
					
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
		DenseTree relationshipTree=new DenseTree(tree);
		
		//Pass similarity file to create similarity matrix object
		DenseSimilarityMatrix matrix=new DenseSimilarityMatrix(sim, relationshipTree);
		
		//Add parent node similarity percentages to each node in the tree.
		//addRelationSims(relationshipTree, matrix);
		
		//Traverse the tree and add levels to all nodes.
		//Hardcoded to start at node "0" or "life" node.
		relationshipTree.beginTraverse("0");
		
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
//	void addRelationSims(DenseTree tree, DenseSimilarityMatrix matrix){
//		
//		//Iterate over organisms/nodes in the tree.
//		for ( String keyOrg : tree.keySet() ) {
//			
//			TreeNode keyNode = tree.getNode(keyOrg);
//			
//			//Identify parent node.
//			String parentName = keyNode.getParentName();
//			
//			//Get descendant nodes.
//			HashSet<String> childNames = keyNode.getChildren();
//			
//			if(!parentName.equals("0")) {
//				
//				double parSim = matrix.getSimilarity(keyOrg, parentName);
//				
//				keyNode.addParSim(parSim);
//			}
//			
//			for(String kid : childNames) {
//				
//				if(!tree.getNode(kid).parentName.equals("0")) {
//					
//					double kidSim = matrix.getSimilarity(keyOrg, kid);
//					
//					keyNode.addChildSim(kid, kidSim);
//				}
//			}
//		}
//	}
	
	void checkSimilarities(DenseTree tree, DenseSimilarityMatrix matrix) {

		//Iterate over organisms/nodes in the tree.
		for ( String keyOrgName : tree.keySet() ) {

			//If the organism isn't the life/0 node.
			if(!keyOrgName.equals("0")) {
				//System.out.println("key org " + keyOrg);

				//Get the node from the tree
				TreeNode keyNode = tree.getNode(keyOrgName);
				
				int keyNodeID = keyNode.getNodeId();

				tree.root.resetIdentity();
				
				tree.setIdentity(keyNode, matrix);
				
				tree.root.percolateIdentityUp(keyNode.nodeId);
				
				//Prevents analysis of "empty" nodes that don't contain sequences (genus/phylum/etc).
				//TODO: if there are no sibling nodes, the parent sim could be 0.
				if(keyNode.parentNode.averageIdentity() != 0.0) {

					//Identify parent node.
					String parentName = keyNode.getParentName();
					
					double parentSimilarity = matrix.getSimilarity(keyOrgName, parentName);

					//Get the child node names.
					HashSet<String> childNameSet = keyNode.getChildren();

					//Iterate over the organisms in the matrix.
					for(String matrixOrg : tree.keySet()) {

						TreeNode matrixOrgNode = tree.getNode(matrixOrg);

						//if we aren't comparing similarities of the node to itself and
						//if we aren't examining a child node and
						//if we aren't examining a parent node
						if(!matrixOrgNode.isDescendantOf(keyNode) && !matrixOrgNode.isAncestorOf(keyNode) && !matrixOrg.equals(parentName) ) {

							double matrixOrgSim = matrix.getSimilarity(keyOrgName, matrixOrg);


							if(matrixOrgSim > parentSimilarity) {


								System.out.println();
								System.out.println("problem");
								System.out.println("key org " + keyOrgName);
								//System.out.println("kid name " + minChildName);
								System.out.println("par name " + parentName);
								System.out.println("other org " + matrixOrg);
								System.out.println("par sim " + parentSimilarity);
								//System.out.println("child sim " + minChildSim);
								System.out.println("matrix sim " + matrixOrgSim);

								keyNode.flagRelation(matrixOrg, matrixOrgSim);
								System.out.println(keyNode.getFlaggedRelations());
							}

						}
					}
				}
			}
		}
	}

	
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private String sim=null;
	private String tree=null;
	private String out=null;
	
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
