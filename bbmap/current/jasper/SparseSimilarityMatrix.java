package jasper;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class SparseSimilarityMatrix {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Takes in a file of sketch similarity percentages from SketchCompare.
	 * Returns a sparse matrix object containing each percentage
	 * 
	 * @param inputFile The file containing pairwise comparisons of each sketch
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public SparseSimilarityMatrix(String inputFile, SparseTree tree_) throws FileNotFoundException, IOException {

		//Assigns the input tree object to the tree variable.
		tree = tree_;

		//Take file name as input for building tree of related nodes
		in = inputFile;

		//Read in file, add header line and add to header variable
		try (BufferedReader br = new BufferedReader(new FileReader(in))) {
			String line;

			//while line isn't empty, process
			while ((line = br.readLine()) != null) {

				//if line is the header line, split and assign to variable.
				//may be used when header becomes more complex
				if(line.startsWith("#")) {header=line.split("\t");
				} else {

					//If not a header line, split on tab.
					String[] data = line.split("\t");

					//Query organism is column 0.
					String queryName = data[0];
					//String refName = data[1];

				}
			}
		}
		
		//Get the total number of organisms in the tree.
		orgCount = tree.getOrgCount();
		
		//Initialize the matrix with the appropriate size of all nodes.
		sparseMatrix = new ArrayList[orgCount + 1];

		//Iterate over the matrix and add an ArrayList<Comparison> to each ArrayList.
		for(int i=0; i<sparseMatrix.length; i++) {
			
			sparseMatrix[i] = new ArrayList<Comparison>();
			
		}
		
		//Begin reading the file a second time.
		try (BufferedReader br = new BufferedReader(new FileReader(in))) {
			String line;

			//while line isn't empty, process
			while ((line = br.readLine()) != null) {

				//If line is the header line, split and assign to variable.
				//may be used when header becomes more complex
				if(line.startsWith("#")) {assert true;
				} else {
					
					//If not a header line, split on tab.
					String[] data = line.split("\t");
					
					//Column 0 is query name.
					String queryName = data[0];
					
					//Column 1 is reference name.
					String refName = data[1];
					
					//Column 2 is the similarity percentage.
					double similarity = Double.parseDouble(data[2]);
					
					//Check that both names are in the HashMap (too slow?)
					if(tree.containsName(queryName)==true && tree.containsName(refName)) {
						
						//Get the positions assigned to both organisms.
						int queryPos = nameToNodeId(queryName);
						int refPos = nameToNodeId(refName);
						
						Comparison currentComparison = new Comparison(queryPos, refPos, similarity);
						
						//Add the similarity percentage to the appropriate matrix position.
						sparseMatrix[queryPos].add(currentComparison);
					}
				}
			}
		}	
	}
	
	/**
	 * Method for taking the node name and returning the node ID value
	 * @param orgName the organism node name (String).
	 * @return int The node ID of the organism name taken as input.
	 */
	public int nameToNodeId(String orgName) {
		
		//Get the node associated with the input name.
		TreeNode org = tree.getNode(orgName);
		
		//Asserts the org nod is in the tree.
		assert(org != null) : orgName;
		
		//Return the int node ID.
		return org.nodeId;
	}
	
	
	/**
	 * Prints out the entire matrix.
	 * Impractical in cases of large input datasets.
	 * 
	 */
	public String toString() {
		StringBuilder sb=new StringBuilder();
		for (int i = 0; i < sparseMatrix.length; i++) {
		    for (int j = 0; j < sparseMatrix[i].size(); j++) {
		        sb.append(sparseMatrix[i].get(j) + " ");
		    }
		    sb.append('\n');
		}
		return sb.toString();
	}
	
	
//TODO: This method is slow and doesnt work, need something better.	
//	/**
//	 * Returns the similarity of two specified organisms.
//	 * Both organisms must have been compared using SketchCompare.
//	 * 
//	 * @param org1 The Name of an organism.
//	 * @param org2 The name of a second organism.
//	 * @return similarity The Double percentage similarity between the two sketches.
//	 */
//	public Comparison getComparison(String org1, String org2) {
//		int orgName1 = nameToNodeId(org1);
//		int orgName2 = nameToNodeId(org2);
//		
//		return sparseMatrix[orgName1].get(orgName2);
//	}

	
	public int getSize() {
		return orgCount;
	}
	
	
	public ArrayList<Comparison> getOrgRow(String orgName) {
		int rowNum = tree.nodeMap.get(orgName).nodeId;
		return sparseMatrix[rowNum];
	}
	

	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
    
	/**
	 * A SparseTree object that contains taxonomic information relevant to this matrix.
	 */
	final SparseTree tree;
	
	/**
	 * An arraylist containing comparisons between nodes in the tree.
	 */
	private final ArrayList<Comparison>[] sparseMatrix;
	
	/**
	 * The number of sketches being analyzed.
	 */
	private int orgCount;
	
	/**
	 * ArrayList that will hold the lines of the input file.
	 */
	ArrayList<String> lines = new ArrayList<String>();
	
	/**
	 * Header line of the comparison input file.
	 */
	private String[] header;
	
	/**
	 * Input file name.
	 */
	private String in=null;
	
	/**
	 * Number of lines processed from the sketch comparison file.
	 */
	private long linesProcessed=0;
	
}
