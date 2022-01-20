package jasper;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;

import shared.Tools;


public class NCBISparseTree {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Takes in an input file with 2 columns (organism, parent organism) and adds these to TreeNodes
	 * TreeNodes are then added to a HashMap and children node values are added to each node if applicable
	 * 
	 * @param inputFile The input file you wish to have values added to the Tree object
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public NCBISparseTree(String inputFile) throws FileNotFoundException, IOException {
		
		int nodeId = 0;
		in = inputFile;
		
		//parse file. create each node and place in 
		try (BufferedReader br = new BufferedReader(new FileReader(in))) {
	        String line;
	        
	        while ((line = br.readLine()) != null) {
	        	
	        	
	        	//if line is the header line, split and assign to variable.
	        	if(line.startsWith("#")) {header=line.split("\t");
	        	} else {
	        		//String[] data = line.split("\t");
	        		String[] data = Tools.tabPattern.split(line);

	        		//Make sure you're not adding the header line to any data structure
	        		//if(!Arrays.asList(header).contains(data[0])) {

	        		//Create a TreeNode containing the name of the organism and the parent node/organism
	        		//System.out.println(data[0]);

	        		int taxID = Integer.valueOf(data[0]);
	        		int parentTaxID = Integer.valueOf(data[2]);
	        		String taxonomicRank = data[4];

	        		NCBITreeNode orgNode = new NCBITreeNode(taxID, data[0], parentTaxID, nodeId, taxonomicRank);

	        		if(nodeId == 0) {root = orgNode;}

	        		nodeId++;

	        		//Add node to HashMap nodes with the name of the organism as the key
	        		nodeMap.put(taxID, orgNode);

	        		nodeList.add(orgNode);

	        		//Add line to lines list for further processing
	        		lines.add(line);

	        		//Increment linesProcessed
	        		linesProcessed++;
	        			
	        		
	        		//}
	        	}
	        }
		}
	
		//Run method to add children nodes to each node if applicable
		addChild(nodeMap, lines);
		
		//return nodes;
	}
	
	
	/**
	 * Adds children node names to each node if applicable
	 * 
	 * @param treeNodeMap HashMap of TreeNode objects
	 * @param lineList ArrayList<String> of lines from the input file
	 */
	void addChild(HashMap<Integer, NCBITreeNode> treeNodeMap, ArrayList<String> lineList) {
		int par;
		int org;
		String rank;
		
		//iterate over lines from the file and split into the organism and the parent node
		for(String line : lineList) {
			//String[] split=line.split("\t");
			String[] split = Tools.tabPattern.split(line);
			
			//isolate the organism and the parent from the split
			org = Integer.valueOf(split[0]);
			par = Integer.valueOf(split[2]);
			rank = split[4];
			
			//get the organism node and parent node
			NCBITreeNode orgNode = treeNodeMap.get(org);
			NCBITreeNode parNode = treeNodeMap.get(par);
			
			//Assert parent node isn't empty or parent node is the 0/life node.
			assert(parNode != null || par == 1);
			
			//Assert the query organism node isn't empty, if it is, return node name.
			assert(orgNode != null): org;
			
			//Add the child node name to the query node.
			parNode.addChildren(org);
			
			//add query node to its parent node's list of children nodes.
			parNode.childNodes.add(orgNode);
			
			orgNode.parentNode = parNode;
			
		}
			
	}
	
	/**
	 * Returns a StringBuilder of names of organisms/nodes along with
	 * the parent node and the names of children nodes.
	 * 
	 * @return StringBuilder
	 */
	public String toString() {
		StringBuilder sb=new StringBuilder();
		for(Entry<Integer, NCBITreeNode> e : nodeMap.entrySet()) {
			
			NCBITreeNode tn = e.getValue();
	    	
	    	sb.append(tn);
	    	sb.append('\n');
	    }
		return sb.toString();
	}
	
	/**
	 * Returns Set<String> of node keys for the tree.
	 * @return Set<String>
	 */
	public Set<Integer> keySet() {
		return nodeMap.keySet();
	}
	
	/**
	 * Starting point for adding levels to nodes in the tree.
	 * @param nodeID_ Lowest node name, corresponding to "Life"
	 */
	public void beginTraverse(int nodeID_) {
		NCBITreeNode firstNode = nodeMap.get(nodeID_);
		firstNode.assignLevels(1, "life");
	}
	
	/**
	 * Returns TreeNode from tree based on String node name.
	 * @param TaxonID String name of node (organism/file)
	 * @return TreeNode.
	 */
	public NCBITreeNode getNodeByTaxID(int TaxonID) {
		return nodeMap.get(TaxonID);
	}
	
	/**
	 * Returns TreeNode from the tree based on the node ID
	 * @param nodeID int nodeId
	 * @return TreeNode.
	 */
	public NCBITreeNode getNodeByNodeID(int nodeID) {
		return nodeList.get(nodeID);
	}
	
	/**
	 * Currently takes node and adds all descendant node names to a HashSet<String>.
	 * Returned with the getDescendentNames method.
	 * @param nodeName
	 */
	public void beginAddDescendants(int taxID_) {
		
		//Place target node in variable
		NCBITreeNode earliestNode = nodeMap.get(taxID_);
		
		//Run the method to add descendant names to HashSet in each node.
		earliestNode.nodeAddDescendantNames(nodeMap.get(taxID_).descendentIDs);
	}
	
	
	public void assignMatrixIdentity(NCBISparseSimilarityMatrix matrix, TreeNode node) {
		
	}
	
	/**
	 * Returns boolean of whether the String organism name/node name is found in the tree.
	 * @param orgName String organism/node name.
	 * @return boolean
	 */
	public boolean containsTaxID(int orgTaxID) {
		return nodeMap.containsKey(orgTaxID);
	}
	
//	//TODO: make a more efficient method of getting the total node count.
//	public int getOrgCount() {
//		int max = 0;
//		for(String node : nodeMap.keySet()) {
//			int id = this.getNode(node).getNodeId();
//			if(id > max) {
//				max = id;
//			}
//		}
//		return max;
//	}
	
	/**
	 * Returns the size of the HashMap containing all nodes in the tree.
	 * @return int TreeNode count.
	 */
	public int getOrgCount() {
		
		//Return the size of the Set of keys in the nodeMap HashMap.
		return nodeList.size();
	}
	
	/**
	 * Sets the identity of all other TreeNodes in relation to the input TreeNode name.
	 * @param keyNode TreeNode query node.
	 * @param matrix SparseSimilarityMatrix name containing similarity Comparison objects.
	 */
	public void setIdentity(NCBITreeNode keyNode, NCBISparseSimilarityMatrix matrix) {
		
		//Get the row containing all Comparisons for the query node.
		ArrayList<NCBIComparison> row = matrix.getOrgRowByTaxonID(keyNode.taxID);
		
		int votes = NCBISparseTreeValidate.MAX_VOTES;
		
		//Iterate over the row.
		for(int i=0; i<row.size(); i++) {
			
			//Get a Comparison object from the row.
			NCBIComparison c = row.get(i);
			
			//Get the other nodes ID from the Comparison object.
			int otherNodeId = c.refNodeID;
			
			//Get the TreeNode of the node being compared to the query node.
			NCBITreeNode otherNode = nodeList.get(otherNodeId);
			
			//Set the other nodes identity to the value in the comparison.
			otherNode.identity = c.identity;
			
			otherNode.votes = votes;
			votes = Tools.max(votes - 1, 0);
			
		}
	}
	
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Root node of the tree. Usually the "life" node or node 1 in NCBI.
	 */
	NCBITreeNode root;
	
	/**
	 * HashMap holding the taxon ID of the organisms as keys and the organism node as values.
	 * 
	 */
	HashMap<Integer, NCBITreeNode> nodeMap = new HashMap<Integer, NCBITreeNode>();
	
	/**
	 * ArrayList of NCBITreeNode objects.
	 */
	ArrayList<NCBITreeNode> nodeList = new ArrayList<NCBITreeNode>();
		
	/**
	 * ArrayList of all lines in input file. Need these to fill in values for children nodes
	 */
	ArrayList<String> lines = new ArrayList<String>();
	
	/**
	 * Header line of input file.
	 */
	private String[] header;
	
	/**
	 * Input file name.
	 */
	private String in=null;
	
	/**
	 * Number of lines processed for data from input file.
	 */
	private long linesProcessed=0;
	
	/**
	 * Node level counter.
	 */
	int orgLvl = 0;
	
	//int orgCount = 0;
}
