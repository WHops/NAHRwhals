package jasper;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;


public class DenseTree {
	
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
	public DenseTree(String inputFile) throws FileNotFoundException, IOException {
		
		int nodeId = 0;
		in = inputFile;
		
		//parse file. create each node and place in 
		try (BufferedReader br = new BufferedReader(new FileReader(in))) {
	        String line;
	        
	        while ((line = br.readLine()) != null) {
	        	
	        	
	        	//if line is the header line, split and assign to variable.
	        	if(line.startsWith("#")) {header=line.split("\t");
	        	} else {
	        		String[] data = line.split("\t");
	        		
	        		//Make sure you're not adding the header line to any data structure
	        		if(!Arrays.asList(header).contains(data[0])) {
	        			
	        			//Create a TreeNode containing the name of the organism and the parent node/organism
	        			//System.out.println(data[0]);
	        			TreeNode orgNode = new TreeNode(data[0], data[1], nodeId);
	        			
	        			if(nodeId == 0) {root = orgNode;}
	        			
	        			nodeId++;
	        			
	        			//Add node to HashMap nodes with the name of the organism as the key
	        			nodeMap.put(data[0], orgNode);
	        			
	        			nodeList.add(orgNode);
	        			
	        			//Add line to lines list for further processing
	        			lines.add(line);
	        			
	        			//Increment linesProcessed
	        			linesProcessed++;
	        		}
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
	void addChild(HashMap<String, TreeNode> treeNodeMap, ArrayList<String> lineList) {
		String par;
		String org;
		
		//iterate over lines from the file and split into the organism and the parent node
		for(String line : lineList) {
			String[] split=line.split("\t");
			
			//isolate the organism and the parent from the split
			org = split[0];
			par = split[1];
			
			//get the organism node and parent node
			TreeNode orgNode = treeNodeMap.get(org);
			TreeNode parNode = treeNodeMap.get(par);
			
			//Assert parent node isn't empty or parent node is the 0/life node.
			assert(parNode != null || par.equals("0"));
			
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
		for(Entry<String, TreeNode> e : nodeMap.entrySet()) {
			
			TreeNode tn = e.getValue();
	    	
	    	sb.append(tn);
	    	sb.append('\n');
	    }
		return sb.toString();
	}
	
	/**
	 * Returns Set<String> of node keys for the tree.
	 * @return Set<String>
	 */
	public Set<String> keySet() {
		return nodeMap.keySet();
	}
	
	/**
	 * Starting point for adding levels to nodes in the tree.
	 * @param nodeName Lowest node name, corresponding to "Life"
	 */
	public void beginTraverse(String nodeName) {
		TreeNode firstNode = nodeMap.get(nodeName);
		firstNode.traverse(0);
	}
	
	/**
	 * Returns node from tree.
	 * @param nodeName
	 * @return
	 */
	public TreeNode getNode(String nodeName) {
		return nodeMap.get(nodeName);
	}
	
	public TreeNode getNode(int nodeId) {
		return nodeList.get(nodeId);
	}
	
	/**
	 * Currently takes node and adds all descendant node names to a HashSet<String>.
	 * Returned with the getDescendentNames method.
	 * @param nodeName
	 */
	public void beginAddDescendants(String nodeName) {
		//Place target node in variable
		TreeNode earliestNode = nodeMap.get(nodeName);
		
		//Run the method to add descendant names to HashSet in each node.
		earliestNode.nodeAddDescendantNames(nodeMap.get(nodeName).descendentNames);
	}
	
	
	public void assignMatrixIdentity(DenseSimilarityMatrix matrix, TreeNode node) {
		
	}
	
	public boolean containsName(String orgName) {
		return nodeMap.containsKey(orgName);
	}
	
	public int getOrgCount() {
		int max = 0;
		for(String node : nodeMap.keySet()) {
			int id = this.getNode(node).getNodeId();
			if(id > max) {
				max = id;
			}
		}
		return max;
	}
	
	public void setIdentity(TreeNode node, DenseSimilarityMatrix matrix) {
		
		double[] row = matrix.getOrgRow(node.orgName);
		for(int i=0; i<row.length; i++) {
			double id = row[i];
			nodeList.get(i).identity = id;
		}
	}
	
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	TreeNode root;
	
	//HashMap holding the names of the organisms as keys and the organism node as values
	HashMap<String, TreeNode> nodeMap = new HashMap<String, TreeNode>();
	
	ArrayList<TreeNode> nodeList = new ArrayList<TreeNode>();
		
	//ArrayList of all lines in input file. need these later to fill in values for children nodes
	ArrayList<String> lines = new ArrayList<String>();
	
	//Header line of input file
	private String[] header;
	
	//Input file name
	private String in=null;
	
	//Number of lines processed for data from input file
	private long linesProcessed=0;
	
	//Node level counter
	int orgLvl = 0;
	
	int orgCount = 0;
}
