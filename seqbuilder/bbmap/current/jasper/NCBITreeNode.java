package jasper;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class NCBITreeNode {
	
	/**
	 * Class of object NCBITreeNode. Contains the taxon ID number, name, parent node, node ID within the tree
	 * and the Taxonomic rank of the node.
	 * @param taxID int Taxon ID.
	 * @param name String Name of node.
	 * @param parentNode_ int Parent node ID.
	 * @param nodeId int Node ID within the tree, separate from taxon ID.
	 * @param taxonomicRank_ String Corresponds to something like genus, phylum, kingdom, etc.
	 */
	public NCBITreeNode(int taxID, String name, int parentNode_, int nodeId, String taxonomicRank_) {
	    //this.taxId = id;
	    this.taxID = taxID;
	    //this.children = kids;
	    this.parentID = parentNode_;
	    this.nodeId = nodeId;
	    this.orgName = name;
	    this.taxonomicRank = taxonomicRank_;
	    
    }
	
	/**
	 * Add a child node to the HashSet of children nodes for this node.
	 * 
	 * @param org Name of child node/organism.
	 */
	public void addChildren(int org) {
	      childIDs.add(org);
	   }
	
	/**
	 * Add a descendant node to the HashSet of descendant nodes for this node.
	 * 
	 * @param kid Name of child node/organism.
	 */
	public void nodeAddDescendantNames(HashSet<Integer> desNames) {
		//Iterate over child nodes
		for(NCBITreeNode childNode : childNodes) {
			
			//If the name of the descendant node is not equal to this node
			//This is only applicable for the "life"/"0" node
			if(childNode.taxID != taxID) {
				
				//Add the descendant node name to the input HashSet of names.
				desNames.add(childNode.taxID);
				
				//Run this function on the child nodes of each child node of the starting node.
				childNode.nodeAddDescendantNames(desNames);
			}
		}
	}
	
	/**
	 * Return HashSet containing the names of descendant nodes.
	 * @return descendantNames HashSet<String> of descendant nodes.
	 */
	public HashSet getDescendantNames() {
		return descendentIDs;
	}
	
	/**
	 * Returns HashSet of child nodes for this node.
	 * 
	 * @return children HashSet of child nodes.
	 */
	public HashSet<Integer> getChildren() {
	      return childIDs;
	   }
	
	/**
	 * Returns the parent node of this node.
	 * 
	 * @return parent The name of the parent node.
	 */
	public int getParentTaxonID() {
	      return parentID;
	   }
	
	/**
	 * Returns a string of the structure <Organism name>, <Parent Organism/node name>, <Child node names if any>.
	 */
	public String toString() {
		return "Taxon ID = " + taxID + ", Parent ID = " + parentID + ", Child IDs = " + childIDs + ", Level = " + level
				+ ", Nodes with identity = " + nodesWithIdentity +  ", votes = "+ votes + ", Identity = " + identity + 
				", Average identity = " + averageIdentity();
	}
	
//	/**
//	 * Add similarity values and child node names to a HashMap.
//	 * 
//	 * @param childName Name of direct child node.
//	 * @param similarity Similarity percentage between child node and this node.
//	 */
//	public void addChildSim(String childName, double similarity) {
//		childSims.put(childName, similarity);
//	}
	
//	/**
//	 * Returns the similarity percentage between this node and a child node.
//	 * 
//	 * @param childName The name of the child node.
//	 * @return double similarity percentage
//	 */
//	public double getChildSim(String childName) {
//		double sim = childSims.get(childName);
//		return sim;
//	}
	
//	/**
//	 * Takes a double similarity percentage and associates that with the parent node.
//	 * 
//	 * @param similarity Percentage similarity between node and parent of type double
//	 */
//	public void addParSim(double similarity) {
//		parSim = similarity;
//	}
	
//	/**
//	 * Returns the minimum similarity percentage of all descendants.
//	 * 
//	 * @return Similarity between node and child node with lowest similarity.
//	 */
//	public double minimumDescendantSim() {
//
//		for(String childName : childNames) {
//			if(childSims.get(childName) < minChildSim) {
//				minChildSim = childSims.get(childName);
//				//minChildName = childName;
//			}
//		}
//		return minChildSim;
//	}
	
//	/**
//	 * Returns name of child node with lowest similarity to this node.
//	 * 
//	 * @return Node name of child with lowest similarity (type String).
//	 */
//	public String minimumDescendantName() {
//		for(String childName : childNames) {
//			
//			if(childSims.get(childName) < minChildSim) {
//				minChildName = childName;
//			}
//		}
//		return minChildName;
//	}
	
//	/**
//	 * Adds name and similarity to HashMap if flagged as higher than a parent or child similarity. 
//	 * 
//	 * @param orgName
//	 * @param sim
//	 */
//	public void flagRelation(String orgName, double sim) {
//		flaggedRelationships.put(orgName, sim);
//	}
//	
//	/**
//	 * Returns HashMap holding all flagged relationships (nodes with higher similarity than
//	 * this nodes parent or the lowest similarity child).
//	 * @return HashMap
//	 */
//	public HashMap<String, Double> getFlaggedRelations(){
//		return flaggedRelationships;
//	}
	
	/**
	 * Add level value to node in the tree
	 * @param lvl Level of node in tree (type int).
	 */
	public void addLevel(int lvl) {
		level = lvl;
	}
	
	/**
	 * Adds hierarchical levels to nodes recursively through the tree.
	 * @param level_
	 */
	public void assignLevels(int level_, String parentRankName) {
		
		level = level_;
		if(taxonomicRank == null || taxonomicRank.equalsIgnoreCase("no rank") || taxonomicRank.equalsIgnoreCase("clade")) {
			taxonomicRank = parentRankName + ".1";
		}
		level_++;
		
		
		for(NCBITreeNode childNode : childNodes) {
			
			if(childNode.taxID != taxID) {
				childNode.assignLevels(level_, taxonomicRank);
				
			}
			//childNode.traverse(level_);
		}
		
	}
	
	/**
	 * Returns hierarchical level of the node.
	 * @return int level
	 */
	public int getLevel() {
		return level;
	}
	
	
	/**
	 * Tests to see if argument is descendant of this.
	 * A node is considered a descendant of itself.
	 * @param nodeB 
	 * @return true if nodeB is a descendant of this
	 */
	public boolean isDescendantOf(final NCBITreeNode nodeB) {
		if(this == nodeB) {return true;}
		else if(this.parentNode == this){return false;}
		else {return parentNode.isDescendantOf(nodeB);}
	}
	
	/**
	 * Tests if the input node is an ancestor of this node.
	 * Calls the isDescendantOf method on the input node
	 * @param nodeB NCBITreeNode the node that might be an ancestor of this node.
	 * @return boolean
	 */
	public boolean isAncestorOf(final NCBITreeNode nodeB) {
		return nodeB.isDescendantOf(this);
	}
	
	/**
	 * Resets the identity, nodesWithIdentity and votes values for all descendant nodes
	 * of this node.
	 */
	public void resetRecursively() {
		identity = 0;
		identitySum = 0;
		nodesWithIdentity = 0;
		sizeSum = 0;
		votes = 0;
		flaggedNode = false;
		color = "white";
		

		for(NCBITreeNode childNode : childNodes) {

			if(childNode.taxID != taxID) {
				childNode.resetRecursively();

			}
		}
	}
	
	/**
	 * Set nodes identity to the average identity of its descendants.
	 * 
	 * @param queryNode int ID of node relative to this node.
	 */
	public void percolateIdentityUp(int queryNode) {

		//If THIS nodes identity is greater than 0 (has an identity at all),
		//and if THIS nodes ID isn't equal to the input query node
		//nodes with identity becomes 1,
		//sizeSum (number of nodes accounted for) is set to size (usually 0),
		//identitySum (total of all identities) is set to THIS nodes identity.
		if(identity > 0 && nodeId != queryNode) {nodesWithIdentity = 1; sizeSum = size; identitySum = identity;}

		//Iterate over the child nodes of THIS node.
		for(NCBITreeNode childNode : childNodes) {

			//If this node isn't a child of itself (handles the life or "0" node).
			if(childNode.taxID != taxID) {

				//Calls percolateIdentityUp method recursively, pulling identities up to this node.
				childNode.percolateIdentityUp(queryNode);

				nodesWithIdentity+=childNode.nodesWithIdentity;
				identitySum+=childNode.identitySum;
				sizeSum+=childNode.sizeSum;
				votes += childNode.votes;

			}
		}
	}
	
	
	/**
	 * Returns the average identity of all descendant nodes to this node.
	 * @return int Average identity.
	 */
	public double averageIdentity() {
		if(nodesWithIdentity < 1) {return 0;}
		else {return identitySum / nodesWithIdentity;}
	}
	
	/**
	 * Returns the nodeID.
	 * @return int.
	 */
	public int getNodeId() {
		return nodeId;
	}
	
//	/**
//	 * Returns the similarity of this node to its parent or
//	 * the average identity of the parent node. Whichever value is higher.
//	 * @return double Similarity value.
//	 */
//	public double parentSimilarity() {
//		if(parentNode.identity > parentNode.averageIdentity()) {
//			return parentNode.identity;
//		} else {
//			return parentNode.averageIdentity();
//		}
//	}
	
	/**
	 * Returns the current averageIdentity value of this node to its parent.
	 * @return double
	 */
	public double parentSimilarity() {
		return parentNode.averageIdentity();
		//return parentNode.identity;
	}
	
	/**
	 * Method to call toDot method without requiring input.
	 * @return StringBuilder from toDot.
	 */
	public StringBuilder toDot(boolean printAllNodes_) {
		return toDot(null, printAllNodes_);
	}
	
	/**
	 * Method to implement tree structure in GraphViz .dot format.
	 * @param sb StringBuilder.
	 * @return StringBuilder with structure of tree in .dot format.
	 */
	public StringBuilder toDot(StringBuilder sb, boolean printAllNodes) {
		
		//If the input StringBuilder is null, start a new String Builder
		if(sb==null) {sb = new StringBuilder();}
		
		//Initialize first as true when sb is first initialized.
		boolean first = sb.length() == 0;
		
		//If first == true, begin adding node and edge information for the current node
		//to the StringBuilder.
		if(first) {
			
			//First line of a .dot file.
			sb.append("digraph g{\n");
		} 


		if(votes > 0 || printAllNodes) {//Node information for the .dot file.
			sb.append("\t" + nodeId + " [style=filled fillcolor=" + this.color + " label=\" Node ID= " + nodeId +"\\n"
					+ "Taxon ID= " + taxID + "\\n"
					+ "Tax Rank= " + taxonomicRank + "\\n"
					+ "ID= " + String.format("%.2f", identity) +"\\n"
					+ "Avg= " + String.format("%.2f", averageIdentity()) 
					+ "\\n" + "Votes = " + votes + "\"]\n");


			//Iterate over child nodes and recursively call this method
			//adding node connection to the StringBuilder.
			for(NCBITreeNode childNode : childNodes) {
				if(childNode != this && (printAllNodes || childNode.votes > 0)) {
					childNode.toDot(sb, printAllNodes);
				}
			}

			//iterate over child nodes and add edge information to the StringBuilder.
			for(NCBITreeNode childNode : childNodes) {
				if((printAllNodes || childNode.votes > 0)) {
					sb.append("\t" + nodeId + " -> " + childNode.nodeId + "\n");
				}
			}
		}
		//If first == true, append a final closing curly brace.
		if(first) {
			sb.append("}\n");
		}
		return sb;
	}
	
//	digraph g{
//		a [label="a\nidentity=98\navgident=5"]
//
//		a -> b
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * The children/descendant nodes of this node.
	 */
	List<NCBITreeNode> childNodes = new ArrayList<NCBITreeNode>();
	
	/**
	 * The parent node of this node.
	 */
	NCBITreeNode parentNode = null;
	
	/**
	 * HashMap holding node names and similarities flagged as higher than similarities with
	 * direct children.
	 */
	boolean flaggedNode = false;
	
	//Minimum similarity of all child nodes and this node.
	//double minChildSim = 100;
	
	//Name of child node with minimum similarity.
	//String minChildName = null;
	
	//HashMap holding the names and similarity values between any direct children nodes
	//and this node.
	//HashMap<String, Double> childSims = new HashMap<>();
	
	//Similarity percentage to the parent node.
	//double parSim = -1;
	
	/**
	 * Organisms name associated with this node.
	 */
	String orgName;
	
	/**
	 * Taxanomic ID of this organism.
	 */
	int taxID;
	/**
	 * int Node ID number.
	 */
	final int nodeId;
	
	String taxonomicRank;
	
	/**
	 * HashSet of direct children of this node.
	 */
	HashSet<Integer> childIDs=new HashSet<Integer>();
	
	/**
	 * The names of descendant nodes of this node.
	 */
	HashSet<Integer> descendentIDs=new HashSet<Integer>();
	
	/**
	 * Name of the parent of this node.
	 */
	int parentID;
	
	/**
	 * Node level within the tree. Important for taxonomic classification.
	 */
	int level;
	
	/**
	 * Sequence similarity between this node and another, query node.
	 */
	double identity = 0.0;
	
	long size = 0;
	
	long descendantSize = 0;
	
	int numDescendants = 0;
	
	/**
	 * The identities of all children nodes to this node.
	 */
	double identitySum = 0;
	
	/**
	 * The number of children nodes to this node that have identities.
	 */
	long nodesWithIdentity = 0;
	
	long sizeSum = 0;
	
	int votes = 0;
	
	boolean printAllNodes;
	
	String color = "white";
}
