package jasper;

public class NCBIComparison implements Comparable<NCBIComparison> {

	/**
	 * Object for storing sequence similarity values between two nodes
	 * @param queryID_ int Node ID of the primary query sequence.
	 * @param refID_ int Node ID of the reference sequence.
	 * @param identity_ double similarity value between both nodes.
	 */
	public NCBIComparison(int queryID_, int refID_, double identity_) {
		
		//Set the queryID as the input value
		this.queryID = queryID_;
		
		//Set the refID as the input value
		this.refNodeID = refID_;
		
		//Set the identity as the input value
		this.identity = identity_;
	}
	
	
	/**
	 * toString method to return the queryID, the refID and the similarity with some formatting.
	 */
	public String toString() {
		return "Query node ID = " + queryID + ", Reference node ID = " + refNodeID + 
				", Similarity identity = " + identity;
	}
	
	/**
	 * Allows sorting of Comparison objects in an ArrayList<Comparison>
	 */
	public int compareTo(NCBIComparison comparison){
		
		//If the identity of this comparison is equal to another comparison, return 0.
		if(identity == comparison.identity)  {
			return 0;  
		
		//but if the identity of this Comparison is greater than another, return 1.
		}else if(identity < comparison.identity)  {
			return 1;  
		
		//The only remaining option is that the identity of this Comparison is less than another Comparison.
		//Return -1
		}else  {
			return -1;  
		}  
	}
	
	/**
	 * Current node/organism of focus.
	 */
	int queryID;
	
	/**
	 * Node being compared to the queryID node.
	 */
	int refNodeID;
	
	/**
	 * Similarity between the queryID node and the refID node.
	 */
	double identity;
	
}
