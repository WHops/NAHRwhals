package jasper;

public class Comparison {

	/**
	 * Object for storing sequence similarity values between two nodes
	 * @param queryID_ int Node ID of the primary query sequence.
	 * @param refID_ int Node ID of the reference sequence.
	 * @param identity_ double similarity value between both nodes.
	 */
	public Comparison(int queryID_, int refID_, double identity_) {
		
		//Set the queryID as the input value
		this.queryID = queryID_;
		
		//Set the refID as the input value
		this.refID = refID_;
		
		//Set the identity as the input value
		this.identity = identity_;
	}
	/**
	 * toString method to return the queryID, the refID and the similarity with some formatting.
	 */
	public String toString() {
		return "Query node ID = " + queryID + ", Reference node ID = " + refID + 
				", Similarity identity = " + identity;
	}
	
	/**
	 * Current node/organism of focus.
	 */
	int queryID;
	
	/**
	 * Node being compared to the queryID node.
	 */
	int refID;
	
	/**
	 * Similarity between the queryID node and the refID node.
	 */
	double identity;
	
}
