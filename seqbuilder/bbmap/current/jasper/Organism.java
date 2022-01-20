package jasper;


public class Organism {

	String orgName;
	int taxId;
	
	// This is the constructor of the class Organism
	public Organism(int id, String name) {
        this.taxId = id;
        this.orgName = name;
	    }
	
	//public void addName(String addName) {
	//    orgName = addName;
	//    }
	
		
	public void printOrg() {
	    System.out.println("ID:"+ taxId );
	    System.out.println("Name:" + orgName );
	   }

}
