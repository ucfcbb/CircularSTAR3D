
public class ResID implements Comparable<ResID>{
	String chainID;
	int seqnum;
	char icode;
	
	ResID(String chainID, int seqnum, char icode){
		this.chainID=chainID; this.seqnum=seqnum; this.icode=icode;
	}

	public String toString(){
		return this.chainID+":"+Integer.toString(this.seqnum)+String.valueOf(this.icode).trim();
	}
	
	public boolean equals(Object o){
		return (o instanceof ResID) &&
                (this.chainID.equals(((ResID)o).chainID) && this.seqnum==((ResID)o).seqnum && this.icode==((ResID)o).icode);
	}
	
	public int hashCode(){
        int tempChainID = 0;
        for(Character c : chainID.toCharArray()) {
            tempChainID += c;
        }
        return tempChainID*1000000+this.icode*10000+this.seqnum;
	}
	
	public int compareTo(ResID o){
		int cmp = o == null ? 1 : this.seqnum - o.seqnum;
        return cmp == 0 ? this.icode - o.icode : cmp;
	}
}
