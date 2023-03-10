
public class ResID implements Comparable<ResID>{
	String chainID;
	int seqnum;
	char icode;
	
	ResID(String chainID, int seqnum, char icode){
		this.chainID=chainID; this.seqnum=seqnum; this.icode=icode;
	}
	/*
	ResID(String MCA_res_index){
		String[] token;
		String chainID_seqnum;
		
		if(MCA_res_index.contains(".")){
			token=MCA_res_index.split("\\.");
			chainID_seqnum=token[0];
			this.icode=token[1].charAt(0);
		}
		else{
			chainID_seqnum=MCA_res_index;
			this.icode=' ';
		}
		
		if(chainID_seqnum.contains("\'")){
			this.chainID=chainID_seqnum.charAt(1);
			this.seqnum=Integer.parseInt(chainID_seqnum.substring(3));
		}
		else{
			this.chainID=chainID_seqnum.charAt(0);
			this.seqnum=Integer.parseInt(chainID_seqnum.substring(1));
		}
	}
	*/
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
