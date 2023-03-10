
public class Residue {
	ResID rid;
	String symbol;
	
	Residue(ResID rid, String symbol){
		this.rid=rid;
		this.symbol=symbol;
	}
	public boolean equals(Object o){
		return (o instanceof Residue) &&
				(this.rid.equals(((Residue)o).rid) && this.symbol.equals(((Residue)o).symbol));
	}

	@Override
	public String toString(){
		return rid.toString()+"("+symbol+")";
	}
}
