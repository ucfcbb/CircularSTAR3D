
public class Atom {
	int sn;
	Residue res;
	String atom;
	Point coord;
	float occupancy;
	float B_iso_or_equiv;

	public Atom(int sn, Residue res, String atom, Point coord, float occupancy, float B_iso_or_equiv){
		this.sn=sn;
		this.res=res;
		this.atom=atom;
		this.coord=coord;
		this.occupancy=occupancy;
		this.B_iso_or_equiv=B_iso_or_equiv;
	}

	public String toString(){
		return sn+"\t"+res.toString()+"\t"+atom+':'+coord;
	}

}
