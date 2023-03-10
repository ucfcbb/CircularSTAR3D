import java.io.Serializable;

public class Stackmap implements Comparable<Stackmap>, Serializable {
	int i1, i2, j1, j2, size;
	double rmsd;
	Stackmap(int i1, int i2, int j1, int j2, int size, double rmsd){
		this.i1=i1; this.i2=i2; this.j1=j1; this.j2=j2; this.size=size; this.rmsd=rmsd;
	}
	
	public int compareTo(Stackmap o){
		return this.rmsd-0.1*this.size > o.rmsd-0.1*o.size ? 1: (this.rmsd-0.1*this.size < o.rmsd-0.1*o.size ? -1:0);
	}
	
	public String toString(){
		return "("+STAR3D.ResID1_list.get(this.i1)+","+STAR3D.ResID1_list.get(this.j1)+")"+"<->("+
				STAR3D.ResID2_list.get(this.i2)+","+STAR3D.ResID2_list.get(this.j2)+"):"+this.size;
	}
}
