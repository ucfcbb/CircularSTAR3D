import java.io.Serializable;
import java.lang.*;

public class Pair<T extends Comparable<? super T>, V extends Comparable<? super V>> implements Comparable<Pair<T, V>>, Serializable {
	public final T left;
	public final T right;
	public final V v;
	public Pair(T left, T right, V v) {this.left=left; this.right=right; this.v=v; }
	@Override
	public String toString(){
		return "("+left.toString()+","+right.toString()+"):"+v;
	}

	public String toString(int a) {
		if ((Integer) left==-1 || (Integer) right==-1)
			return "("+left.toString()+","+right.toString()+"):"+v;
		if(a==1)
			return "(" + STAR3D.ResID1_list.get((Integer) left).toString() + "," + STAR3D.ResID1_list.get((Integer) right).toString() + "):" + v;
		else
			return "(" + STAR3D.ResID2_list.get((Integer) left).toString() + "," + STAR3D.ResID2_list.get((Integer) right).toString() + "):" + v;
	}

	public int compareTo(Pair<T, V> o) {
        int cmp = o == null ? 1 : (this.left).compareTo(o.left);
        return cmp == 0 ? (this.right).compareTo(o.right) : cmp;
    }

	public boolean equals(Object o){
//		return (o instanceof Pair) &&
//				(this.left.equals(((Pair)o).left) && this.right.equals(((Pair)o).right) && this.v.equals(((Pair)o).v));
		return (o instanceof Pair) &&
				(this.left.equals(((Pair)o).left) && this.right.equals(((Pair)o).right));
	}

//	public int hashCode(){
//		return this.left.hashCode()*1000000 + this.right.hashCode()*10000+ this.v.hashCode();
//	}
	public int hashCode(){
		return this.left.hashCode()*1000 + this.right.hashCode();
	}

}

