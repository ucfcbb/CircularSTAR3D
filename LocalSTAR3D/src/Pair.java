public class Pair<T extends Comparable<? super T>, V extends Comparable<? super V>> implements Comparable<Pair<T, V>>{
	public final T left;
	public final T right;
	public final V v;
	public Pair(T left, T right, V v) {this.left=left; this.right=right; this.v=v; }
	@Override
	public String toString(){
		return "("+left.toString()+","+right.toString()+"):"+v;
	}
	public int compareTo(Pair<T, V> o) {
        int cmp = o == null ? 1 : (this.left).compareTo(o.left);
        return cmp == 0 ? (this.right).compareTo(o.right) : cmp;
    }
}

