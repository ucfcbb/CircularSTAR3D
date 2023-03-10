
public class Point  implements Comparable<Point>{
	public double x;
	public double y;
	public double z;
	
	Point(double x, double y, double z){
		this.x=x; this.y=y; this.z=z;
	}
	
	@Override
	public String toString(){
		return "("+this.x+","+this.y+","+this.z+")";
	}
	
	@Override
	public boolean equals(Object o){
		if(!(o instanceof Point)) return false;
		Point other=(Point) o;
		return (this.x==other.x && this.y==other.y && this.z==other.z);
	}
	
	@Override
	public int compareTo(Point o){
		double rdist1=this.x*this.x+this.y*this.y+this.z*this.z;
		double rdist2=o.x*o.x+o.y*o.y+o.z*o.z;
		return Double.compare(rdist1, rdist2);
	}
}
