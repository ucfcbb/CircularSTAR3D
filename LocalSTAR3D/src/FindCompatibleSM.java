import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.ejml.data.DenseMatrix64F;


public class FindCompatibleSM implements Callable<List<Pair<Integer, Integer>>>{
	private int id;
	private int num_jobs;
	private List<Stackmap> sm_list;
	private List<Point> coord1, coord2;
	
	public FindCompatibleSM(int id, int num_jobs, List<Stackmap>sm_list, List<Point> coord1, List<Point> coord2){
		this.id=id; this.num_jobs=num_jobs; this.sm_list=sm_list; this.coord1=coord1; this.coord2=coord2;
	}
	
	private boolean overlap(int s1, int l1, int s2, int l2){
		if((s2>=s1 && s2<=s1+l1-1) || (s1>=s2 && s1<=s2+l2-1)) 
			return true;
		else 
			return false;
	}
	
	private int relation(int s1, int e1, int s2, int e2){
		if(e1<s2) return 1;					// 1 <J 2
		else if (e2<s1) return 2; 			//2 <J 1
		else if (s2>s1 && e2<e1) return 3; 	// 1 <I 2
		else return 4;						// 2 <I 1
	}
	
	private int SMC(Stackmap sm1, Stackmap sm2, List<Point>coord1, List<Point> coord2){
		if(overlap(sm1.i1, sm1.size, sm2.i1, sm2.size) || overlap(sm1.i1, sm1.size, sm2.j1-sm2.size+1, sm2.size) ||
		   overlap(sm1.j1-sm1.size+1, sm1.size, sm2.i1, sm2.size) || overlap(sm1.j1-sm1.size+1, sm1.size, sm2.j1-sm2.size+1, sm2.size)){
			return 0;
		}
		if(overlap(sm1.i2, sm1.size, sm2.i2, sm2.size) || overlap(sm1.i2, sm1.size, sm2.j2-sm2.size+1, sm2.size) ||
		   overlap(sm1.j2-sm1.size+1, sm1.size, sm2.i2, sm2.size) || overlap(sm1.j2-sm1.size+1, sm1.size, sm2.j2-sm2.size+1, sm2.size)){
			return 0;
		}
//		if(relation(sm1.i1, sm1.j1, sm2.i1, sm2.j1)!=relation(sm1.i2, sm1.j2, sm2.i2, sm2.j2)){
//			return 0;
//		}

		int relation1 = relation(sm1.i1, sm1.j1, sm2.i1, sm2.j1);
		int relation2 = relation(sm1.i2, sm1.j2, sm2.i2, sm2.j2);
		if(relation1 !=relation2) return 0;

		List<Point> s1=new ArrayList<Point>(coord1.subList(sm1.i1, sm1.i1+sm1.size));
		s1.addAll(coord1.subList(sm1.j1-sm1.size+1, sm1.j1+1));
		s1.addAll(coord1.subList(sm2.i1, sm2.i1+sm2.size));
		s1.addAll(coord1.subList(sm2.j1-sm2.size+1, sm2.j1+1));

		List<Point> s2=new ArrayList<Point>(coord2.subList(sm1.i2, sm1.i2+sm1.size));
		s2.addAll(coord2.subList(sm1.j2-sm1.size+1, sm1.j2+1));
		s2.addAll(coord2.subList(sm2.i2, sm2.i2+sm2.size));
		s2.addAll(coord2.subList(sm2.j2-sm2.size+1, sm2.j2+1));
		
		DenseMatrix64F X=Geom.List2Matrix(s1);
		DenseMatrix64F Y=Geom.List2Matrix(s2);
		double rmsd=Geom.superimpose(X, Y);
		if(rmsd>STAR3D.rmsd_cutoff) return 0;
//		else return 1;
		else return relation1;

//        return 0;
	}
	
	public List<Pair<Integer, Integer>> call(){
		List<Pair<Integer, Integer>> ret=new ArrayList<Pair<Integer, Integer>>();
		
		int count=0;
		for(int i=0; i<sm_list.size(); i++){
			if(count%num_jobs==id){
				for(int j=i; j<sm_list.size(); j++){
					int v = SMC(sm_list.get(i), sm_list.get(j), coord1, coord2);
                    			if(v!=0)
                        			ret.add(new Pair<Integer, Integer>(i, j, v));
				}
			}
			count++;
		}
			
		return ret;
	}

}
