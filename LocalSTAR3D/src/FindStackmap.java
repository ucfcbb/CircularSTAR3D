import java.util.*;
import java.util.concurrent.*;

import org.ejml.data.*;

public class FindStackmap implements Callable<List<Stackmap>>{
	private int id;
	private int num_jobs;
	private List<Pair<Integer, Integer>> bp1, bp2;
	private List<Point> coord1, coord2;
	
	public FindStackmap(int id, int num_jobs, List<Pair<Integer, Integer>> bp1, List<Pair<Integer, Integer>> bp2, List<Point> coord1, List<Point> coord2){
		this.id=id; this.num_jobs=num_jobs; this.bp1=bp1; this.bp2=bp2; this.coord1=coord1; this.coord2=coord2;
	}
	
	private List<Stackmap> diagonal(int i, int j){
		int i1, j1, i2, j2;
		int i3, j3, i4, j4;
		double rmsd;
		List<Stackmap> min_len_stackmap=new ArrayList<Stackmap>();
		while(i<bp1.size()-STAR3D.min_stack_size+1 && j<bp2.size()-STAR3D.min_stack_size+1){
			i1=bp1.get(i).left;
			j1=bp1.get(i).right;
			i2=bp2.get(j).left;
			j2=bp2.get(j).right;

			i3=bp1.get(i+STAR3D.min_stack_size-1).left;
			j3=bp1.get(i+STAR3D.min_stack_size-1).right;
			i4=bp2.get(j+STAR3D.min_stack_size-1).left;
			j4=bp2.get(j+STAR3D.min_stack_size-1).right;

			if((i3-i1+1)==STAR3D.min_stack_size && (j1-j3+1)==STAR3D.min_stack_size &&
			(i4-i2+1)==STAR3D.min_stack_size && (j2-j4+1)==STAR3D.min_stack_size){
				List<Point> s1=new ArrayList<Point>(coord1.subList(i1, i1+STAR3D.min_stack_size));
				s1.addAll(coord1.subList(j1-STAR3D.min_stack_size+1, j1+1));
				List<Point> s2=new ArrayList<Point>(coord2.subList(i2, i2+STAR3D.min_stack_size));
				s2.addAll(coord2.subList(j2-STAR3D.min_stack_size+1, j2+1));
				
				DenseMatrix64F X=Geom.List2Matrix(s1);
				DenseMatrix64F Y=Geom.List2Matrix(s2);
				rmsd=Geom.superimpose(X, Y);

				//System.out.println(i1+" "+j1+" "+i2+" "+j2);
				if(rmsd<STAR3D.rmsd_cutoff) min_len_stackmap.add(new Stackmap(i1, i2, j1, j2, STAR3D.min_stack_size, rmsd));
			}
			i++;
			j++;	
		}

		// to-do: if longest e-stack rmsd >= cutoff, shorter e-stack will be discard by mistake.
		List<Stackmap> cont_stackmap=new ArrayList<Stackmap>();
		List<Stackmap> ret=new ArrayList<Stackmap>();
		for(Stackmap SM: min_len_stackmap){
			if(cont_stackmap.isEmpty()==true) cont_stackmap.add(SM);
			else{
				Stackmap preSM=cont_stackmap.get(cont_stackmap.size()-1);
				if((preSM.i1+1==SM.i1 && SM.j1+1==preSM.j1) && (preSM.i2+1==SM.i2 && SM.j2+1==preSM.j2)) cont_stackmap.add(SM);
				else{
					Stackmap firstSM=cont_stackmap.get(0);
					firstSM.size+=cont_stackmap.size()-1;
					
					List<Point> s1=new ArrayList<Point>(coord1.subList(firstSM.i1, firstSM.i1+firstSM.size));
					s1.addAll(coord1.subList(firstSM.j1-firstSM.size+1, firstSM.j1+1));
					List<Point> s2=new ArrayList<Point>(coord2.subList(firstSM.i2, firstSM.i2+firstSM.size));
					s2.addAll(coord2.subList(firstSM.j2-firstSM.size+1, firstSM.j2+1));
					
					DenseMatrix64F X=Geom.List2Matrix(s1);
					DenseMatrix64F Y=Geom.List2Matrix(s2);
					rmsd=Geom.superimpose(X, Y);
					firstSM.rmsd=rmsd;
					if(rmsd<STAR3D.rmsd_cutoff)
						ret.add(firstSM);
					cont_stackmap.clear();
					cont_stackmap.add(SM);
				}
			}
		}
		if(cont_stackmap.isEmpty()==false){
			Stackmap firstSM=cont_stackmap.get(0);
			firstSM.size+=cont_stackmap.size()-1;
			List<Point> s1=new ArrayList<Point>(coord1.subList(firstSM.i1, firstSM.i1+firstSM.size));
			s1.addAll(coord1.subList(firstSM.j1-firstSM.size+1, firstSM.j1+1));
			List<Point> s2=new ArrayList<Point>(coord2.subList(firstSM.i2, firstSM.i2+firstSM.size));
			s2.addAll(coord2.subList(firstSM.j2-firstSM.size+1, firstSM.j2+1));
			
			DenseMatrix64F X=Geom.List2Matrix(s1);
			DenseMatrix64F Y=Geom.List2Matrix(s2);
			rmsd=Geom.superimpose(X, Y);
			firstSM.rmsd=rmsd;

			if(rmsd<STAR3D.rmsd_cutoff)
				ret.add(firstSM);
			cont_stackmap.clear();
		}
		return ret;
	}
	
	public List<Stackmap> call(){
		List<Stackmap> ret=new ArrayList<Stackmap>();
		int count=0;
		for(int i=bp1.size(); i>=0; i--){
			if(count%num_jobs==id) ret.addAll(diagonal(i, 0));
			count++;
		}

		for(int i=1; i<bp2.size(); i++){
			if(count%num_jobs==id) ret.addAll(diagonal(0, i));
			count++;
		}
		return ret;
	}
}
