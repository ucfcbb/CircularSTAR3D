import java.util.*;
import java.util.concurrent.Callable;
import org.ejml.data.DenseMatrix64F;

public class FindCompatibleSM implements Callable<List<Pair<Integer, Integer>>>{
	private int id;
	private int num_jobs;
    private List<Integer> sm_index_list;
	private List<Stackmap> sm_list;
	private List<Point> coord1, coord2;
    public FindCompatibleSM(int id, int num_jobs, List<Integer> sm_index_list, List<Point> coord1, List<Point> coord2){
        this.id=id; this.num_jobs=num_jobs; this.sm_index_list = sm_index_list; this.coord1=coord1; this.coord2=coord2;
    }

    private boolean overlap(int s1, int l1, int s2, int l2){
        if((s2>=s1 && s2<=s1+l1-1) || (s1>=s2 && s1<=s2+l2-1))
            return true;
        else
            return false;
    }

    private boolean overlap_sm(Stackmap sm1, Stackmap sm2){
        if(overlap(sm1.i1, sm1.size, sm2.i1, sm2.size) ||
                overlap(sm1.i1, sm1.size, sm2.j1-sm2.size+1, sm2.size) ||
                overlap(sm1.j1-sm1.size+1, sm1.size, sm2.i1, sm2.size) ||
                overlap(sm1.j1-sm1.size+1, sm1.size, sm2.j1-sm2.size+1, sm2.size)){
            return true;
        }
        if(overlap(sm1.i2, sm1.size, sm2.i2, sm2.size) ||
                overlap(sm1.i2, sm1.size, sm2.j2-sm2.size+1, sm2.size) ||
                overlap(sm1.j2-sm1.size+1, sm1.size, sm2.i2, sm2.size) ||
                overlap(sm1.j2-sm1.size+1, sm1.size, sm2.j2-sm2.size+1, sm2.size)){
            return true;
        }
        return false;
    }

    private List<String> get_pattern(Pair<Integer, Integer> stk1, Pair<Integer, Integer> stk2){
        HashMap<Integer, String> pos_to_str = new HashMap<Integer, String>();
        Integer a1 = stk1.left;
        Integer a2 = stk1.right;
        Integer b1 = stk2.left;
        Integer b2 = stk2.right;

        pos_to_str.put(a1, "a1");
        pos_to_str.put(a2, "a2");
        pos_to_str.put(b1, "b1");
        pos_to_str.put(b2, "b2");

        List<Integer> sorted_stk = Arrays.asList(a1 ,a2, b1, b2);
        Collections.sort(sorted_stk);
        List<String> pattern = new ArrayList<String>();
        for(Integer i : sorted_stk)
            pattern.add(pos_to_str.get(i));

        return pattern;
    }

    private int pattern_match(List<String> p1, List<String> p2){
        if(p1.equals(p2)) {
            return 1; //match without rotate
        }

        for(int i=1;i<p1.size();i++) {
            List<String> rotated_list = new ArrayList<String>(p1.subList(i, p1.size()));
            rotated_list.addAll(p1.subList(0,i));

            if(rotated_list.equals(p2)) {
                return 2; //match after rotate
            }
        }

//        p1.addAll(p1);
//        if(p1.contains(p2)) {
//            return 2; //match after rotate
//        }

        return 0; //not match
    }

    private int SMC(Stackmap sm1, Stackmap sm2, List<Point>coord1, List<Point> coord2){
        if(overlap_sm(sm1, sm2) ){
            return 0;
        }

        Pair<Integer, Integer> s11 = new Pair<Integer, Integer>(sm1.i1, sm1.j1, sm1.size);
        Pair<Integer, Integer> s21 = new Pair<Integer, Integer>(sm1.i2, sm1.j2, sm1.size);

        Pair<Integer, Integer> s12 = new Pair<Integer, Integer>(sm2.i1, sm2.j1, sm2.size);
        Pair<Integer, Integer> s22 = new Pair<Integer, Integer>(sm2.i2, sm2.j2, sm2.size);

        List<String> p1 = get_pattern(s11, s12);
        List<String> p2 = get_pattern(s21, s22);

        int pattern_match_res = pattern_match(p1, p2);
        // pattern unmatch
        if(pattern_match_res == 0) {
            return 0;
        }

        List<Point> s1=new ArrayList<Point>(coord1.subList(sm1.i1, sm1.i1+sm1.size));
        s1.addAll(coord1.subList(sm1.j1-sm1.size+1, sm1.j1+1));
        s1.addAll(coord1.subList(sm2.i1, sm2.i1+sm2.size));
        s1.addAll(coord1.subList(sm2.j1-sm2.size+1, sm2.j1+1));

        List<Point> s2=new ArrayList<Point>(coord2.subList(sm1.i2, sm1.i2+sm1.size));
        s2.addAll(coord2.subList(sm1.j2-sm1.size+1, sm1.j2+1));
        s2.addAll(coord2.subList(sm2.i2, sm2.i2+sm2.size));
        s2.addAll(coord2.subList(sm2.j2-sm2.size+1, sm2.j2+1));

//        List<ResID> r1=new ArrayList<ResID>(STAR3D.ResID1_list.subList(sm1.i1, sm1.i1+sm1.size));
//        r1.addAll(STAR3D.ResID1_list.subList(sm1.j1-sm1.size+1, sm1.j1+1));
//        r1.addAll(STAR3D.ResID1_list.subList(sm2.i1, sm2.i1+sm2.size));
//        r1.addAll(STAR3D.ResID1_list.subList(sm2.j1-sm2.size+1, sm2.j1+1));
//
//        List<ResID> r2=new ArrayList<ResID>(STAR3D.ResID2_list.subList(sm1.i2, sm1.i2+sm1.size));
//        r2.addAll(STAR3D.ResID2_list.subList(sm1.j2-sm1.size+1, sm1.j2+1));
//        r2.addAll(STAR3D.ResID2_list.subList(sm2.i2, sm2.i2+sm2.size));
//        r2.addAll(STAR3D.ResID2_list.subList(sm2.j2-sm2.size+1, sm2.j2+1));

        DenseMatrix64F X=Geom.List2Matrix(s1);
        DenseMatrix64F Y=Geom.List2Matrix(s2);
        double rmsd=Geom.superimpose(X, Y);

//        System.out.println(rmsd);

        if(rmsd>STAR3D.rmsd_cutoff){
            return 0;
        }
//        return pattern_match_res;
        return (int)Math.ceil(rmsd);//8_20
    }

    public List<Pair<Integer, Integer>> call(){
        List<Pair<Integer, Integer>> ret=new ArrayList<Pair<Integer, Integer>>();

//        int count=0;
//        for(int i=0; i<sm_list.size(); i++){
//            if(count%num_jobs==id){
//                for(int j=i+1; j<sm_list.size(); j++){
//                    int v = SMC(sm_list.get(i), sm_list.get(j), coord1, coord2);
//                    if(v!=0)
//                        ret.add(new Pair<Integer, Integer>(i, j, v));
//                }
//            }
//            count++;
//        }

        int count=0;
        for(int i=0; i<sm_index_list.size(); i++){
            if(count%num_jobs==id){
                for(int j=i+1; j<sm_index_list.size(); j++){
                    Integer sm_index_i = sm_index_list.get(i);
                    Integer sm_index_j = sm_index_list.get(j);
                    if(sm_index_i == sm_index_j)
                        continue;

                    int v = SMC(STAR3D.SM_top.get(sm_index_i), STAR3D.SM_top.get(sm_index_j), coord1, coord2);
                    if(v!=0)
                        ret.add(new Pair<Integer, Integer>(sm_index_i, sm_index_j, v));
                }
            }
            count++;
        }


        return ret;
    }

}
