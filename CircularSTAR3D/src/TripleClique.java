import java.util.*;
import java.util.concurrent.Callable;

public class TripleClique implements Callable<List<Map<Set<Integer>, Integer>>> {

    private Map<Integer, Set<Integer>> SM_connected_graph;
    private Map<Pair<Integer, Integer>, Map<Integer, Integer>> SMC_2_to_1;
    private List<Set<Integer>> SM_clique;
    private List<Integer> rotation;
    private List<Map<Set<Integer>, Integer>> clique_to_rotation;

    public TripleClique(Map<Integer, Set<Integer>> sm_conn_graph, Map<Pair<Integer, Integer>, Map<Integer, Integer>> SMC_2_to_1, List<Set<Integer>> SM_clique, List<Integer> rotation){
        this.SM_connected_graph = sm_conn_graph;
        this.SMC_2_to_1 = SMC_2_to_1;
        this.SM_clique = SM_clique;
        this.rotation = rotation;
        this.clique_to_rotation = new ArrayList<Map<Set<Integer>, Integer>>();
    }

    //for each new v, construct key (v, i) for i in R, search all value of the key in SMC_2_to_1
    public static Set<Integer> GN(Map<Pair<Integer, Integer>, Map<Integer, Integer>> SMC_2_to_1, Map<Integer, Integer> R, Integer v){
        Set<Integer> ret=new HashSet<Integer>();
        Pair<Integer, Integer> p;
        for(Integer i : R.keySet()){
            if(i<v)
                p = new Pair<Integer, Integer>(i,v,0);
            else
                p = new Pair<Integer, Integer>(v,i,0);

            if(!SMC_2_to_1.containsKey(p)) continue;

            if(ret.isEmpty())
                ret = SMC_2_to_1.get(p).keySet();
            else
                ret.retainAll(SMC_2_to_1.get(p).keySet());
        }
        return ret;
    }

    public static void get_triplet_clique_search(
            Map<Integer, Set<Integer>> sm_conn_graph,
            Map<Pair<Integer, Integer>, Map<Integer, Integer>> SMC_2_to_1,
            Map<Integer, Integer> R,
            Map<Integer, Integer> P,
            Set<Integer> Q,
            Set<Integer> X,
            List<Set<Integer>> SM_clique,
            List<Integer> rotation){
        Map<Integer, Integer> PIQ = new HashMap<Integer, Integer>(P);
        PIQ.keySet().retainAll(Q);

        Set<Integer> XIQ = new HashSet<Integer>(X);
        XIQ.retainAll(Q);

        if (PIQ.isEmpty() && XIQ.isEmpty()) {
            SM_clique.add(R.keySet());
            //add rotation info
            for(Integer i : R.values()){
                if(i>0) {
                    rotation.add(1);
                    return;
                }
            }
            rotation.add(0);
            return;
        }

        for (Integer v : PIQ.keySet()) {
            //if(R.contains(v)) continue;
            Map<Integer, Integer> RV = new HashMap<Integer, Integer>(R);
            RV.put(v, PIQ.get(v));

            Set<Integer> GNv = GN(SMC_2_to_1, R,v);

            Set<Integer> QINV = new HashSet<Integer>(Q);
            QINV.addAll(sm_conn_graph.get(v));

            Map<Integer, Integer> PINV = new HashMap<Integer, Integer>(P);
            PINV.keySet().retainAll(GNv);

            Set<Integer> XINV = new HashSet<Integer>(X);
            XINV.retainAll(GNv);

            get_triplet_clique_search(sm_conn_graph, SMC_2_to_1, RV, PINV, Q, XINV, SM_clique, rotation);
            P.remove((Integer) v);
            X.add((Integer) v);
        }
    }

    public static void get_triplet_clique(
            Map<Integer, Set<Integer>> sm_conn_graph,
            Map<Pair<Integer, Integer>, Map<Integer, Integer>> SMC_2_to_1,
            List<Set<Integer>> SM_clique,
            List<Integer> rotation){
        for(Pair<Integer, Integer> p : SMC_2_to_1.keySet()){
            // initiate R as two connected nodes
            if(!(sm_conn_graph.get(p.left).contains(p.right)))
                continue;
            Map<Integer, Integer> R = new HashMap<Integer, Integer>();
            R.put(p.left, 0);
            R.put(p.right,0);

            // initial P as a set of the 3rd stack in the smc
            Map<Integer, Integer> P = SMC_2_to_1.get(p);

            //initial Q
            Set<Integer> Q = new HashSet<Integer>();
            for(Integer i : R.keySet())
                Q.addAll(sm_conn_graph.get(i));

            get_triplet_clique_search(sm_conn_graph, SMC_2_to_1, R, P, Q, new HashSet<Integer>(), SM_clique, rotation);
        }
    }

    public List<Map<Set<Integer>, Integer>> call(){
        get_triplet_clique(SM_connected_graph, SMC_2_to_1, SM_clique, rotation);
//        System.out.println("SM_clique.size() "+ SM_clique.size());
//        System.out.println("rotation.size() "+ rotation.size());
        for(int i = 0; i<SM_clique.size();i++){
            Map<Set<Integer>, Integer> tmp = new HashMap<Set<Integer>, Integer>();
            tmp.put(SM_clique.get(i), rotation.get(i));
//            System.out.println("get tmp");
            clique_to_rotation.add(tmp);
        }
        //return SM_clique;
        return clique_to_rotation;
    }
}


