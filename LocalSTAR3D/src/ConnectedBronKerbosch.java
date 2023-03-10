import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.Callable;

public class ConnectedBronKerbosch implements Callable<List<Set<Integer>>>{

    private int[][] SMC_graph;
    private int[][] SM_connected_graph;
    private Set<Integer> SMC_graph_vertex;
    private List<Set<Integer>> SM_clique;


    public ConnectedBronKerbosch(int[][] SMC_graph, int[][] SM_connected_graph, Set<Integer> SMC_graph_vertex,List<Set<Integer>> SM_clique){
        this.SMC_graph = SMC_graph; this.SM_connected_graph = SM_connected_graph; this.SMC_graph_vertex = SMC_graph_vertex;this.SM_clique = SM_clique;
    }

    //find the neighbor of i in graph G
    private static Set<Integer> GN(int[][] G, int i){
        Set<Integer> ret=new HashSet<Integer>();
        for(int j=0; j<G[i].length; j++){
            if(G[i][j]!=0 && j!=i) ret.add(j);
        }
        return ret;
    }
    public static void connected_Bron_Kerbosch(int[][] G, int[][] GC, Set<Integer> R, Set<Integer> P, Set<Integer> X, Set<Integer> Q, List<Set<Integer>> ret) {
        //GC:connecte graph
        //Q:vertex connected to R in 1 step
        Set<Integer> P_cur = new HashSet<Integer>(P);

        Set<Integer> PX = new HashSet<Integer>(P);
        PX.addAll(X);

        int i=new Random().nextInt(PX.size());
        int j=0, u=0;
        for(Integer I: PX){
            if(j==i) u=I;
            j++;
        }

        Set<Integer> PGNU=new HashSet<Integer>(P_cur);

        PGNU.removeAll(GN(G, u));
        for(Integer v : PGNU){
        //for (Integer v : P_cur) {

            Set<Integer> RV = new HashSet<Integer>(R);
            RV.add((Integer)v);

            Set<Integer> QINV = new HashSet<Integer>(Q);
            for (Integer I : GN(GC, v)){
                if(!R.contains(I))
                    QINV.add(I);
            }

            Set<Integer> PINV = new HashSet<Integer>();
            for (Integer I : GN(G, v))
                if (P.contains(I)) PINV.add(I);

            Set<Integer> XINV = new HashSet<Integer>();
            for (Integer I : GN(G, v))
                if (X.contains(I) /*&& GN(GC,v).contains(I)*/) XINV.add(I);

            connected_Bron_Kerbosch_helper(G, GC, RV, PINV, XINV, QINV, ret);
            P.remove((Integer) v);
            X.add((Integer) v);
        }
    }

    public static void connected_Bron_Kerbosch_helper(int[][] G, int[][] GC, Set<Integer> R, Set<Integer> P, Set<Integer> X, Set<Integer> Q, List<Set<Integer>> ret) {
        //GC:connecte graph
        //Q:vertex connected to R in 1 step

        Set<Integer> PIQ = new HashSet<Integer>();
        for (Integer I : Q){
            if(P.contains(I))
                PIQ.add(I);
        }
        Set<Integer> XIQ = new HashSet<Integer>();
        for (Integer I : Q){
            if(X.contains(I))
                XIQ.add(I);
        }

        if (PIQ.isEmpty() && XIQ.isEmpty()) {
            ret.add(R);
            return;
        }

        for (Integer v : PIQ) {

            Set<Integer> RV = new HashSet<Integer>(R);
            RV.add((Integer)v);

            Set<Integer> QINV = new HashSet<Integer>(Q);
            for (Integer I : GN(GC, v)){
                if(!R.contains(I))
                    QINV.add(I);
            }

            Set<Integer> PINV = new HashSet<Integer>();
            for (Integer I : GN(G, v))
                if (P.contains(I)) PINV.add(I);

            Set<Integer> XINV = new HashSet<Integer>();
            for (Integer I : GN(G, v))
                if (X.contains(I) /*&& GN(GC,v).contains(I)*/) XINV.add(I);

            connected_Bron_Kerbosch_helper(G, GC, RV, PINV, XINV, QINV, ret);
            P.remove((Integer) v);
            X.add((Integer) v);
        }
    }
    public List<Set<Integer>> call(){
        connected_Bron_Kerbosch(SMC_graph, SM_connected_graph, new HashSet<Integer>(), SMC_graph_vertex, new HashSet<Integer>(), new HashSet<Integer>(),SM_clique);
        return SM_clique;
    }
}
