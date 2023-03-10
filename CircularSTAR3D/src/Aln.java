//import java.util.*;
//
//public class Aln implements Comparable<Aln>{
//    List<Integer> index1= new ArrayList<Integer>();
//    List<Integer> index2= new ArrayList<Integer>();
//    double rmsd;
//
//    public Aln(List<Integer> index1, List<Integer> index2, double rmsd) {this.index1=index1; this.index2=index2; this.rmsd=rmsd; }
//
//    public int compareTo(Aln o) {
//        return this.rmsd < o.rmsd ? 0 : 1;
//    }
//}

import java.util.*;
import org.ejml.data.*;

public class Aln implements Comparable<Aln>{

    List<Integer> rna1_index;
    List<Integer> rna2_index;

    List<Integer> rna1_loop_index;
    List<Integer> rna2_loop_index;

    List<Integer> rna1_stack_index;
    List<Integer> rna2_stack_index;

    Set<Integer> aligned_stack;
    double score;
    double loop_score;
    double stack_score;

    double rmsd;

    boolean contain_junction;
    List<List<Integer>> info;
    Set<Integer> seed_loop_sm;
//    Integer seed_len = 0;

    public Aln(List<Integer> rna1_index, List<Integer> rna2_index, List<Integer> rna1_loop_index,
               List<Integer> rna2_loop_index, List<Integer> rna1_stack_index, List<Integer> rna2_stack_index,
               double rmsd, double score, double loop_score, double stack_score, Set<Integer> aligned_stack,
               boolean contain_junction, List<List<Integer>> info, Set<Integer> seed_loop_sm
               ){
        this.rna1_index=rna1_index;
        this.rna2_index=rna2_index;
        this.rna1_loop_index = rna1_loop_index;
        this.rna2_loop_index = rna2_loop_index;
        this.rna1_stack_index = rna1_stack_index;
        this.rna2_stack_index = rna2_stack_index;
        this.rmsd=rmsd;
        this.score= score;
        this.loop_score = loop_score;
        this.stack_score = stack_score;
        this.aligned_stack = aligned_stack;
//        this.seed_len = seed_len;
        this.contain_junction = contain_junction;
        this.info = info;
        this.seed_loop_sm = seed_loop_sm;
    }
    @Override
    public int compareTo(Aln o) {
        if(this.score > o.score)
            return -1;
        else {
            if(this.score < o.score)
                return 1;
            else
                return this.rmsd < o.rmsd ? -1 : 1;
        }
    }

    public String toString() {
        String tmp = "";
        for (int i : this.rna2_index) {
            tmp += Integer.toString(i);
            tmp += ",";
        }
        tmp += "\n";
        return tmp;
    }
}

