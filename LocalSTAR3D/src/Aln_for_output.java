import java.util.*;

public class Aln_for_output implements Comparable<Aln_for_output>{

    List<Integer> rna1_index= new ArrayList<Integer>();
    List<Integer> rna2_index= new ArrayList<Integer>();

    List<Integer> rna1_loop_index= new ArrayList<Integer>();
    List<Integer> rna2_loop_index= new ArrayList<Integer>();

    List<Integer> rna1_stack_index= new ArrayList<Integer>();
    List<Integer> rna2_stack_index= new ArrayList<Integer>();

    Set<Integer> aligned_stack = new HashSet<>();
    double score = 0.;
    double loop_score = 0.;
    double stack_score = 0.;

    double rmsd = 0.;

    public Aln_for_output(List<Integer> rna1_index, List<Integer> rna2_index, List<Integer> rna1_loop_index, List<Integer> rna2_loop_index,
                          List<Integer> rna1_stack_index, List<Integer> rna2_stack_index,double rmsd,double score,double loop_score, double stack_score, Set<Integer> aligned_stack){
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
    }
    @Override
    public int compareTo(Aln_for_output o) {
        if(this.score > o.score)
            return -1;
        else {
            if(this.score < o.score)
                return 1;
            else
                return this.rmsd < o.rmsd ? -1 : 1;
        }
    }

}
