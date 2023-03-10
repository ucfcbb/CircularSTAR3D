import java.util.ArrayList;
import java.util.List;

public class LoopAlnRes{
    public List<Pair<Integer, Integer>> ret_nt;
    public double ret_score;

    public LoopAlnRes(){
        this.ret_nt=new ArrayList<Pair<Integer, Integer>>();
        this.ret_score=0.;
    }

    public LoopAlnRes(List<Pair<Integer, Integer>> ret_nt, double ret_score) {
        this.ret_nt = ret_nt;
        this.ret_score = ret_score;
    }
}