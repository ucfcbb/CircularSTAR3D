import java.util.*;
import java.util.concurrent.Callable;

import org.ejml.data.*;
import org.ejml.ops.CommonOps;

public class LoopAlign implements Callable<LoopAlnRes>{
	private int id, num_jobs;
	List<Pair<Integer, Integer>> loop1, loop2;
	private Map<Integer, ArrayList<Integer>> bp1_map, bp2_map;
	private List<Point> coord1, coord2;
	private DenseMatrix64F R;
	private Point XC, YC;


	LoopAlign(int id, int num_jobs, List<Pair<Integer, Integer>> loop1, List<Pair<Integer, Integer>> loop2, Map<Integer, ArrayList<Integer>> bp1_map, Map<Integer, ArrayList<Integer>> bp2_map, List<Point> coord1, List<Point> coord2, DenseMatrix64F R, Point XC, Point YC){
		this.id=id;
		this.num_jobs=num_jobs;
		this.loop1=loop1;
		this.loop2=loop2;
		this.bp1_map=bp1_map;
		this.bp2_map=bp2_map;
		this.coord1=coord1;
		this.coord2=coord2;
		this.R=R;
		this.XC=XC;
		this.YC=YC;
	}

	private Pair<Integer, Integer> get_base_neighbor(int len1, int i1, int len2, int i2){
		int L=Math.min(Math.min(i1, 1), Math.min(i2, 1));
		int R=Math.min(Math.min(len1-1-i1, 1), Math.min(len2-1-i2, 1));
		return new Pair<Integer, Integer>(L, R, 0);
	}

	private double base_superimpose(int i, int j){
		Pair<Integer, Integer> Reg=get_base_neighbor(coord1.size(), i, coord2.size(), j);
		int RegL=Reg.left;
		int RegR=Reg.right;

		List<Point> s1=new ArrayList<Point>(coord1.subList(i-RegL, i+RegR+1));
		List<Point> s2=new ArrayList<Point>(coord2.subList(j-RegL, j+RegR+1));

		DenseMatrix64F X=Geom.List2Matrix(s1);
		DenseMatrix64F Y=Geom.List2Matrix(s2);

		DenseMatrix64F XM=Geom.translation(X, XC);
		DenseMatrix64F YM=Geom.translation(Y, YC);

		DenseMatrix64F XT=new DenseMatrix64F(XM.numCols, XM.numRows);
		CommonOps.transpose(XM, XT);
		DenseMatrix64F RX=new DenseMatrix64F(3, XT.numCols);
		CommonOps.mult(R, XT, RX);

		CommonOps.transpose(RX);
		return Geom.RMSD(RX, YM);
	}

	private double pair_superimpose(int i1, int j1, int i2, int j2){
		Pair<Integer, Integer> Reg1=get_base_neighbor(coord1.size(), i1, coord2.size(), i2);
		Pair<Integer, Integer> Reg2=get_base_neighbor(coord1.size(), j1, coord2.size(), j2);

		int Reg1L=Reg1.left;
		int Reg1R=Reg1.right;
		int Reg2L=Reg2.left;
		int Reg2R=Reg2.right;

		List<Point> s1=new ArrayList<Point>(coord1.subList(i1-Reg1L, i1+Reg1R+1));
		s1.addAll(coord1.subList(j1-Reg2L, j1+Reg2R+1));
		List<Point> s2=new ArrayList<Point>(coord2.subList(i2-Reg1L, i2+Reg1R+1));
		s2.addAll(coord2.subList(j2-Reg2L, j2+Reg2R+1));

		DenseMatrix64F X=Geom.List2Matrix(s1);
		DenseMatrix64F Y=Geom.List2Matrix(s2);

		DenseMatrix64F XM=Geom.translation(X, XC);
		DenseMatrix64F YM=Geom.translation(Y, YC);

		DenseMatrix64F XT=new DenseMatrix64F(XM.numCols, XM.numRows);
		CommonOps.transpose(XM, XT);
		DenseMatrix64F RX=new DenseMatrix64F(3, XT.numCols);
		CommonOps.mult(R, XT, RX);

		CommonOps.transpose(RX);
		return Geom.RMSD(RX, YM);
	}

	private int multi_pairing_align(int ix, int iy){
		List<Integer> bp1_list=bp1_map.get(ix);
		List<Integer> bp2_list=bp2_map.get(iy);
		int bp1_num=bp1_list.size();
		int bp2_num=bp2_list.size();

		int[][] M=new int[1+bp1_num][1+bp2_num];
		M[0][0]=0;
		for(int i=1; i<=bp1_num; i++) M[i][0]=0;
		for(int j=1; j<=bp2_num; j++) M[0][j]=0;

		int M1, M2, M3;
		for(int i=1; i<=bp1_num; i++)
			for(int j=1; j<=bp2_num; j++){
				M1=M[i-1][j-1]+ (pair_superimpose(ix, bp1_list.get(i-1), iy, bp2_list.get(j-1))<STAR3D.rmsd_cutoff ? 1:0);
				M2=M[i-1][j];
				M3=M[i][j-1];
				M[i][j]=Math.max(M1, Math.max(M2, M3));
			}
		return M[bp1_num][bp2_num];
	}
	//...i....j...
//...x....x...
//...[....]...
//...y....y...
//...i....j...
	private LoopAlnRes affine_gap_align(int ix, int jx, int iy, int jy){
//		if(jx != -1 && iy != -1 && jx!= -1 && iy != -1) {
//			System.out.println("breakpoint");
//		}

		List<Pair<Integer, Integer>> ret=new ArrayList<Pair<Integer, Integer>>();
		double optimal_score=0.;
		LoopAlnRes optimal_ret = new LoopAlnRes(ret,optimal_score);
		if(ix==-1 && jx==-1){
			if (jy != iy)
				optimal_ret.ret_score = STAR3D.gap_open_penalty+STAR3D.gap_extend_penalty*(jy - iy);
			return optimal_ret;
		}
		if(iy==-1 && jy==-1){
			if (jx != ix)
				optimal_ret.ret_score = STAR3D.gap_open_penalty+STAR3D.gap_extend_penalty*(jy - iy);
			return optimal_ret;
		}

		int seq1_len=jx-ix+1;
		int seq2_len=jy-iy+1;

		double[][] M=new double[1+seq1_len][1+seq2_len];
		double[][] IX=new double[1+seq1_len][1+seq2_len];
		double[][] IY=new double[1+seq1_len][1+seq2_len];

		int[][] TBM=new int[1+seq1_len][1+seq2_len];
		int[][] TBIX=new int[1+seq1_len][1+seq2_len];
		int[][] TBIY=new int[1+seq1_len][1+seq2_len];

		M[0][0]=IX[0][0]=IY[0][0]=0;

		for(int i=1; i<=seq1_len; i++){
			M[i][0]=-65535.;
			IX[i][0]=STAR3D.gap_open_penalty+STAR3D.gap_extend_penalty*(i-1);
			IY[i][0]=-65535.;
		}

		for(int j=1; j<=seq2_len; j++){
			M[0][j]=-65535.;
			IX[0][j]=-65535.;
			IY[0][j]=STAR3D.gap_open_penalty+STAR3D.gap_extend_penalty*(j-1);
		}

		double base_rmsd=0.;
		double base_score=0.;
		int num_of_pair=0;
		double M1, M2, M3, IX1, IX2, IX3, IY1, IY2, IY3;

		int tb_i=0;
		for(int i=1; i<=seq1_len; i++){
			for(int j=1; j<=seq2_len; j++){
				base_rmsd=base_superimpose(ix+i-1, iy+j-1);
				if(base_rmsd>=2*STAR3D.rmsd_cutoff) base_score=-65535;
				else if(base_rmsd<0.5*STAR3D.rmsd_cutoff) base_score=STAR3D.match_score;
				else if(base_rmsd<STAR3D.rmsd_cutoff) base_score=0.5*STAR3D.match_score;
				else base_score=STAR3D.mismatch_score;

				if(bp1_map.get(ix+i-1)!=null && bp2_map.get(iy+j-1)!=null)
					num_of_pair=multi_pairing_align(ix+i-1, iy+j-1);

				M1=M[i-1][j-1]+(base_score+num_of_pair);
				M2=IX[i-1][j-1]+(base_score+num_of_pair);
				M3=IY[i-1][j-1]+(base_score+num_of_pair);

				M[i][j]=Math.max(M1, Math.max(M2, M3));

				if(M[i][j]==M1) TBM[i][j]=1;
				else if(M[i][j]==M2) TBM[i][j]=2;
				else TBM[i][j]=3;

				IX1=M[i-1][j]+STAR3D.gap_open_penalty;
				IX2=IX[i-1][j]+STAR3D.gap_extend_penalty;
				IX3=IY[i-1][j]+STAR3D.gap_open_penalty;

				IX[i][j]=Math.max(IX1, Math.max(IX2, IX3));
				if(IX[i][j]==IX1) TBIX[i][j]=1;
				else if(IX[i][j]==IX2) TBIX[i][j]=2;
				else TBIX[i][j]=3;

				IY1=M[i][j-1]+STAR3D.gap_open_penalty;
				IY2=IX[i][j-1]+STAR3D.gap_open_penalty;
				IY3=IY[i][j-1]+STAR3D.gap_extend_penalty;

				IY[i][j]=Math.max(IY1, Math.max(IY2, IY3));
				if(IY[i][j]==IY1) TBIY[i][j]=1;
				else if(IY[i][j]==IY2) TBIY[i][j]=2;
				else TBIY[i][j]=3;
			}
		}
		optimal_score=Math.max(M[seq1_len][seq2_len], Math.max(IX[seq1_len][seq2_len], IY[seq1_len][seq2_len]));

//		System.out.print("loopalign ");
//		System.out.println(ix +" "+ jx +" "+ iy +" "+jy);

		if(optimal_score==M[seq1_len][seq2_len]) tb_i=1;
		else if(optimal_score==IX[seq1_len][seq2_len]) tb_i=2;
		else tb_i=3;

		List<Pair<Integer, Integer>> aln_map=affine_gap_TB(TBM, TBIX, TBIY, tb_i);
		for(Pair<Integer, Integer> P: aln_map){
			if(P.left!=-1 && P.right!=-1) ret.add(new Pair<Integer, Integer>(P.left+ix, P.right+iy, 0));
		}

		optimal_ret.ret_nt = ret;
		optimal_ret.ret_score = optimal_score;
//        System.out.print(ix +"\t"+jx+"\t");
//		System.out.println(optimal_score);
		return optimal_ret;
	}

	private List<Pair<Integer, Integer>> affine_gap_TB(int[][] TBM, int[][] TBIX, int[][] TBIY, int tb_i){
		int i=TBM.length-1;
		int j=TBM[0].length-1;

		List<Pair<Integer, Integer>> aln_map=new ArrayList<Pair<Integer, Integer>>();
		while(i!=0 && j!=0){
			if(tb_i==1) {
				aln_map.add(0, new Pair<Integer, Integer>(i-1, j-1, 0));
				tb_i=TBM[i][j];
				i--;
				j--;
			}
			else if(tb_i==2) {
				aln_map.add(0, new Pair<Integer, Integer>(i-1, -1, 0));
				tb_i=TBIX[i][j];
				i--;
			}
			else{
				aln_map.add(0, new Pair<Integer, Integer>(-1, j-1, 0));
				tb_i=TBIY[i][j];
				j--;
			}
		}
		if(i==0){
			for(int k=0; k<j; k++) aln_map.add(0, new Pair<Integer, Integer>(-1, j-k-1, 0));
		}
		if(j==0){
			for(int k=0; k<i; k++) aln_map.add(0, new Pair<Integer, Integer>(i-k-1, -1, 0));
		}
		return aln_map;
	}

	public LoopAlnRes call(){
		LoopAlnRes optimal_rets = new LoopAlnRes();
		int count=0;
		for(int i=0; i<loop1.size(); i++){
			if(count%num_jobs==id){
				LoopAlnRes optimal_ret = new LoopAlnRes();
				if(STAR3D.dangle_end == false && /*loop1.get(i).right-loop1.get(i).left>STAR3D.connect_distance  || loop2.get(i).right-loop2.get(i).left>STAR3D.connect_distance
							||*/ (loop1.get(i).v == 1 || loop2.get(i).v == 1)){
						optimal_ret = affine_gap_align(-1 ,-1,-1,-1);
				}else {
					optimal_ret = affine_gap_align(loop1.get(i).left, loop1.get(i).right, loop2.get(i).left, loop2.get(i).right);
				}
				optimal_rets.ret_nt.addAll(optimal_ret.ret_nt);
				optimal_rets.ret_score += Math.max(optimal_ret.ret_score,STAR3D.loop_penalty_limit);
			}
			count++;
		}
		return optimal_rets;
	}
}
