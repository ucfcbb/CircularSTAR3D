import java.io.*;
import java.util.*;
import org.ejml.data.DenseMatrix64F;

public class Lib {
	//get the base pair (index based) on a sequence coming from a PDB file
    public static HashSet<Pair<Integer, String>> get_seq_bp(PDBParser pdb, HashSet<Pair<ResID, String>> basepair, String chainID){
        ArrayList<ResID> rid_list=new ArrayList<ResID>();
        for(Residue R: pdb.get_chain_res().get(chainID)) {
            rid_list.add(R.rid);
            //System.out.println("R=" + R);
        }

        HashSet<Pair<Integer, String>> bp=new HashSet<Pair<Integer, String>>();

        int i=0, j=0;
        for(Pair<ResID, String> P: basepair){
            //System.out.println("BP=" + P);
            if(P.left.chainID.equals(chainID) && P.right.chainID.equals(chainID)){
                i=rid_list.indexOf(P.left);
                j=rid_list.indexOf(P.right);
                if(i != -1 && j != -1)
                    bp.add(new Pair<Integer, String>(i, j, P.v));
            }
        }
        return bp;
    }
	
	//get the pairing partners of bases on a sequence coming from a PDB file
	public static HashMap<Integer, ArrayList<Integer>> get_seq_pairing_map(HashSet<Pair<Integer, String>> bp){
		HashMap<Integer, ArrayList<Integer>> pairing_map=new HashMap<Integer, ArrayList<Integer>>();
		for(Pair<Integer, String> P: bp){
			if(! "WHS".contains(P.v.substring(0,1)) || ! "WHS".contains(P.v.substring(1,2))) continue;
			if(pairing_map.get(P.left)==null) pairing_map.put(P.left, new ArrayList<Integer>());
			if(pairing_map.get(P.right)==null) pairing_map.put(P.right, new ArrayList<Integer>());
			pairing_map.get(P.left).add(P.right);
			pairing_map.get(P.right).add(P.left);
		}
		return pairing_map;
	}	
	
	public static List<Pair<Integer, Integer>> get_ct_bp(File ct_fn) throws IOException{
		List<Pair<Integer, Integer>> BP=new ArrayList<Pair<Integer, Integer>>();
		BufferedReader in=new BufferedReader(new FileReader(ct_fn));
		String line;
		String[] token;
		
		in.readLine();
		int i, j;
		while((line=in.readLine())!=null){
			token=line.trim().split("\\s+");
			i=Integer.valueOf(token[0]);
			j=Integer.valueOf(token[4]);
			if(i<j) BP.add(new Pair<Integer, Integer>(i-1, j-1, 1));
		}
		in.close();
		return BP;
	}
	
	//find the neighbor of i in graph G
	private static Set<Integer> GN(int[][] G, int i){
		Set<Integer> ret=new HashSet<Integer>();
		for(int j=0; j<G[i].length; j++){
			//if(G[i][j]==1 && j!=i) ret.add(j);
			if(G[i][j]!=0 && j!=i) ret.add(j);
		}
		return ret;
	}
	
	public static void Bron_Kerbosch(int[][] G, Set<Integer> R, Set<Integer> P, Set<Integer> X, List<Set<Integer>> ret){
		if(P.isEmpty() && X.isEmpty()){
			ret.add(R);
			return;
		}
		Set<Integer> PX = new HashSet<Integer>(P);
		PX.addAll(X);

		int i=new Random().nextInt(PX.size());
		int j=0, u=0;
		for(Integer I: PX){
			if(j==i) u=I;
			j++;
		}

		Set<Integer> PGNU=new HashSet<Integer>(P);

		PGNU.removeAll(GN(G, u));
		for(Integer v: PGNU){
			Set<Integer> RV=new HashSet<Integer>(R);
			RV.add((Integer)v);
			
			Set<Integer> PGNV=new HashSet<Integer>();
			for(Integer I: GN(G, v))
				if(P.contains(I)) PGNV.add(I);
			
			Set<Integer> XGNV=new HashSet<Integer>();
			for(Integer I: GN(G, v))
				if(X.contains(I)) XGNV.add(I);
			
			Bron_Kerbosch(G, RV, PGNV, XGNV, ret);
			P.remove((Integer)v);
			X.add((Integer)v);
		}
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

	//		System.out.print("RV");
	//		System.out.println(RV);
			Set<Integer> PINV = new HashSet<Integer>();
			for (Integer I : GN(G, v))
				if (P.contains(I)) PINV.add(I);

			Set<Integer> XINV = new HashSet<Integer>();
			for (Integer I : GN(G, v))
				if (X.contains(I) /*&& GN(GC,v).contains(I)*/) XINV.add(I);

	//		System.out.print("XINV");
	//		System.out.println(XINV);

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

	//		System.out.print("RV");
	//		System.out.println(RV);
			Set<Integer> PINV = new HashSet<Integer>();
			for (Integer I : GN(G, v))
				if (P.contains(I)) PINV.add(I);

			Set<Integer> XINV = new HashSet<Integer>();
			for (Integer I : GN(G, v))
				if (X.contains(I) /*&& GN(GC,v).contains(I)*/) XINV.add(I);

	//		System.out.print("XINV");
	//		System.out.println(XINV);

			connected_Bron_Kerbosch_helper(G, GC, RV, PINV, XINV, QINV, ret);
			P.remove((Integer) v);
			X.add((Integer) v);
		}
	}

    public static List<Pair<Integer, Integer>> fix_broken_stack_v0(List<Pair<Integer, Integer>> npk_bp, ArrayList<Residue> Res_list){
        HashMap<String, Integer> symbol_to_num= new HashMap<String, Integer>();
        symbol_to_num.put("A",1);
        symbol_to_num.put("C",2);
        symbol_to_num.put("G",3);
        symbol_to_num.put("U",4);

        List<Pair<Integer, Integer>> to_add = new ArrayList<Pair<Integer, Integer>>();
        for (int i = 0; i < npk_bp.size()-1; i++){
            Pair<Integer, Integer> first = npk_bp.get(i);
            Pair<Integer, Integer> second = npk_bp.get(i+1);
            if(second.left - first.left != first.right - second.right || second.left - first.left>=4)
                continue;

            boolean need_add = true;
            for(int j = 1; j<second.left - first.left; j++){
                int bp_sum = symbol_to_num.get(Res_list.get(first.left+j).symbol)+ symbol_to_num.get(Res_list.get(first.right-j).symbol);
                if( bp_sum!= 5 && bp_sum!=7) {
                    need_add = false;
                    break;
                }
            }
            if(need_add){
                for(int j = 1; j<second.left - first.left; j++){
                    to_add.add(new Pair<Integer, Integer>(first.left+j, first.right-j, 1));
                }
            }
        }

        npk_bp.addAll(to_add);
        Collections.sort(npk_bp);
        return npk_bp;
    }

    //v1
    public static List<Pair<Integer, Integer>> fix_broken_stack(List<Pair<Integer, Integer>> npk_bp, ArrayList<Residue> Res_list){
        HashMap<String, Integer> symbol_to_num= new HashMap<String, Integer>();
        symbol_to_num.put("A",1);
        symbol_to_num.put("C",2);
        symbol_to_num.put("G",3);
        symbol_to_num.put("U",4);

        HashMap<Integer, Pair<Integer, Integer>> to_add = new HashMap<Integer, Pair<Integer, Integer>>();
        for (int i = 0; i < npk_bp.size(); i++){
            Pair<Integer, Integer> cur = npk_bp.get(i);
            while((i==0 && cur.left>1 && cur.right+1< Res_list.size()) || (i>0 && cur.left>npk_bp.get(i-1).left+1 && cur.right<npk_bp.get(i-1).right-1 && cur.right+1< Res_list.size())){
                int prev_left = cur.left - 1;
                int prev_right = cur.right + 1;

                if(!symbol_to_num.containsKey(Res_list.get(prev_left).symbol) || !symbol_to_num.containsKey(Res_list.get(prev_right).symbol))
                    break;

                int bp_sum = symbol_to_num.get(Res_list.get(prev_left).symbol)+ symbol_to_num.get(Res_list.get(prev_right).symbol);
                if( bp_sum!= 5 && bp_sum!=7)
                    break;
                cur = new Pair<Integer, Integer>(prev_left, prev_right, 1);
                if(to_add.containsKey(prev_left))
                    break;
                to_add.put(prev_left, cur);
            }

            cur = npk_bp.get(i);
            while((i==npk_bp.size()-1 && cur.right-cur.left>2) || (i<npk_bp.size()-1 && cur.left+1<npk_bp.get(i+1).left && cur.right-cur.left>2)){
                int next_left = cur.left + 1;
                int next_right = cur.right - 1;

//                if(next_left > next_right)
//                    System.out.println("next " + next_left +">"+next_right);

                if(!symbol_to_num.containsKey(Res_list.get(next_left).symbol) || !symbol_to_num.containsKey(Res_list.get(next_right).symbol))
                    break;
                int bp_sum = symbol_to_num.get(Res_list.get(next_left).symbol)+symbol_to_num.get(Res_list.get(next_right).symbol);
                if( bp_sum!= 5 && bp_sum!=7)
                    break;
                cur = new Pair<Integer, Integer>(next_left, next_right, 1);
                if(to_add.containsKey(next_left))
                    break;
                to_add.put(next_left, cur);
            }
        }

        npk_bp.addAll(to_add.values());
        Collections.sort(npk_bp);
        return npk_bp;
    }

	public static void reverse_index(List<ResID> ResID, Map<Integer, Integer> ResID_to_star3d_index){
		for(Integer i=0;i<ResID.size();i++){
			ResID_to_star3d_index.put(ResID.get(i).seqnum, i);
		}
	}

	public static void get_matrices(
			List<Integer> Map1_index,
			List<Integer> Map2_index,
			Point stack_XC,
			Point stack_YC,
			DenseMatrix64F stack_R
	){
		DenseMatrix64F Map1_stack_X = new DenseMatrix64F(Map1_index.size(), 3);
		DenseMatrix64F Map2_stack_Y = new DenseMatrix64F(Map2_index.size(), 3);

		for (int i = 0; i < Map1_index.size(); i++) {
			Map1_stack_X.set(i, 0, STAR3D.coord1.get(Map1_index.get(i)).x);
			Map1_stack_X.set(i, 1, STAR3D.coord1.get(Map1_index.get(i)).y);
			Map1_stack_X.set(i, 2, STAR3D.coord1.get(Map1_index.get(i)).z);

			Map2_stack_Y.set(i, 0, STAR3D.coord2.get(Map2_index.get(i)).x);
			Map2_stack_Y.set(i, 1, STAR3D.coord2.get(Map2_index.get(i)).y);
			Map2_stack_Y.set(i, 2, STAR3D.coord2.get(Map2_index.get(i)).z);
		}

		Point _XC = Geom.centroid(Map1_stack_X);
		stack_XC.x = _XC.x;
		stack_XC.y = _XC.y;
		stack_XC.z = _XC.z;

		Point _YC = Geom.centroid(Map2_stack_Y);
		stack_YC.x = _YC.x;
		stack_YC.y = _YC.y;
		stack_YC.z = _YC.z;

		DenseMatrix64F _R = Geom.Kabsch(Geom.translation(Map1_stack_X, stack_XC), Geom.translation(Map2_stack_Y, stack_YC));
		stack_R.set(_R.numRows, _R.numCols, true, _R.data);
	}
}

