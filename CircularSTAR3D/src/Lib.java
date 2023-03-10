import org.ejml.data.DenseMatrix64F;

import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.*;


public class Lib {
	//get the base pair (index based) on a sequence coming from a PDB file
	public static HashSet<Pair<Integer, String>> get_seq_bp(PDBParser pdb, HashSet<Pair<ResID, String>> basepair, String chainID){
		ArrayList<ResID> rid_list=new ArrayList<ResID>();
		for(Residue R: pdb.get_chain_res().get(chainID)) {
			rid_list.add(R.rid);
		}

		HashSet<Pair<Integer, String>> bp=new HashSet<Pair<Integer, String>>();

		int i=0, j=0;
		for(Pair<ResID, String> P: basepair){
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

	public static List<Pair<Integer, Integer>> get_ct_bp(File ct_fn, int which) throws IOException{
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
			if(which == 2) {
				if (i < j) BP.add(new Pair<Integer, Integer>(i - 1, j - 1, 1));
			} else {
				if (0 < j)
					BP.add(new Pair<Integer, Integer>(i - 1, j - 1, 1)); //for rotate sm in cLocalSTAR3D
			}
		}
		in.close();
		return BP;
	}

	public static boolean rotate_match(Integer v, Set<Integer> R){
		if(R.isEmpty())
			return true;

		List<Integer> stk1  = new ArrayList<Integer>();
		Map<Integer, Integer> stk1_stk2  = new HashMap<Integer, Integer>();
		List<Integer> stk2  = new ArrayList<Integer>();
		for(Integer i : R){
			Stackmap sm = STAR3D.SM_top.get(i);
			stk1.add(sm.i1);stk1.add(sm.j1);
			stk1_stk2.put(sm.i1, sm.i2);stk1_stk2.put(sm.j1, sm.j2);
		}
		Collections.sort(stk1);

		boolean before_rev = true;
		Integer prev = -1;

		for(int i=0;i<stk1.size();i++) {
			Integer stk_to_add = stk1_stk2.get(stk1.get(i));
			if(prev > stk_to_add)
				before_rev = false;
			if(before_rev)
				stk2.add(stk_to_add-1000000);
			else
				stk2.add(stk_to_add);
			prev = stk_to_add;
		}

		Stackmap stk_insert = STAR3D.SM_top.get(v);

		int to_insert_i2 = stk_insert.i2;
		int to_insert_j2 = stk_insert.j2;

		if(before_rev == false){
			if(to_insert_i2> stk2.get(stk2.size()-1))
				to_insert_i2 -= 1000000;
			if(to_insert_j2> stk2.get(stk2.size()-1))
				to_insert_j2 -= 1000000;
		}else{
			to_insert_i2 -= 1000000;
			to_insert_j2 -= 1000000;
		}

		int pos_i1 = -(Collections.binarySearch(stk1, stk_insert.i1)+1);
		int pos_j1 = -(Collections.binarySearch(stk1, stk_insert.j1)+1);
		int pos_i2 = -(Collections.binarySearch(stk2, to_insert_i2)+1);
		int pos_j2 = -(Collections.binarySearch(stk2, to_insert_j2)+1);

		if(pos_i1 == stk1.size())	pos_i1 = 0;
		if(pos_i2 == stk1.size())	pos_i1 = 0;
		if(pos_j1 == stk1.size())	pos_j1 = 0;
		if(pos_j2 == stk1.size())	pos_j2 = 0;

		if(pos_i1 == pos_i2 && pos_j1 == pos_j2) {
			return true;
		}

		return false;
	}


	public static void Bron_Kerbosch(Map<Integer, Set<Integer>> G, Set<Integer> R, Set<Integer> P, Set<Integer> X, List<Set<Integer>> ret){
		// remove the node not rotate match in X
		Set<Integer> X_rotate_match = new HashSet<>(X);
		Set<Integer> P_rotate_match = new HashSet<>(P);
		for(Integer i : X){
			if(!rotate_match(i, R)){
				X_rotate_match.remove(i);
			}
		}
		for(Integer i : P){
			if(!rotate_match(i, R)){
				P_rotate_match.remove(i);
			}
		}

		if(P_rotate_match.isEmpty() && X_rotate_match.isEmpty()){
			ret.add(R);
			return;
		}

		Set<Integer> Pcur = new HashSet<>(P_rotate_match);

		for(Integer v: Pcur){
			if(rotate_match(v,R)) {
				Set<Integer> RV = new HashSet<Integer>(R);
				RV.add((Integer) v);

				Set<Integer> PGNV = new HashSet<Integer>(P_rotate_match);
				PGNV.retainAll(G.get(v));

				Set<Integer> XGNV = new HashSet<Integer>(X_rotate_match);
				XGNV.retainAll(G.get(v));

				Bron_Kerbosch(G, RV, PGNV, XGNV, ret);
			}
			P.remove((Integer)v);
			X.add((Integer)v);
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

	public static void extend_sm_set(Set<Integer> seed, Map<Integer, Set<Integer>> i1_sm, Map<Integer, Set<Integer>> j1_sm,
									 Map<Integer, Set<Integer>> i2_sm, Map<Integer, Set<Integer>> j2_sm,
									 Set<Integer> cur_all, Set<Integer> cur_ext,  Set<Set<Integer>> all_sm) {
		// cur_all are the stacks in current set
		// cur_ext are the current stacks that can be extended
		Set<Integer> extending = new HashSet<>(cur_all);
		// will change to Set<Set<Integer>> extending, one set for each sm
		for (Integer i_sm : cur_ext) {
			Stackmap _sm = STAR3D.SM_top.get(i_sm);

			Integer l1_end = _sm.i1 + _sm.size;
			Integer r1_end = _sm.j1 - _sm.size;
			Integer l2_end = _sm.i2 + _sm.size;
			Integer r2_end = _sm.j2 - _sm.size;

			Set<Integer> nxt_sm_1 = new HashSet<>();
			Set<Integer> nxt_sm_2 = new HashSet<>();
			Set<Integer> to_add;

			//left connected
			for (Integer tmp_i1 = l1_end; tmp_i1 < l1_end + STAR3D.connect_distance; tmp_i1++) {
				if (tmp_i1 >= r1_end || tmp_i1 >= STAR3D.ResID1_list.size()) break;
				if (i1_sm.keySet().contains(tmp_i1)) {
					nxt_sm_1.addAll(i1_sm.get(tmp_i1));
				}
			}
			//8_21
			for (Integer tmp_i1 = l1_end; tmp_i1 < l1_end + STAR3D.connect_distance; tmp_i1++) {
				if (tmp_i1 >= r1_end || tmp_i1 >= STAR3D.ResID1_list.size()) break;
				if (j1_sm.keySet().contains(tmp_i1)) {
					nxt_sm_1.addAll(j1_sm.get(tmp_i1));
				}
			}

			for (Integer tmp_i2 = l2_end; tmp_i2 < l2_end + STAR3D.connect_distance; tmp_i2++) {
				if (tmp_i2 >= r2_end || tmp_i2 >= STAR3D.ResID2_list.size()) break;
				if (i2_sm.keySet().contains(tmp_i2)) {
					nxt_sm_2.addAll(i2_sm.get(tmp_i2));
				}
			}
			//8_21
			for (Integer tmp_i2 = l2_end; tmp_i2 < l2_end + STAR3D.connect_distance; tmp_i2++) {
				if (tmp_i2 >= r2_end || tmp_i2 >= STAR3D.ResID2_list.size()) break;
				if (j2_sm.keySet().contains(tmp_i2)) {
					nxt_sm_2.addAll(j2_sm.get(tmp_i2));
				}
			}

			to_add = new HashSet<>(nxt_sm_1);
			to_add.retainAll(nxt_sm_2);

			extending.addAll(to_add);

			nxt_sm_1.clear();
			nxt_sm_2.clear();

			//right connected
			for (Integer tmp_j1 = r1_end; tmp_j1 > r1_end - STAR3D.connect_distance; tmp_j1--) {
				if (tmp_j1 <= 0 || tmp_j1 <= l1_end) break;
				if (j1_sm.keySet().contains(tmp_j1)) {
					nxt_sm_1.addAll(j1_sm.get(tmp_j1));
				}
			}

			for (Integer tmp_j1 = r1_end; tmp_j1 > r1_end - STAR3D.connect_distance; tmp_j1--) {
				if (tmp_j1 <= 0 || tmp_j1 <= l1_end) break;
				if (i1_sm.keySet().contains(tmp_j1)) {
					nxt_sm_1.addAll(i1_sm.get(tmp_j1));
				}
			}

			for (Integer tmp_j2 = r2_end; tmp_j2 > r2_end - STAR3D.connect_distance; tmp_j2--) {
				if (tmp_j2 <= 0 || tmp_j2 <= l2_end) break;
				if (j2_sm.keySet().contains(tmp_j2)) {
					nxt_sm_2.addAll(j2_sm.get(tmp_j2));
				}
			}

			for (Integer tmp_j2 = r2_end; tmp_j2 > r2_end - STAR3D.connect_distance; tmp_j2--) {
				if (tmp_j2 <= 0 || tmp_j2 <= l2_end) break;
				if (i2_sm.keySet().contains(tmp_j2)) {
					nxt_sm_2.addAll(i2_sm.get(tmp_j2));
				}
			}

			to_add = new HashSet<>(nxt_sm_1);
			to_add.retainAll(nxt_sm_2);

			extending.addAll(to_add);

		}

		Set<Integer> set_of_sm_index = new HashSet<>(cur_all);
		set_of_sm_index.addAll(extending);

		//for sm in each loop, find smc
		List<Pair<Integer, Integer>> SMC_loop = new ArrayList<Pair<Integer, Integer>>();
		ExecutorService exec = Executors.newFixedThreadPool(STAR3D.thread_num);
		List<Future<List<Pair<Integer, Integer>>>> SMC_multi_thread_ret = new ArrayList<Future<List<Pair<Integer, Integer>>>>();
		List<Integer> list_of_sm_index = new ArrayList<>(set_of_sm_index);
		for (int i = 0; i < STAR3D.thread_num; i++)
			SMC_multi_thread_ret.add(exec.submit(new FindCompatibleSM(i, STAR3D.thread_num, list_of_sm_index, STAR3D.coord1, STAR3D.coord2)));

		for (Future<List<Pair<Integer, Integer>>> fs : SMC_multi_thread_ret) {
			try {
				SMC_loop.addAll(fs.get());
			} catch (InterruptedException e) {
				System.out.println(e);
			} catch (ExecutionException e) {
				System.out.println(e);
			} finally {
				exec.shutdown();
			}
		}

		Map<Integer, Set<Integer>> SMC_graph = new HashMap<>();
		for (Pair<Integer, Integer> P : SMC_loop) {
			if (P.v != 0) {
				if(!SMC_graph.keySet().contains(P.left))
					SMC_graph.put(P.left, new HashSet<>());
				if(!SMC_graph.keySet().contains(P.right))
					SMC_graph.put(P.right, new HashSet<>());

				if(P.left!=P.right) {
					SMC_graph.get(P.left).add(P.right);
					SMC_graph.get(P.right).add(P.left);
				}
			}
		}

		Set<Integer> P_for_clique = new HashSet<Integer>(SMC_graph.keySet());
		P_for_clique.removeAll(seed);

		List<Set<Integer>> SM_clique = new ArrayList<Set<Integer>>();

		for(Integer i : seed)
			P_for_clique.retainAll(SMC_graph.get(i));

		Bron_Kerbosch(SMC_graph, seed, P_for_clique, new HashSet<Integer>(), SM_clique);

//		for (Set<Integer> C : SM_clique) {
//			if (C.size() > cur_all.size()) {
//				Set<Integer> nxt_nxt_sm = new HashSet<>(C);
//				nxt_nxt_sm.removeAll(cur_all);
//				extend_sm_set(seed, i1_sm, j1_sm, i2_sm, j2_sm, new HashSet<>(C), nxt_nxt_sm, all_sm);
//			} else { // cannot extent
//				if(C.size() >= seed.size()) // can be removed after updating the clique function
//					all_sm.add(C);
//			}
//		}

		int bp_in_cur_all = 0;
		for (Integer I : cur_all){
			bp_in_cur_all += STAR3D.SM_top.get(I).size;
		}

		Set<Set<Integer>> max_cliques_for_cur_loop = new HashSet<>();
		int max_bp = 0;
		for (Set<Integer> C : SM_clique) {
			int BP_size = 0;
			for (Integer I : C)
				BP_size += STAR3D.SM_top.get(I).size;
			if (BP_size > max_bp) {
				max_bp = BP_size;
				max_cliques_for_cur_loop.clear();
				max_cliques_for_cur_loop.add(C);
			}
		}
		if(max_bp == bp_in_cur_all){// can not extend anymore
			all_sm.add(cur_all);
//			return;
		}else{
			for(Set<Integer> C : max_cliques_for_cur_loop){
				Set<Integer> nxt_nxt_sm = new HashSet<>(C);
				nxt_nxt_sm.removeAll(cur_all);
				extend_sm_set(seed, i1_sm, j1_sm, i2_sm, j2_sm, new HashSet<>(C), nxt_nxt_sm, all_sm);
			}
		}
	}

	public static void get_bp_close_loop(DSSR DSSR_parser, HashMap<Pair<Integer, Integer>, Integer> bp_to_loop){
		for(int i=0 ;i<DSSR_parser.loops_close_bp.size();i++) {
			List<Integer> cur_loop = DSSR_parser.loops_close_bp.get(i);
			bp_to_loop.put(new Pair<Integer, Integer>(cur_loop.get(0), cur_loop.get(cur_loop.size() - 1), 0), i);
			for (int j = 1; j < cur_loop.size() - 1; j += 2) {
				bp_to_loop.put(new Pair<Integer, Integer>(cur_loop.get(j), cur_loop.get(j + 1), 0), i);
			}
		}
	}

	public static void shrink_loop(List<Pair<Integer, Integer>> npk_bp, HashMap<Pair<Integer, Integer>, Integer> bp_to_loop_raw,
								   HashMap<Pair<Integer, Integer>, Integer> bp_to_loop, Map<Integer, Integer> ResID_to_star3d_index){
		Set<Pair<Integer, Integer>> bp_set = new HashSet<>(npk_bp);
		for(Pair<Integer, Integer> p : bp_to_loop_raw.keySet()){
			if(!(ResID_to_star3d_index.containsKey(p.left) && ResID_to_star3d_index.containsKey(p.right)))
				continue;
			Integer left = ResID_to_star3d_index.get(p.left); Integer right = ResID_to_star3d_index.get(p.right);
			while(bp_set.contains(new Pair<Integer, Integer>(left+1, right-1, 1))) {
				left++;right--;
			}
			bp_to_loop.put(new Pair<Integer, Integer>(left, right, 0), bp_to_loop_raw.get(p));
		}
	}

	public static void reverse_index(List<ResID> ResID, Map<Integer, Integer> ResID_to_star3d_index){
		for(Integer i=0;i<ResID.size();i++){
			ResID_to_star3d_index.put(ResID.get(i).seqnum, i);
		}
	}

	public static void get_bp_close_loop_and_shrink(DSSR DSSR_parser,
													Map<Pair<Integer, Integer>, Integer> bp_to_loop,
													List<Pair<Integer, Integer>> npk_bp,
													Map<Integer, Integer> ResID_to_star3d_index,
													Set<Integer> seed_loop_nt){
		List<Integer> seed_loop_boundary = new ArrayList<>();

		Set<Pair<Integer, Integer>> bp_set = new HashSet<>(npk_bp);

		for(int i=0 ;i<DSSR_parser.loops_close_bp.size();i++) {
			// the outmost bp
			List<Integer> cur_loop = DSSR_parser.loops_close_bp.get(i);
			int left_index = cur_loop.get(0);
			int right_index = cur_loop.get(cur_loop.size() - 1);
			if(!(ResID_to_star3d_index.containsKey(left_index) && ResID_to_star3d_index.containsKey(right_index)))
				continue;
			int left = ResID_to_star3d_index.get(left_index);
			int right = ResID_to_star3d_index.get(right_index);

			while(bp_set.contains(new Pair<Integer, Integer>(left+1, right-1, 1))) {
				left++;
				right--;
			}
			bp_to_loop.put(new Pair<Integer, Integer>(left, right, 0), i);

			seed_loop_boundary.add(left);
			seed_loop_boundary.add(right);

			// other bp
			for (int j = 1; j < cur_loop.size() - 1; j += 2) {
				left_index = cur_loop.get(j);
				right_index = cur_loop.get(j + 1);
				if(!(ResID_to_star3d_index.containsKey(left_index) && ResID_to_star3d_index.containsKey(right_index)))
					continue;
				left = ResID_to_star3d_index.get(left_index);
				right = ResID_to_star3d_index.get(right_index);

				while(bp_set.contains(new Pair<Integer, Integer>(left-1, right+1, 1))) {
					left--;
					right++;
				}
				bp_to_loop.put(new Pair<Integer, Integer>(left, right,1), i);

				seed_loop_boundary.add(left);
				seed_loop_boundary.add(right);
			}

			Collections.sort(seed_loop_boundary);
			for(int j = 0;j < seed_loop_boundary.size();j+=2){
				for (Integer k = seed_loop_boundary.get(j)+1;k < seed_loop_boundary.get(j+1);k += 1){
					seed_loop_nt.add(k);
				}
			}
		}
	}

	public static void get_bp_list_ordered(Set<Pair<Integer,Integer>> bp_list, int i, int j, int dir){
		bp_list.add(new Pair<Integer, Integer>(i, j, 1));

		bp_list.add(new Pair<Integer, Integer>(i+dir, j-dir, 1));
		bp_list.add(new Pair<Integer, Integer>(i, j-dir, 1));
		bp_list.add(new Pair<Integer, Integer>(i+dir, j, 1));

		bp_list.add(new Pair<Integer, Integer>(i+dir, j-2*dir, 1));
		bp_list.add(new Pair<Integer, Integer>(i+2*dir, j-dir, 1));
		bp_list.add(new Pair<Integer, Integer>(i+2*dir, j-2*dir, 1));
		bp_list.add(new Pair<Integer, Integer>(i+2*dir, j, 1));
		bp_list.add(new Pair<Integer, Integer>(i, j-2*dir, 1));
	}

	public static void get_bp_list(Set<Pair<Integer,Integer>> sm_bp_list, int i, int j, int sm_size){
		if (i < j) {
			get_bp_list_ordered(sm_bp_list, i, j,-1);
			get_bp_list_ordered(sm_bp_list, i, j,1); // to rm
			i += sm_size;j -= sm_size;
			get_bp_list_ordered(sm_bp_list, i, j,1);
			get_bp_list_ordered(sm_bp_list, i, j,-1); // to rm
		} else {
			get_bp_list_ordered(sm_bp_list, j, i,1);
			get_bp_list_ordered(sm_bp_list, j, i,-1); // to rm
			i += sm_size;j -= sm_size;
			get_bp_list_ordered(sm_bp_list, j, i,-1);
			get_bp_list_ordered(sm_bp_list, j, i,1); // to rm
		}
	}

	public static void assign_sm_to_loop(List<Stackmap> SM_top, Map<Pair<Integer, Integer>, Set<Integer>>  looppair_to_sm,
										 Map<Pair<Integer, Integer>, Integer> bp_to_loop1, Map<Pair<Integer, Integer>, Integer> bp_to_loop2, int i){
		Stackmap sm = SM_top.get(i);

		int i1 = sm.i1;
		int j1 = sm.j1;
		int i2 = sm.i2;
		int j2 = sm.j2;

		Set<Pair<Integer,Integer>> sm_bp1_list = new HashSet<>();
		Set<Pair<Integer,Integer>> sm_bp2_list = new HashSet<>();

		// get the bp belongs to a sm
		get_bp_list(sm_bp1_list, i1, j1, sm.size);
		get_bp_list(sm_bp2_list, i2, j2, sm.size);

		// get sm_bp close to loop
		sm_bp1_list.retainAll(bp_to_loop1.keySet());
		sm_bp2_list.retainAll(bp_to_loop2.keySet());

		for(Pair<Integer,Integer> p1 : sm_bp1_list){
			for(Pair<Integer,Integer> p2 : sm_bp2_list){
				Pair<Integer, Integer> p1_p2 = new Pair(new Integer(bp_to_loop1.get(p1)), new Integer(bp_to_loop2.get(p2)),0);
				if(!looppair_to_sm.containsKey(p1_p2))
					looppair_to_sm.put(p1_p2, new HashSet<Integer>());
				looppair_to_sm.get(p1_p2).add(i);
			}
		}
	}

	public static boolean overlap(int s1, int l1, int s2, int l2){
		if((s2>=s1 && s2<=s1+l1-1) || (s1>=s2 && s1<=s2+l2-1))
			return true;
		else
			return false;
	}

	public static boolean overlap_sm(Stackmap sm1, Stackmap sm2){
		if(overlap(sm1.i1, sm1.size, sm2.i1, sm2.size) &&
				overlap(sm1.j1, sm1.size, sm2.j1, sm2.size) &&
				overlap(sm1.i2, sm1.size, sm2.i2, sm2.size) &&
				overlap(sm1.j2, sm1.size, sm2.j2, sm2.size))
			return true;
		return false;
	}

	public static Aln overlapped_with_old_alns(Aln new_aln, PriorityQueue<Aln> Alns){
		for (Aln old_aln : Alns) {
			for(Integer i : old_aln.aligned_stack){
				for(Integer j: new_aln.aligned_stack){
					if(overlap_sm(STAR3D.SM_top.get(i), STAR3D.SM_top.get(j))){
						return old_aln;
					}
				}
			}
		}
		return null;
	}

	/*** generate strand from sm set ****/
	public static void get_stack_stand(
			Set<Integer> C,
			List<Pair<Integer, Integer>> Map1_stack_strand,
			List<Pair<Integer, Integer>> Map2_stack_strand,
			Map<Pair<Integer, Integer>, Pair<Integer, Integer>>stk1_to_stk2
	){
		int stk_index = 0;
		for (Integer I : C) {
			Stackmap SM = STAR3D.SM_top.get(I);

			Pair<Integer, Integer> stk1_strand1 = new Pair<Integer, Integer>(SM.i1, SM.i1 + SM.size - 1, stk_index);
			Pair<Integer, Integer> stk1_strand2 = new Pair<Integer, Integer>(SM.j1 - SM.size + 1, SM.j1, stk_index);
			Pair<Integer, Integer> stk2_strand1 = new Pair<Integer, Integer>(SM.i2, SM.i2 + SM.size - 1, stk_index);
			Pair<Integer, Integer> stk2_strand2 = new Pair<Integer, Integer>(SM.j2 - SM.size + 1, SM.j2, stk_index);

			Map1_stack_strand.add(stk1_strand1);
			Map1_stack_strand.add(stk1_strand2);
			stk1_to_stk2.put(stk1_strand1, stk2_strand1);
			stk1_to_stk2.put(stk1_strand2, stk2_strand2);

			stk_index += 1;
		}

		Collections.sort(Map1_stack_strand);
		for (Pair<Integer, Integer> stk1 : Map1_stack_strand) {
			Map2_stack_strand.add(stk1_to_stk2.get(stk1));
		}
	}

	public static double run_loop_aln(
			List<Pair<Integer, Integer>> Map_loop_align,
			List<Pair<Integer, Integer>> Map1_loop,
			List<Pair<Integer, Integer>> Map2_loop,
			Map<Integer, ArrayList<Integer>> bp1_map,
			Map<Integer, ArrayList<Integer>> bp2_map,
			List<Point> coord1, List<Point> coord2,
			DenseMatrix64F stack_R,
			Point stack_XC,
			Point stack_YC,
			boolean print_flag
	){
		double Map_loop_score = 0.;

		ExecutorService exec = Executors.newFixedThreadPool(STAR3D.thread_num);
		List<Future<LoopAlnRes>> Align_multi_thread_ret = new ArrayList<Future<LoopAlnRes>>();

		for (int i = 0; i < STAR3D.thread_num; i++)
			Align_multi_thread_ret.add(exec.submit(new LoopAlign(i, STAR3D.thread_num, Map1_loop, Map2_loop, bp1_map, bp2_map, coord1, coord2, stack_R, stack_XC, stack_YC, print_flag)));

		for (Future<LoopAlnRes> fs : Align_multi_thread_ret) {
			try {
				if (fs.get().ret_nt != null) {
					Map_loop_align.addAll(fs.get().ret_nt);
					Map_loop_score += fs.get().ret_score;
				}
			} catch (InterruptedException e) {
				System.out.println(e);
				exec.shutdown();
			} catch (ExecutionException e) {
				System.out.println(e);
				exec.shutdown();
			} finally {
				exec.shutdown();
			}
		}

		return Map_loop_score;
	}

	public static void add_loop_to_map_index(
			List<Pair<Integer, Integer>> Map_loop_align,
			ArrayList<Integer> Map1_loop_index,
			ArrayList<Integer> Map2_loop_index,
			ArrayList<Integer> Map1_index,
			ArrayList<Integer> Map2_index
	){
		// nt in stacks
		Map<Integer, Integer> Map1_to_Map2_index = new HashMap<Integer, Integer>();
		for(int i = 0;i< Map1_index.size();i++)
			Map1_to_Map2_index.put(Map1_index.get(i), Map2_index.get(i));

		// add nt in loop
		for (Pair<Integer, Integer> P : Map_loop_align) {
			if (P.left != -1 && P.right != -1) {
				Map1_loop_index.add(P.left);
				Map2_loop_index.add(P.right);
				Map1_to_Map2_index.put(P.left, P.right);
			}
		}

		Map1_index.addAll(Map1_loop_index);

		Map2_index.clear();
		Collections.sort(Map1_index);

		for (Integer i : Map1_index)
			Map2_index.add(Map1_to_Map2_index.get(i));
	}

	public static void print_result(
			PriorityQueue<Aln> Alns,
			String PDBID1,
			String chainID1,
			String PDBID2,
			String chainID2
	){
		List<Aln> aln_list = new ArrayList<Aln>();
		while(!Alns.isEmpty()) {
			aln_list.add(Alns.poll());
		}

		PrintWriter out= null;
		try {
			out = new PrintWriter(new BufferedWriter(new FileWriter(STAR3D.output_fn)));
		} catch (IOException e) {
			e.printStackTrace();
		}

		DateFormat dateFormat=new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Calendar cal=Calendar.getInstance();
		out.println(String.format("#CircularSTAR3D alignment for %s and %s (%s)", PDBID1+"_"+chainID1, PDBID2+"_"+chainID2, dateFormat.format(cal.getTime())));
		out.println("############");
		out.println("#Parameters#");
		out.println("############");
		out.println(String.format("#RMSD cutoff: %.1fA", STAR3D.rmsd_cutoff));
		out.println(String.format("#Minimum stack size: %d", STAR3D.min_stack_size));
		out.println(String.format("#Gap open penalty: %.1f", STAR3D.gap_open_penalty));
		out.println(String.format("#Gap extension penalty: %.1f", STAR3D.gap_extend_penalty));
		out.println(String.format("#Match score: %.1f", STAR3D.match_score));
		out.println(String.format("#Mismatch score: %.1f", STAR3D.mismatch_score));
		out.println(String.format("#Number of alignments: %d", STAR3D.map_num));
		out.println("#########");
		out.println("#Results#");
		out.println("#########");
		for(int i=aln_list.size()-1; i>=0; i--){
			Aln cur_aln = aln_list.get(i);
			out.println(String.format("#Alignment: %d", aln_list.size()-i));
			out.println(String.format("#Alignment score: %.2f", cur_aln.score));
			out.println(String.format("#Aligned nucleotide: %d", cur_aln.rna1_index.size()));
			out.println(String.format("#Alignment RMSD: %.2fA", cur_aln.rmsd));

			out.println(String.format("#contain junction: %b", cur_aln.contain_junction));
			out.print(String.format("#Seed loop1: "));
			for(Integer in :  cur_aln.info.get(0))
				out.print(in + " ");
			out.println();
			out.print(String.format("#Seed loop2: "));
			for(Integer in :  cur_aln.info.get(1))
				out.print(in + " ");
			out.println();

			out.println("#Nucleotide mapping:");
			for(int nt_index =0; nt_index < cur_aln.rna1_index.size();nt_index++) {
				out.println(STAR3D.ResID1_list.get(cur_aln.rna1_index.get(nt_index))+"<->"+STAR3D.ResID2_list.get(cur_aln.rna2_index.get(nt_index)));
			}

//			out.println("#for pymol");
//
//			out.print("0	"+PDBID1+"	"+chainID1+"	");
//			String aln1 = "";
//			String aln2 = "";
//			for(int nt_index = 0; nt_index < cur_aln.rna1_index.size();nt_index++) {
//				aln1 += STAR3D.ResID1_list.get(cur_aln.rna1_index.get(nt_index)).seqnum;
//				if(nt_index != cur_aln.rna1_index.size()-1)
//					aln1 += ",";
//			}
//			out.println(aln1 + "\t" + aln1);
//
//			out.print("1	"+PDBID2+"	"+chainID2+"	");
//			for(int nt_index = 0; nt_index < cur_aln.rna2_index.size();nt_index++) {
//				aln2 += STAR3D.ResID2_list.get(cur_aln.rna2_index.get(nt_index)).seqnum;
//				if(nt_index != cur_aln.rna2_index.size()-1)
//					aln2 += ",";
//			}
//			out.println(aln2 + "\t" + aln2);

			out.println("#########");
		}

		STAR3D.endTime = System.nanoTime();
		STAR3D.totalTime = (STAR3D.endTime - STAR3D.startTime) / 1000000;
		System.out.println(String.format("All time: %d.%d s", STAR3D.totalTime / 1000, STAR3D.totalTime % 1000));

		out.println(String.format("Total time: %d.%d s", STAR3D.totalTime/1000, STAR3D.totalTime%1000));
		out.close();
	}

	public static Boolean get_rotation_matrix(
			List<Integer> Map1_stack_index,
			List<Integer> Map2_stack_index,
			Point stack_XC,
			Point stack_YC,
			DenseMatrix64F stack_R,
			double rmsd
	){
		DenseMatrix64F Map1_stack_X = new DenseMatrix64F(Map1_stack_index.size(), 3);
		DenseMatrix64F Map2_stack_Y = new DenseMatrix64F(Map2_stack_index.size(), 3);

		for (int i = 0; i < Map1_stack_index.size(); i++) {
			Map1_stack_X.set(i, 0, STAR3D.coord1.get(Map1_stack_index.get(i)).x);
			Map1_stack_X.set(i, 1, STAR3D.coord1.get(Map1_stack_index.get(i)).y);
			Map1_stack_X.set(i, 2, STAR3D.coord1.get(Map1_stack_index.get(i)).z);

			Map2_stack_Y.set(i, 0, STAR3D.coord2.get(Map2_stack_index.get(i)).x);
			Map2_stack_Y.set(i, 1, STAR3D.coord2.get(Map2_stack_index.get(i)).y);
			Map2_stack_Y.set(i, 2, STAR3D.coord2.get(Map2_stack_index.get(i)).z);
		}

		double stack_rmsd = Geom.superimpose(Map1_stack_X, Map2_stack_Y);

		if (rmsd >0 && stack_rmsd > rmsd) return false;

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
		return true;
	}

	public static void get_junction_id(
			DSSR DSSR_parser,
			Set<Integer> junction_id
	){
		for(int i=0;i<DSSR_parser.loops_close_bp.size();i++){
			if(DSSR_parser.loops_close_bp.get(i).size()>4){
				junction_id.add(i);
			}
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
