import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

import org.apache.commons.cli.*;
import org.ejml.data.*;
import org.ejml.ops.CommonOps;

public class STAR3D{
	/*
	 * The global parameters.
	 */
	public static String output_prefix="output";
	public static String output_fn = null;
	public static int pdb_output=0;
	public static double rmsd_cutoff=4.0;
	public static int min_stack_size=3;
	public static double gap_open_penalty=-5.;
	public static double gap_extend_penalty=-2.;
	public static int thread_num=8;
	public static double match_score=3.;
	public static double mismatch_score=0.;
	public static int map_num=5;
	public static long startTime = 0;
	public static long totalTime = 0;
	public static long endTime = 0;
	public static int connect_distance = 15;
	public static int candidate_after_map_num = 100;
	public static String annotate_tool = null;
	public static boolean dangle_end = true;
	public static double loop_penalty_limit = gap_open_penalty + gap_extend_penalty*14;
	public static boolean consider_all_cliques = false;
	public static boolean fix_stacks = true;
	public static int clique_time = 5;

	public static ArrayList<ResID> ResID1_list = new ArrayList<ResID>();
	public static ArrayList<ResID> ResID2_list = new ArrayList<ResID>();

	public static List<Point> coord1 = new ArrayList<>();
	public static List<Point> coord2 = new ArrayList<>();

	public static boolean debug = false;

	public static void main(String[] args)  throws Exception {

		startTime = System.nanoTime();
		totalTime = 0;
		endTime = 0;

		Options options = new Options();

		options.addOption("o", true, "Output alignment file");
		options.addOption("r", true, "RMSD cutoff");
		options.addOption("s", true, "Minimum stack size");
		options.addOption("g", true, "Gap open penalty");
		options.addOption("e", true, "Gap extension penalty");
		options.addOption("t", true, "Number of threads");
		options.addOption("h", false, "help information");
		options.addOption("m", true, "Match score");
		options.addOption("i", true, "Mismatch score");
		options.addOption("d", true, "Dangling end");
		options.addOption("c", true, "E-stack distance cutoff");
		options.addOption("n", true, "number of alignments");
		options.addOption("f", true, "fix stacks");
		options.addOption("l", true, "clique finding timeout");
		options.addOption("p", true, "the number of top alignments with PDB output");

		CommandLineParser parser = new BasicParser();
		CommandLine cmd = null;
		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.err.println(e.getMessage());
			System.exit(0);
		}

		HelpFormatter formatter = new HelpFormatter();

		if (cmd.hasOption("h") == true) {
			formatter.printHelp("java -jar LocalSTAR3D.jar <options> PDB1 chain1 PDB2 chain2", options);
			System.exit(0);
		}


		if (cmd.getArgList().size() != 4) {
			formatter.printHelp("java -jar LocalSTAR3D.jar <options> PDB1 chain1 PDB2 chain2", options);
			System.exit(0);
		}

		String PDBID1 = cmd.getArgList().get(0).toString().toLowerCase();
		String chainID1 = cmd.getArgList().get(1).toString();
		String PDBID2 = cmd.getArgList().get(2).toString().toLowerCase();
		String chainID2 = cmd.getArgList().get(3).toString();

		if (cmd.getOptionValue("o") != null) STAR3D.output_prefix = cmd.getOptionValue("o");
		if (cmd.getOptionValue("r") != null) STAR3D.rmsd_cutoff = Double.parseDouble(cmd.getOptionValue("r"));
		if (cmd.getOptionValue("s") != null) STAR3D.min_stack_size = Integer.parseInt(cmd.getOptionValue("s"));
		if (cmd.getOptionValue("g") != null) STAR3D.gap_open_penalty = Double.parseDouble(cmd.getOptionValue("g"));
		if (cmd.getOptionValue("e") != null) STAR3D.gap_extend_penalty = Double.parseDouble(cmd.getOptionValue("e"));
		if (cmd.getOptionValue("t") != null) STAR3D.thread_num = Integer.parseInt(cmd.getOptionValue("t"));
		if (cmd.getOptionValue("m") != null) STAR3D.match_score = Double.parseDouble(cmd.getOptionValue("m"));
		if (cmd.getOptionValue("i") != null) STAR3D.mismatch_score = Double.parseDouble(cmd.getOptionValue("i"));
		if (cmd.getOptionValue("c") != null) STAR3D.connect_distance = Integer.parseInt(cmd.getOptionValue("c"));
		if (cmd.getOptionValue("a") != null) STAR3D.annotate_tool = cmd.getOptionValue("a");
		if (cmd.getOptionValue("n") != null) STAR3D.map_num = Integer.parseInt(cmd.getOptionValue("n"));
		if (cmd.getOptionValue("f") != null) STAR3D.fix_stacks = Boolean.parseBoolean(cmd.getOptionValue("f"));
		if (cmd.getOptionValue("d") != null) STAR3D.dangle_end = Boolean.parseBoolean(cmd.getOptionValue("d"));
		if (cmd.getOptionValue("l") != null) STAR3D.clique_time = Integer.parseInt(cmd.getOptionValue("l"));
		if(cmd.getOptionValue("p") != null) STAR3D.pdb_output=Integer.parseInt(cmd.getOptionValue("p"));
		/*
		 * Parse MCA files
		 */

		output_fn = output_prefix + ".aln";

		String STAR3D_PATH = new File(STAR3D.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParentFile().getPath();
		String ROOT_PATH = new File(STAR3D_PATH).getParentFile().getPath();

		File PDB_DATA_PATH=new File(ROOT_PATH, "PDB");
		File DSSR_ANNO_PATH=new File(ROOT_PATH, "DSSR_annotation");
		File SI_DATA_PATH=new File(STAR3D_PATH, "STAR3D_struct_info");
//parse the PDB file
		File PDB1_fn = new File(PDB_DATA_PATH, PDBID1 + ".pdb");
		File PDB2_fn = new File(PDB_DATA_PATH, PDBID2 + ".pdb");

		File PDBx1_fn = new File(PDB_DATA_PATH, PDBID1 + ".cif");
		File PDBx2_fn = new File(PDB_DATA_PATH, PDBID2 + ".cif");


		int pdb1_type = -1, pdb2_type = -1;

		if (PDB1_fn.exists() && PDB1_fn.isFile()) pdb1_type = 0;
		else if (PDBx1_fn.exists() && PDBx1_fn.isFile()) pdb1_type = 1;

		if (PDB2_fn.exists() && PDB2_fn.isFile()) pdb2_type = 0;
		else if (PDBx2_fn.exists() && PDBx2_fn.isFile()) pdb2_type = 1;

		PDBParser PDB1_parser = null;
		if (pdb1_type == 0)
			PDB1_parser = new PDB(PDB1_fn, chainID1);
		else
			PDB1_parser = new PDBx(PDBx1_fn, chainID1);

		PDBParser PDB2_parser = null;
		if (pdb2_type == 0)
			PDB2_parser = new PDB(PDB2_fn, chainID2);
		else
			PDB2_parser = new PDBx(PDBx2_fn, chainID2);

		ArrayList<Residue> Res1_list = PDB1_parser.get_chain_res().get(chainID1);
		ArrayList<Residue> Res2_list = PDB2_parser.get_chain_res().get(chainID2);

		for (Residue R : Res1_list) ResID1_list.add(R.rid);
		for (Residue R : Res2_list) ResID2_list.add(R.rid);

		Map<Integer, Integer> ResID_to_star3d_index1 = new HashMap<>();
		Map<Integer, Integer> ResID_to_star3d_index2 = new HashMap<>();
		Lib.reverse_index(ResID1_list, ResID_to_star3d_index1);
		Lib.reverse_index(ResID2_list, ResID_to_star3d_index2);

		String seq1 = PDB1_parser.get_chain_seq(chainID1);
		String seq2 = PDB2_parser.get_chain_seq(chainID2);

		coord1 = PDB1_parser.get_chain_centroid(chainID1);
		coord2 = PDB2_parser.get_chain_centroid(chainID2);

		File anno1_file = null, anno2_file = null;
		DSSR DSSR1_parser = null, DSSR2_parser = null;
		HashSet<Pair<Integer, String>> bp1 = null, bp2 = null;

		anno1_file = new File(DSSR_ANNO_PATH, PDBID1 + ".dssr");
		anno2_file = new File(DSSR_ANNO_PATH, PDBID2 + ".dssr");

//		System.out.println(SI_DATA_PATH +"/"+ PDBID1 + ".dssr");

		DSSR1_parser = new DSSR(anno1_file);
		DSSR2_parser = new DSSR(anno2_file);
		bp1 = Lib.get_seq_bp(PDB1_parser, DSSR1_parser.basepair, chainID1);
		bp2 = Lib.get_seq_bp(PDB2_parser, DSSR2_parser.basepair, chainID2);

//get pair info for a specific chain (chainID)
		HashMap<Integer, ArrayList<Integer>> bp1_map = Lib.get_seq_pairing_map(bp1);
		HashMap<Integer, ArrayList<Integer>> bp2_map = Lib.get_seq_pairing_map(bp2);

		/*
		 * get base pairs in the npk secondary structure
		 */
		File ct1_fn = new File(SI_DATA_PATH, PDBID1 + "_" + chainID1 + ".npk.ct");
		File ct2_fn = new File(SI_DATA_PATH, PDBID2 + "_" + chainID2 + ".npk.ct");

		List<Pair<Integer, Integer>> npk_bp1 = Lib.get_ct_bp(ct1_fn);
		List<Pair<Integer, Integer>> npk_bp2 = Lib.get_ct_bp(ct2_fn);

		if(fix_stacks) {
			npk_bp1 = Lib.fix_broken_stack(npk_bp1, Res1_list);
			npk_bp2 = Lib.fix_broken_stack(npk_bp2, Res2_list);
		}
		/*
		 * use multiple threading to find ungapped stack mapping
		 */
		ExecutorService exec = Executors.newFixedThreadPool(thread_num);
		List<Future<List<Stackmap>>> Stackmap_multi_thread_ret = new ArrayList<Future<List<Stackmap>>>();

		for (int i = 0; i < thread_num; i++)
			Stackmap_multi_thread_ret.add(exec.submit(new FindStackmap(i, thread_num, npk_bp1, npk_bp2, coord1, coord2)));

		List<Stackmap> SM_all = new ArrayList<Stackmap>();
		for (Future<List<Stackmap>> fs : Stackmap_multi_thread_ret) {
			try {
				SM_all.addAll(fs.get());
			} catch (InterruptedException e) {
				System.out.println(e);
			} catch (ExecutionException e) {
				System.out.println(e);
			} finally {
				exec.shutdown();
			}
		}

		if (SM_all.size() == 0) {
			System.out.println("No similar stacks. You can try to reduce min_stack_size.");
			System.exit(0);
		}

		Collections.sort(SM_all);

		List<Stackmap> SM_top = SM_all.subList(0, Math.min(SM_all.size(), 10000));
//
//		for(Stackmap sm : SM_top)
//			System.out.println("sm	"+sm);

		/*
		 * find stack map configuration
		 */
		exec = Executors.newFixedThreadPool(thread_num);
		List<Future<List<Pair<Integer, Integer>>>> SMC_multi_thread_ret = new ArrayList<Future<List<Pair<Integer, Integer>>>>();

		for (int i = 0; i < thread_num; i++)
			SMC_multi_thread_ret.add(exec.submit(new FindCompatibleSM(i, thread_num, SM_top, coord1, coord2)));

		List<Pair<Integer, Integer>> SMC_all = new ArrayList<Pair<Integer, Integer>>();
		for (Future<List<Pair<Integer, Integer>>> fs : SMC_multi_thread_ret) {
			try {
				SMC_all.addAll(fs.get());
			} catch (InterruptedException e) {
				System.out.println(e);
			} catch (ExecutionException e) {
				System.out.println(e);
			} finally {
				exec.shutdown();
			}
		}

		int[][] SMC_graph = new int[SM_top.size()][SM_top.size()];
		for (Pair<Integer, Integer> P : SMC_all) {
			if (P.v != 0) {
				SMC_graph[P.left][P.right] = P.v;
				SMC_graph[P.right][P.left] = P.v;
			}
		}

		Set<Integer> SMC_graph_vertex = new HashSet<Integer>();
		for (int i = 0; i < SMC_graph.length; i++)
			SMC_graph_vertex.add(i);

		int[][] SM_connected_graph = new int[SM_top.size()][SM_top.size()];
		for (Pair<Integer, Integer> P : SMC_all) {
			if (P.v == 0) continue;
			Stackmap SM1 = SM_top.get(P.left);
			Stackmap SM2 = SM_top.get(P.right);

			boolean conn = false;
			float loop1, loop2;
			switch (P.v) {
				case 1:
					loop1 = SM2.i1 - SM1.j1;
					loop2 = SM2.i2 - SM1.j2;
					if (loop1 <= connect_distance && loop2 <= connect_distance)
						conn = true;
					break;
				case 2:
					loop1 = SM1.i1 - SM2.j1;
					loop2 = SM1.i2 - SM2.j2;
					if (loop1 <= connect_distance && loop2 <= connect_distance)
						conn = true;
					break;
				case 3:
					//1 encloseing 2 left connect
					if ((SM2.i1 - SM1.i1 - SM1.size <= connect_distance && SM2.i2 - SM1.i2 - SM1.size <= connect_distance)
//					//1 encloseing 2 right connect
							|| (SM1.j1 - SM1.size - SM2.j1 <= connect_distance && SM1.j2 - SM1.size - SM2.j2 <= connect_distance))
						conn = true;
					break;
				case 4:
					if ((SM1.i1 - SM2.i1 - SM2.size <= connect_distance && SM1.i2 - SM2.i2 - SM2.size <= connect_distance)
//					//2 encloseing 1 right connect
							|| (SM2.j1 - SM2.size - SM1.j1 <= connect_distance && SM2.j2 - SM2.size - SM1.j2 <= connect_distance))
						conn = true;
					break;
			}
			if (conn == true) {
				SM_connected_graph[P.left][P.right] = 1;
				SM_connected_graph[P.right][P.left] = 1;
			}
		}

//		System.out.println("Graph is built.");

		PriorityQueue<Aln> Alns = new PriorityQueue<Aln>(map_num, Collections.reverseOrder());
//		Set<Integer> Alns_index1 = new HashSet<Integer>();
//		Set<Integer> Alns_index2 = new HashSet<Integer>();
		int largest_clique_stack_num = 0;
		int i_map_num = 0;
		boolean clique_get_much_smaller = false;
		boolean clique_timeout = false;
		List<Set<Integer>> SM_clique_all=new ArrayList<Set<Integer>>();
		List<Set<Integer>> SM_clique=new ArrayList<Set<Integer>>();
		//Lib.Bron_Kerbosch(SMC_graph, new HashSet<Integer>(), SMC_graph_vertex, new HashSet<Integer>(), SM_clique);

		exec = Executors.newFixedThreadPool(thread_num);
		List<Future<List<Set<Integer>>>> max_cliques_multi_thread_res=new ArrayList<Future<List<Set<Integer>>>>();

        Set<Integer> tmp_vertex = new HashSet<Integer>(SMC_graph_vertex);
		max_cliques_multi_thread_res.add(exec.submit(new ConnectedBronKerbosch(SMC_graph, SM_connected_graph, tmp_vertex, SM_clique)));

		for(Future<List<Set<Integer>>> fs: max_cliques_multi_thread_res){
			try{
				SM_clique_all.addAll(fs.get(clique_time, TimeUnit.SECONDS));
			}catch(InterruptedException e){
				System.out.println(e);
			}catch(ExecutionException e) {
				System.out.println(e);
			}catch(TimeoutException e) {
				exec.shutdown();
				SM_clique_all = SM_clique;
				clique_timeout = true;
				System.out.println("Warning: Clique finding time out. LocalSTAR3D will generate alignments by using current maximum cliques.");
				System.out.println("LocalSTAR3D is not designed to align very similar RNAs. But you can still try with a larger clique_time.");
				//System.exit(0);
			}finally{
				exec.shutdown();
			}
		}

		boolean all_clique_used = false;

		while((i_map_num < map_num+candidate_after_map_num) || (Alns.size() < map_num))//  && all_clique_used==false)// && clique_get_much_smaller==false)
		{
			if(SM_clique_all.isEmpty())
				break;

		    if(debug)
                System.out.println("need aln");
			int max_clique_size = 0;
			List<Set<Integer>> max_cliques = new ArrayList<Set<Integer>>();
			for (Set<Integer> C : SM_clique_all) {
				int BP_size = 0;
				for (Integer I : C) {
					BP_size += SM_top.get(I).size;
				}
				if (BP_size == max_clique_size) max_cliques.add(C);
				if (BP_size > max_clique_size) {
					max_clique_size = BP_size;
					max_cliques.clear();
					max_cliques.add(C);
				}
			}

			/*
			 * find the mapping loops between two structures
			 */

			Set<Integer> optimal_C = new HashSet<Integer>();
			//try all the clique with max size
			for (Set<Integer> C : max_cliques) {
			    if (debug)
			        System.out.println("get new clique");
				List<Pair<Integer, Integer>> Map1_stack = new ArrayList<Pair<Integer, Integer>>();
				List<Pair<Integer, Integer>> Map2_stack = new ArrayList<Pair<Integer, Integer>>();
				//get the stack mapping: Map1_stack, Map2_stack
				for (Integer I : C) {
					Stackmap SM = SM_top.get(I);
					Map1_stack.add(new Pair<Integer, Integer>(SM.i1, SM.j1, SM.size));
					Map2_stack.add(new Pair<Integer, Integer>(SM.i2, SM.j2, SM.size));
				}

				Collections.sort(Map1_stack);
				Collections.sort(Map2_stack);
				//get the stack tree
				Snode root1 = Snode.stack_tree(Map1_stack, seq1.length());
				Snode root2 = Snode.stack_tree(Map2_stack, seq2.length());

				List<Pair<Integer, Integer>> Map1_loop = new ArrayList<Pair<Integer, Integer>>();
				List<Pair<Integer, Integer>> Map2_loop = new ArrayList<Pair<Integer, Integer>>();
				//get the loop mapping from the tree
				Snode.order_loop(root1, Map1_loop, seq1.length());
				Snode.order_loop(root2, Map2_loop, seq2.length());

				//get the mapping between stack residues
				List<Integer> Map1_stack_index = new ArrayList<Integer>();
				List<Integer> Map2_stack_index = new ArrayList<Integer>();
				for (Pair<Integer, Integer> P : Map1_stack) {
					for (int i = P.left; i < P.left + P.v; i++) Map1_stack_index.add(i);
					for (int i = P.right; i > P.right - P.v; i--) Map1_stack_index.add(i);
				}
				for (Pair<Integer, Integer> P : Map2_stack) {
					for (int i = P.left; i < P.left + P.v; i++) Map2_stack_index.add(i);
					for (int i = P.right; i > P.right - P.v; i--) Map2_stack_index.add(i);
				}

				Collections.sort(Map1_stack_index);
				Collections.sort(Map2_stack_index);

				//get the transition and rotate matrix for the stacking mapping
				DenseMatrix64F Map1_stack_X = new DenseMatrix64F(Map1_stack_index.size(), 3);
				DenseMatrix64F Map2_stack_Y = new DenseMatrix64F(Map2_stack_index.size(), 3);

				for (int i = 0; i < Map1_stack_index.size(); i++) {

					Map1_stack_X.set(i, 0, coord1.get(Map1_stack_index.get(i)).x);
					Map1_stack_X.set(i, 1, coord1.get(Map1_stack_index.get(i)).y);
					Map1_stack_X.set(i, 2, coord1.get(Map1_stack_index.get(i)).z);

					Map2_stack_Y.set(i, 0, coord2.get(Map2_stack_index.get(i)).x);
					Map2_stack_Y.set(i, 1, coord2.get(Map2_stack_index.get(i)).y);
					Map2_stack_Y.set(i, 2, coord2.get(Map2_stack_index.get(i)).z);
				}

				Point stack_XC = Geom.centroid(Map1_stack_X);
				Point stack_YC = Geom.centroid(Map2_stack_Y);

				DenseMatrix64F stack_R = new DenseMatrix64F();
				stack_R = Geom.Kabsch(Geom.translation(Map1_stack_X, stack_XC), Geom.translation(Map2_stack_Y, stack_YC));    //stack rotation matrix

				//get the mapping between the loop residues
				exec = Executors.newFixedThreadPool(thread_num);
				List<Future<LoopAlnRes>> Align_multi_thread_ret = new ArrayList<Future<LoopAlnRes>>();

				for (int i = 0; i < thread_num; i++)
					Align_multi_thread_ret.add(exec.submit(new LoopAlign(i, thread_num, Map1_loop, Map2_loop, bp1_map, bp2_map, coord1, coord2, stack_R, stack_XC, stack_YC)));

				List<Pair<Integer, Integer>> Map_loop_align = new ArrayList<Pair<Integer, Integer>>();
				double Map_loop_score = 0.;
				for (Future<LoopAlnRes> fs : Align_multi_thread_ret) {
					try {
						if (fs.get().ret_nt != null) {
							Map_loop_align.addAll(fs.get().ret_nt);
							Map_loop_score += fs.get().ret_score;
						}
					} catch (InterruptedException e) {
						System.out.println(e);
					} catch (ExecutionException e) {
						System.out.println(e);
					} finally {
						exec.shutdown();
					}
				}

				List<Integer> Map1_loop_index = new ArrayList<Integer>();
				List<Integer> Map2_loop_index = new ArrayList<Integer>();
				for (Pair<Integer, Integer> P : Map_loop_align) {
					if (P.left != -1 && P.right != -1) {
						Map1_loop_index.add(P.left);
						Map2_loop_index.add(P.right);
					}
				}

				int aln_stack_nt_number = Map1_stack_index.size();

                double map_stack_score = match_score*aln_stack_nt_number;
				double comb_score = Map_loop_score + map_stack_score;

				//combine the mapped residues
				ArrayList<Integer> Map1_index = new ArrayList<Integer>(Map1_stack_index);
				Map1_index.addAll(Map1_loop_index);
				ArrayList<Integer> Map2_index = new ArrayList<Integer>(Map2_stack_index);
				Map2_index.addAll(Map2_loop_index);

				DenseMatrix64F Map1_X = new DenseMatrix64F(Map1_index.size(), 3);
				DenseMatrix64F Map2_Y = new DenseMatrix64F(Map2_index.size(), 3);

				Collections.sort(Map1_index);
				Collections.sort(Map2_index);

				for (int i = 0; i < Map1_index.size(); i++) {
					Map1_X.set(i, 0, coord1.get(Map1_index.get(i)).x);
					Map1_X.set(i, 1, coord1.get(Map1_index.get(i)).y);
					Map1_X.set(i, 2, coord1.get(Map1_index.get(i)).z);

					Map2_Y.set(i, 0, coord2.get(Map2_index.get(i)).x);
					Map2_Y.set(i, 1, coord2.get(Map2_index.get(i)).y);
					Map2_Y.set(i, 2, coord2.get(Map2_index.get(i)).z);
				}

				double Map_rmsd = Geom.superimpose(Map1_X, Map2_Y);

				if(Map_rmsd < rmsd_cutoff) {
				    if (debug)
				        System.out.println(Map_rmsd);
					boolean both_map_overlap = false;
					boolean no_better_than_existing = false;
					Aln new_aln = new Aln(Map1_index, Map2_index, Map1_loop_index, Map2_loop_index,
							Map1_stack_index, Map2_stack_index, Map_rmsd, comb_score, Map_loop_score, map_stack_score, C);
					i_map_num++;

                    List<Set<Integer>> clique_need_to_rm = new ArrayList<Set<Integer>>();
                    for (Integer I : C) {
                        for (Set<Integer> i_C : SM_clique_all) {
                            if (i_C.contains(I)) {
                                clique_need_to_rm.add(i_C);
                                continue;
                            }
                        }
                    }
                    SM_clique_all.removeAll(clique_need_to_rm);
                    for (Integer I : C) {
                        SMC_graph_vertex.remove(I);
                    }
					if (Alns.size() == map_num && new_aln.compareTo(Alns.peek()) > 0) {
                        // the new aln is worse than the worst old aln
//                        no_better_than_existing = true;
					} else {
                        if (debug)
                            System.out.println("find new or better solution");
						List<Aln> overlapped_alns = new ArrayList<Aln>();
						for (Aln old_aln : Alns) {
							boolean map1_overlap = false;
							boolean map2_overlap = false;
							for (Integer i : Map1_index) {
								if (old_aln.rna1_index.contains(i)) {
									map1_overlap = true;
									break;
								}
							}
							for (Integer i : Map2_index) {
								if (old_aln.rna2_index.contains(i)) {
									map2_overlap = true;
									break;
								}
							}
							if (map1_overlap == true && map2_overlap == true) {
								both_map_overlap = true;
								if (new_aln.compareTo(old_aln) > 0) {
									no_better_than_existing = true;
									break;
								}
								overlapped_alns.add(old_aln);
							}
						}

						if (no_better_than_existing) {

						} else {
							if (both_map_overlap) {
								for (Aln aln : overlapped_alns) {
									if (new_aln.compareTo(aln) < 0)
										Alns.remove(aln);
								}
							} else {
                                if (Alns.size() == map_num && new_aln.compareTo(Alns.peek()) < 0)
                                    Alns.poll();
                            }
							Alns.add(new_aln);
						}
					}
				}
                SM_clique_all.remove(C);
				if(SM_clique_all.size() == 0) {
					if (debug) {
						System.out.println(SMC_graph_vertex.size());
						System.out.println("all clique are used.");
					}
					all_clique_used = true;
					break;
				}
			}

			if((all_clique_used && !clique_timeout) || SMC_graph_vertex.isEmpty())
                break;

			if(all_clique_used && clique_timeout && !SMC_graph_vertex.isEmpty()) {
			    if (debug)
                    System.out.println("search clique");
                exec = Executors.newFixedThreadPool(thread_num);
                max_cliques_multi_thread_res = new ArrayList<Future<List<Set<Integer>>>>();
                tmp_vertex = new HashSet<Integer>(SMC_graph_vertex);
                max_cliques_multi_thread_res.add(exec.submit(new ConnectedBronKerbosch(SMC_graph, SM_connected_graph, tmp_vertex, SM_clique)));
                SM_clique_all=new ArrayList<Set<Integer>>();
                SM_clique=new ArrayList<Set<Integer>>();
                for (Future<List<Set<Integer>>> fs : max_cliques_multi_thread_res) {
                    try {
                        clique_timeout = false;
                        SM_clique_all.addAll(fs.get(clique_time, TimeUnit.SECONDS));
                    } catch (InterruptedException e) {
                        System.out.println(e);
                    } catch (ExecutionException e) {
                        System.out.println(e);
                    } catch (TimeoutException e) {
                        exec.shutdown();
                        //SM_clique_all.addAll(SM_clique);
                        SM_clique_all = SM_clique;
                        clique_timeout = true;
                    } finally {
                        exec.shutdown();
                    }
                }
            }
		}
        if(debug)
		    System.out.println(Alns.size());

		List<Aln> aln_list = new ArrayList<Aln>();
		PriorityQueue<Aln> Alns_for_pdb_print = new PriorityQueue<Aln>(map_num);
		Alns_for_pdb_print.addAll(Alns);
		while(!Alns.isEmpty()) {
			aln_list.add(Alns.poll());
		}
		PrintWriter out=new PrintWriter(new BufferedWriter(new FileWriter(STAR3D.output_fn)));

		DateFormat dateFormat=new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Calendar cal=Calendar.getInstance();
		out.println(String.format("#LocalSTAR3D alignment for %s and %s (%s)", PDBID1+"_"+chainID1, PDBID2+"_"+chainID2, dateFormat.format(cal.getTime())));
		out.println("############");
		out.println("#Parameters#");
		out.println("############");
		out.println(String.format("#RMSD cutoff: %.1fA", rmsd_cutoff));
		out.println(String.format("#Minimum stack size: %d", min_stack_size));
		out.println(String.format("#Gap open penalty: %.1f", gap_open_penalty));
		out.println(String.format("#Gap extension penalty: %.1f", gap_extend_penalty));
		out.println(String.format("#Match score: %.1f", match_score));
		out.println(String.format("#Mismatch score: %.1f", mismatch_score));
		out.println(String.format("#Number of alignments: %d", map_num));
		out.println("#########");
		out.println("#Results#");
		out.println("#########");
		for(int i=aln_list.size()-1; i>=0; i--){
			Aln cur_aln = aln_list.get(i);
			out.println(String.format("#Alignment: %d", aln_list.size()-i));
			out.println(String.format("#Alignment score: %.2f", cur_aln.score));
			out.println(String.format("#Aligned nucleotide: %d", cur_aln.rna1_index.size()));
			out.println(String.format("#Alignment RMSD: %.2fA", cur_aln.rmsd));
//			for(Integer i1 : cur_aln.aligned_stack)
//				out.println(SM_top.get(i1));
			out.println("#Nucleotide mapping:");
			for(int nt_index =0; nt_index < cur_aln.rna1_index.size();nt_index++) {
				out.println(ResID1_list.get(cur_aln.rna1_index.get(nt_index))+"<->"+ResID2_list.get(cur_aln.rna2_index.get(nt_index)));
			}
			out.println("#########");
		}

		endTime = System.nanoTime();
		totalTime = (endTime - startTime)/ 1000000;
		out.println(String.format("Total time: %d.%d s", totalTime/1000, totalTime%1000));
		out.close();
		if(debug)
		    System.out.println("Done!");

		if(STAR3D.pdb_output >0) {
			for (int i = 0; i< STAR3D.pdb_output; i++) {
				//get the transition and rotate matrix for the stacking mapping
				Point Map1_XC = new Point(-1, -1, -1);
				Point Map2_YC = new Point(-1, -1, -1);
				DenseMatrix64F Map_R = new DenseMatrix64F(3,3);

				Aln cur_aln = Alns_for_pdb_print.poll();
				Lib.get_matrices(cur_aln.rna1_index, cur_aln.rna2_index, Map1_XC, Map2_YC, Map_R);

				PrintWriter PDBout = new PrintWriter(new BufferedWriter(new FileWriter(STAR3D.output_fn + Integer.toString(i+1) + ".pdb")));
				DenseMatrix64F Atom_coord = new DenseMatrix64F(1, 3);
				DenseMatrix64F Atom_coord_T = new DenseMatrix64F(3, 1);
				DenseMatrix64F Atom_coord_R = new DenseMatrix64F(3, 1);

				PDBout.println("MODEL        1");
				for (Atom A : PDB1_parser.get_chain_atom().get(chainID1)) {
					if(!cur_aln.rna1_index.contains(ResID_to_star3d_index1.get(A.res.rid.seqnum)))
						continue;
					Atom_coord.set(0, 0, A.coord.x);
					Atom_coord.set(0, 1, A.coord.y);
					Atom_coord.set(0, 2, A.coord.z);
					Atom_coord = Geom.translation(Atom_coord, Map1_XC);

					CommonOps.transpose(Atom_coord, Atom_coord_T);
					CommonOps.mult(Map_R, Atom_coord_T, Atom_coord_R);

					PDBout.println("ATOM " + String.format("%6d", A.sn) + " " + String.format("%-4s", A.atom) + " " +
							String.format("%3s", A.res.symbol) + " " + chainID1 + String.format("%4d", A.res.rid.seqnum) +
							A.res.rid.icode + "   " + String.format("%8.3f", Atom_coord_R.get(0, 0)) +
							String.format("%8.3f", Atom_coord_R.get(1, 0)) +
							String.format("%8.3f", Atom_coord_R.get(2, 0)) +
							String.format("%6.2f", A.occupancy) +
							String.format("%6.2f", A.B_iso_or_equiv)
					);
				}
				PDBout.println("ENDMDL");

				PDBout.println("MODEL        2");
				for (Atom A : PDB2_parser.get_chain_atom().get(chainID2)) {
					if(!cur_aln.rna2_index.contains(ResID_to_star3d_index2.get(A.res.rid.seqnum)))
						continue;
					Atom_coord.set(0, 0, A.coord.x);
					Atom_coord.set(0, 1, A.coord.y);
					Atom_coord.set(0, 2, A.coord.z);
					Atom_coord = Geom.translation(Atom_coord, Map2_YC);

					PDBout.println("ATOM " + String.format("%6d", A.sn) + " " + String.format("%-4s", A.atom) + " " +
							String.format("%3s", A.res.symbol) + " " + chainID2 + String.format("%4d", A.res.rid.seqnum) +
							A.res.rid.icode + "   " + String.format("%8.3f", Atom_coord.get(0, 0)) +
							String.format("%8.3f", Atom_coord.get(0, 1)) +
							String.format("%8.3f", Atom_coord.get(0, 2)) +
							String.format("%6.2f", A.occupancy) +
							String.format("%6.2f", A.B_iso_or_equiv)
					);
				}

				PDBout.println("ENDMDL");
				PDBout.close();
			}
		}

		System.exit(0);
	}
}
