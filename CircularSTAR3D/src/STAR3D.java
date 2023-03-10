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
	public static int thread_num=1;
	public static double gap_open_penalty=-5.;
	public static double gap_extend_penalty=-2.;
	public static double match_score=3.;
	public static double mismatch_score=0.;
	public static double stack_discount = 1;
	public static int map_num=5;
	public static int seed_num=100;
	public static long startTime = 0;
	public static long totalTime = 0;
	public static long endTime = 0;
	public static int connect_distance = 15;
	public static String annotate_tool = null;
	public static boolean no_dangle_end = false;
	public static double loop_penalty_limit = gap_open_penalty + gap_extend_penalty*(connect_distance-1);
	public static boolean fix_stacks = true;
	public static int clique_time = 5;

	public static List<ResID> ResID1_list = new ArrayList<ResID>();
	public static List<ResID> ResID2_list = new ArrayList<ResID>();
	public static List<Stackmap> SM_top = new ArrayList<Stackmap>();
	public static List<Point> coord1 = new ArrayList<>();
	public static List<Point> coord2 = new ArrayList<>();

	public static Set<Integer> seed_loop_nt1 = new HashSet<>();
	public static Set<Integer> seed_loop_nt2 = new HashSet<>();

	public static boolean debug = false;

	public static boolean motif_mode = false;
	public static int minimum_aligned_loop_size = 0;

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
		options.addOption("d", true, "Include dangling end");
		options.addOption("c", true, "E-stack distance cutoff");
		options.addOption("n", true, "Number of alignments");
		options.addOption("f", true, "Fix stacks");
		options.addOption("l", true, "Clique finding timeout");
		options.addOption("p", true, "the number of top alignments with PDB output");
		options.addOption("motif", false, "Motif mode");

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
			formatter.printHelp("java -jar CircularSTAR3D.jar <options> PDB1 chain1 PDB2 chain2", options);
			System.exit(0);
		}

		if (cmd.getArgList().size() != 4) {
			formatter.printHelp("java -jar CircularSTAR3D.jar <options> PDB1 chain1 PDB2 chain2", options);
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
		if (cmd.getOptionValue("d") != null) STAR3D.no_dangle_end = Boolean.parseBoolean(cmd.getOptionValue("d"));
		if (cmd.getOptionValue("l") != null) STAR3D.clique_time = Integer.parseInt(cmd.getOptionValue("l"));
		if(cmd.getOptionValue("p") != null) STAR3D.pdb_output=Integer.parseInt(cmd.getOptionValue("p"));

		output_fn = output_prefix + ".aln";

		if (cmd.hasOption("d") == true) {STAR3D.no_dangle_end = true;}

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

		coord1 = PDB1_parser.get_chain_centroid(chainID1);
		coord2 = PDB2_parser.get_chain_centroid(chainID2);

		File anno1_file, anno2_file;
		DSSR DSSR1_parser, DSSR2_parser;
		HashSet<Pair<Integer, String>> bp1, bp2;

		anno1_file = new File(DSSR_ANNO_PATH, PDBID1 + ".dssr");
		anno2_file = new File(DSSR_ANNO_PATH, PDBID2 + ".dssr");

		DSSR1_parser = new DSSR(anno1_file, chainID1);
		DSSR2_parser = new DSSR(anno2_file, chainID2);
		bp1 = Lib.get_seq_bp(PDB1_parser, DSSR1_parser.basepair, chainID1);
		bp2 = Lib.get_seq_bp(PDB2_parser, DSSR2_parser.basepair, chainID2);

		//get pair info for a specific chain (chainID)
		HashMap<Integer, ArrayList<Integer>> bp1_map = Lib.get_seq_pairing_map(bp1);
		HashMap<Integer, ArrayList<Integer>> bp2_map = Lib.get_seq_pairing_map(bp2);

		/*
		 * get base pairs in the npk secondary structure
		 */
		File ct1_fn = new File(SI_DATA_PATH, PDBID1 + "_" + chainID1 + ".ct");
		File ct2_fn = new File(SI_DATA_PATH, PDBID2 + "_" + chainID2 + ".ct");

		List<Pair<Integer, Integer>> npk_bp1 = Lib.get_ct_bp(ct1_fn, 1);
		List<Pair<Integer, Integer>> npk_bp2 = Lib.get_ct_bp(ct2_fn, 1);

		if (fix_stacks) {
			npk_bp1 = Lib.fix_broken_stack(npk_bp1, Res1_list);
			npk_bp2 = Lib.fix_broken_stack(npk_bp2, Res2_list);
		}

		/*********************** use multiple threading to find ungapped stack mapping*********************************/
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
		SM_top = SM_all.subList(0, Math.min(SM_all.size(), 20000));
		System.out.println("Stack pair number: "+SM_top.size());

		/*********************** assign stack pairs to loop pairs *********************************/
		Map<Pair<Integer, Integer>, Set<Integer>> loopInRNAs_to_sm = new HashMap<>();

		Map<Pair<Integer, Integer>, Integer> bp_to_loop1 = new HashMap<>();
		Map<Pair<Integer, Integer>, Integer> bp_to_loop2 = new HashMap<>();
		Lib.get_bp_close_loop_and_shrink(DSSR1_parser, bp_to_loop1, npk_bp1, ResID_to_star3d_index1, seed_loop_nt1);
		Lib.get_bp_close_loop_and_shrink(DSSR2_parser, bp_to_loop2, npk_bp2, ResID_to_star3d_index2, seed_loop_nt2);

		for (int i = 0; i < SM_top.size(); i++) {
			Lib.assign_sm_to_loop(SM_top, loopInRNAs_to_sm, bp_to_loop1, bp_to_loop2, i);
		}

		Set<Integer> junction_sm = new HashSet<>();
		Set<Integer> junction_id_1 = new HashSet<>();
		Set<Integer> junction_id_2 = new HashSet<>();
		Lib.get_junction_id(DSSR1_parser, junction_id_1);
		Lib.get_junction_id(DSSR2_parser, junction_id_2);
		for (Pair<Integer, Integer> p : loopInRNAs_to_sm.keySet()) {
			if(junction_id_1.contains(p.left) && junction_id_2.contains(p.right))
				junction_sm.addAll(loopInRNAs_to_sm.get(p));
		}

		Set<Set<Integer>> compatible_loop_sm = new HashSet<>();
		Map<Set<Integer>, Pair<Integer, Integer>> loop_sm_to_loop = new HashMap<>();

		/*********************** for each loop pair, find compatible sm *********************************/
		for (Pair<Integer, Integer> p : loopInRNAs_to_sm.keySet()) {

			Set<Integer> set_of_sm = loopInRNAs_to_sm.get(p);

			List<Integer> list_of_sm_index = new ArrayList<Integer>(set_of_sm);

			List<Pair<Integer, Integer>> SMC_loop = new ArrayList<Pair<Integer, Integer>>();
			exec = Executors.newFixedThreadPool(thread_num);
			List<Future<List<Pair<Integer, Integer>>>> SMC_multi_thread_ret = new ArrayList<Future<List<Pair<Integer, Integer>>>>();

			for (int i = 0; i < thread_num; i++)
				SMC_multi_thread_ret.add(exec.submit(new FindCompatibleSM(i, thread_num, list_of_sm_index, coord1, coord2)));

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

			/*********************** get stack sets *********************************/
			//build smc graph
			Map<Integer, Set<Integer>> SMC_graph = new HashMap<>();
			for (Pair<Integer, Integer> sm_pair : SMC_loop) {
				// if v = 0, rotate pattern unmatch
				if(!SMC_graph.keySet().contains(sm_pair.left))
					SMC_graph.put(sm_pair.left, new HashSet<>());
				if(!SMC_graph.keySet().contains(sm_pair.right))
					SMC_graph.put(sm_pair.right, new HashSet<>());

				if(sm_pair.left!=sm_pair.right) {
					SMC_graph.get(sm_pair.left).add(sm_pair.right);
					SMC_graph.get(sm_pair.right).add(sm_pair.left);
				}
			}

			List<Set<Integer>> SM_clique_loop = new ArrayList<Set<Integer>>();
			Lib.Bron_Kerbosch(SMC_graph, new HashSet<Integer>(), new HashSet<>(SMC_graph.keySet()), new HashSet<Integer>(), SM_clique_loop);

			/*********************** keep loop with rotated match sm *********************************/
			for (Set<Integer> C : SM_clique_loop) {
				if (C.size() == 0)
					continue;

				int rev_sm = 0;
				for (Integer i : C) {
					Stackmap cur_sm = SM_top.get(i);
					if ((cur_sm.i1 - cur_sm.j1) * (cur_sm.i2 - cur_sm.j2) < 0)
						rev_sm += 1;
				}
				if (rev_sm == 0)
					continue;

				 compatible_loop_sm.add(C);
				 loop_sm_to_loop.put(C, p);
			}
		}

		if(compatible_loop_sm.size() == 0){
			System.out.println("Can not find any rotated matched loop.");
			System.exit(0);
		}
		System.out.println("Compatible loop number: " + compatible_loop_sm.size());

		/*********************** only keep loops with maximal nts *********************************/
		List<Integer> bp_in_loop_sm = new ArrayList<>();
		for(Set<Integer> C : compatible_loop_sm){
			int BP_size = 0;
			for (Integer I : C)
				BP_size += SM_top.get(I).size;
			bp_in_loop_sm.add(BP_size);
		}
		Collections.sort(bp_in_loop_sm, Collections.reverseOrder());
		Integer bp_cutoff = bp_in_loop_sm.get(Math.min(seed_num, bp_in_loop_sm.size())-1);

		Set<Set<Integer>> seed_to_rm = new HashSet<>();

		for(Set<Integer> C : compatible_loop_sm){
			int BP_size = 0;
			for (Integer I : C)
				BP_size += SM_top.get(I).size;
			if (BP_size < bp_cutoff)
				seed_to_rm.add(C);
		}
		compatible_loop_sm.removeAll(seed_to_rm);

		Set<Set<Integer>> seed_loop_sm = new HashSet<Set<Integer>>();
		Map<Set<Integer>, List<List<Integer>>> seed_loop_to_info = new HashMap<>();
		// from clique to junction loop alignment
		Map<Pair<Pair<Integer, Integer>, Integer>, Pair<Integer, Integer>> stack_strand_to_loop1 = new HashMap<>();
		Map<Pair<Pair<Integer, Integer>, Integer>, Pair<Integer, Integer>> stack_strand_to_loop2 = new HashMap<>();

		System.out.println("Phase I: seed search");

		/***************************** keep stack sets under 4A *********************************/
		for (Set<Integer> C : compatible_loop_sm) {
			// get stack strand
			List<Pair<Integer, Integer>> Map1_stack_strand = new ArrayList<Pair<Integer, Integer>>();
			List<Pair<Integer, Integer>> Map2_stack_strand = new ArrayList<Pair<Integer, Integer>>();
			Map<Pair<Integer, Integer>, Pair<Integer, Integer>> stk1_to_stk2 = new HashMap<Pair<Integer, Integer>, Pair<Integer, Integer>>();

			Lib.get_stack_stand(C, Map1_stack_strand, Map2_stack_strand, stk1_to_stk2);

			// get stack nt
			Map<Integer, Integer> Map1_to_Map2_index = new HashMap<Integer, Integer>();

			List<Integer> Map1_stack_index = new ArrayList<Integer>();
			List<Integer> Map2_stack_index = new ArrayList<Integer>();

			for (int i = 0; i < Map1_stack_strand.size(); i++) {
				for (Integer i1 = Map1_stack_strand.get(i).left; i1 <= Map1_stack_strand.get(i).right; i1++)
					Map1_stack_index.add(i1);
				for (Integer i2 = Map2_stack_strand.get(i).left; i2 <= Map2_stack_strand.get(i).right; i2++)
					Map2_stack_index.add(i2);
			}

			for (int i = 0; i < Map1_stack_index.size(); i++)
				Map1_to_Map2_index.put(Map1_stack_index.get(i), Map2_stack_index.get(i));

			//get the transition and rotate matrix for the stacking mapping
			Point stack_XC = new Point(-1, -1, -1);
			Point stack_YC = new Point(-1, -1, -1);
			DenseMatrix64F stack_R = new DenseMatrix64F(3,3);

			Boolean check_rmsd = Lib.get_rotation_matrix(Map1_stack_index, Map2_stack_index, stack_XC, stack_YC, stack_R, STAR3D.rmsd_cutoff);

			if(!check_rmsd)	continue;

			/*********************** loop aln in seed *********************************/
			List<Pair<Integer, Integer>> Map1_loop = new ArrayList<Pair<Integer, Integer>>();
			List<Pair<Integer, Integer>> Map2_loop = new ArrayList<Pair<Integer, Integer>>();

			// get junction loop from stack in seeds
			for (int i = 0; i < Map1_stack_strand.size() - 1; i++) {
				if (Map1_stack_strand.get(i).v != Map1_stack_strand.get(i + 1).v || //remove hairpin loop, change to || 8_8
						Map2_stack_strand.get(i).v != Map2_stack_strand.get(i + 1).v) { //remove hairpin loop
                    // adjacent stack, no loop in between
                    if((Map1_stack_strand.get(i).right + 1 == Map1_stack_strand.get(i + 1).left) ||
                            (Map2_stack_strand.get(i).right + 1 == Map2_stack_strand.get(i + 1).left)) //8_7
                        continue;

					Pair<Integer, Integer> cur_loop1 = new Pair<Integer, Integer>(Map1_stack_strand.get(i).right + 1, Map1_stack_strand.get(i + 1).left - 1, 0);
					Pair<Integer, Integer> cur_loop2 = new Pair<Integer, Integer>(Map2_stack_strand.get(i).right + 1, Map2_stack_strand.get(i + 1).left - 1, 0);
					Map1_loop.add(cur_loop1);
					Map2_loop.add(cur_loop2);

					Pair<Pair<Integer, Integer>, Integer> cur_stack_strand_pair1 = new Pair<Pair<Integer, Integer>, Integer>(Map1_stack_strand.get(i), Map1_stack_strand.get(i+1),0);
					Pair<Pair<Integer, Integer>, Integer> cur_stack_strand_pair2 = new Pair<Pair<Integer, Integer>, Integer>(Map2_stack_strand.get(i), Map2_stack_strand.get(i+1),0);
					stack_strand_to_loop1.put(cur_stack_strand_pair1, cur_loop1);
					stack_strand_to_loop2.put(cur_stack_strand_pair2, cur_loop2);
				}
			}

			// dangling end, the loop between the last stack_stand and the first stack_stand
            if((Map1_stack_strand.get(Map1_stack_strand.size()-1).v != Map1_stack_strand.get(0).v) ||
                    (Map2_stack_strand.get(Map2_stack_strand.size()-1).v != Map2_stack_strand.get(0).v)){ //8_7 has a loop
				// loop length must >= 1
				if(Map1_stack_strand.get(Map1_stack_strand.size()-1).right + 1 == Map1_stack_strand.get(0).left ||
						(Map2_stack_strand.get(Map2_stack_strand.size()-1).right + 1 == Map2_stack_strand.get(0).left))
					continue;

			    Pair<Integer, Integer> last_strand1 = Map1_stack_strand.get(Map1_stack_strand.size() - 1);
                Pair<Integer, Integer> first_strand1 = Map1_stack_strand.get(0);
                Pair<Integer, Integer> last_strand2 = Map2_stack_strand.get(Map2_stack_strand.size() - 1);
                Pair<Integer, Integer> first_strand2 = Map2_stack_strand.get(0);

                Pair<Integer, Integer> cur_loop1 = new Pair<Integer, Integer>(last_strand1.right + 1, first_strand1.left - 1, 0);
                Pair<Integer, Integer> cur_loop2 = new Pair<Integer, Integer>(last_strand2.right + 1, first_strand2.left - 1, 0);
                Map1_loop.add(cur_loop1);
                Map2_loop.add(cur_loop2);

                Pair<Pair<Integer, Integer>, Integer> cur_stack_strand_pair1 = new Pair<Pair<Integer, Integer>, Integer>(last_strand1, first_strand1, 0);
                Pair<Pair<Integer, Integer>, Integer> cur_stack_strand_pair2 = new Pair<Pair<Integer, Integer>, Integer>(last_strand2, first_strand2, 0);
                stack_strand_to_loop1.put(cur_stack_strand_pair1, cur_loop1);
                stack_strand_to_loop2.put(cur_stack_strand_pair2, cur_loop2);
            }

			//get the mapping between the loop nt
			List<Pair<Integer, Integer>> Map_loop_align = new ArrayList<Pair<Integer, Integer>>();
			double Map_loop_score = Lib.run_loop_aln(Map_loop_align, Map1_loop, Map2_loop, bp1_map, bp2_map, coord1, coord2, stack_R, stack_XC, stack_YC, false);

			// run loop alignment
			ArrayList<Integer> Map1_index = new ArrayList<Integer>(Map1_stack_index);
			ArrayList<Integer> Map2_index = new ArrayList<Integer>(Map2_stack_index);
			ArrayList<Integer> Map1_loop_index = new ArrayList<Integer>();
			ArrayList<Integer> Map2_loop_index = new ArrayList<Integer>();
			Lib.add_loop_to_map_index(Map_loop_align, Map1_loop_index, Map2_loop_index, Map1_index, Map2_index);

			int aln_stack_nt_number = Map1_stack_index.size();
//			double map_stack_score = match_score * aln_stack_nt_number * stack_discount;
//			double comb_score = Map_loop_score + map_stack_score;

			double Map_rmsd = Geom.cal_rmsd_from_index(Map1_index, Map2_index, coord1, coord2);

			if ((Map_rmsd <= STAR3D.rmsd_cutoff) && (Map1_loop_index.size() >= STAR3D.minimum_aligned_loop_size)) {
				seed_loop_sm.add(C);
				List<List<Integer>> info = new ArrayList<>();

				Integer loop1 = loop_sm_to_loop.get(C).left;
				Integer loop2 = loop_sm_to_loop.get(C).right;

				Set<Integer> loop1_nt = new HashSet<>(DSSR1_parser.loops_all_nt.get(loop1));
				Set<Integer> loop2_nt = new HashSet<>(DSSR2_parser.loops_all_nt.get(loop2));
				Integer len_in_loop = 0;
//				for(int i = 0;i<Map1_index.size();i++){
//					Integer nt1 = ResID1_list.get(Map1_index.get(i)).seqnum;
//					Integer nt2 = ResID2_list.get(Map2_index.get(i)).seqnum;
//					if(loop1_nt.contains(nt1) && loop2_nt.contains(nt2))
//						len_in_loop += 1;
//				}
//				info.add(len_in_loop);
				info.add(new ArrayList<>(loop1_nt));
				info.add(new ArrayList<>(loop2_nt));
				seed_loop_to_info.put(C, info);
			}
		}

		endTime = System.nanoTime();
		totalTime = (endTime - startTime) / 1000000;
		System.out.println(String.format("Seed search time: %d.%d s", totalTime / 1000, totalTime % 1000));

		/********************** Phase II: extension *********************/
		System.out.println("Phase II: extension");

		Map<Integer, Set<Integer>> i1_sm = new HashMap<Integer,  Set<Integer>>();
		Map<Integer, Set<Integer>> j1_sm = new HashMap<Integer,  Set<Integer>>();
		Map<Integer, Set<Integer>> i2_sm = new HashMap<Integer,  Set<Integer>>();
		Map<Integer, Set<Integer>> j2_sm = new HashMap<Integer,  Set<Integer>>();

		// construct hashmaps for extension
		// for all sm
		// i1 -> index_sm
		// j1 -> index_sm
		// i2 -> index_sm
		// j2 -> index_sm
		for(Integer index_sm=0;index_sm<SM_top.size();index_sm++) {
			Stackmap sm = SM_top.get(index_sm);

			if (!i1_sm.keySet().contains(sm.i1))
				i1_sm.put(sm.i1, new HashSet<>());
			if (!j1_sm.keySet().contains(sm.j1))
				j1_sm.put(sm.j1, new HashSet<>());
			if (!i2_sm.keySet().contains(sm.i2))
				i2_sm.put(sm.i2, new HashSet<>());
			if (!j2_sm.keySet().contains(sm.j2))
				j2_sm.put(sm.j2, new HashSet<>());

			i1_sm.get(sm.i1).add(index_sm);
			j1_sm.get(sm.j1).add(index_sm);
			i2_sm.get(sm.i2).add(index_sm);
			j2_sm.get(sm.j2).add(index_sm);
		}

		PriorityQueue<Aln> Alns = new PriorityQueue<Aln>(map_num, Collections.reverseOrder());

		for(Set<Integer> C_seed_loop : seed_loop_sm) {
			// extension sm
			Set<Set<Integer>> extended_sm = new HashSet<>();
			Lib.extend_sm_set(C_seed_loop, i1_sm, j1_sm, i2_sm, j2_sm, new HashSet<>(C_seed_loop), new HashSet<>(C_seed_loop), extended_sm);
			if(extended_sm.isEmpty()) continue;

			//get the stack mapping: Map1_stack, Map2_stack
			for (Set<Integer> s : extended_sm) {
				List<Pair<Integer, Integer>> Map1_stack_strand = new ArrayList<Pair<Integer, Integer>>();
				List<Pair<Integer, Integer>> Map2_stack_strand = new ArrayList<Pair<Integer, Integer>>();
				Map<Pair<Integer, Integer>, Pair<Integer, Integer>> stk1_to_stk2 = new HashMap<Pair<Integer, Integer>, Pair<Integer, Integer>>();
				Lib.get_stack_stand(s, Map1_stack_strand, Map2_stack_strand, stk1_to_stk2);

				// get stack nt index
				List<Integer> Map1_stack_index = new ArrayList<Integer>();
				List<Integer> Map2_stack_index = new ArrayList<Integer>();

				for(int i=0;i<Map1_stack_strand.size();i++) {
					for(Integer i1 = Map1_stack_strand.get(i).left; i1 <= Map1_stack_strand.get(i).right; i1++)
						Map1_stack_index.add(i1);
					for(Integer i2 = Map2_stack_strand.get(i).left; i2 <= Map2_stack_strand.get(i).right; i2++)
						Map2_stack_index.add(i2);
				}

				// add loops to extended stk
				List<Pair<Integer, Integer>> Map1_loop = new ArrayList<Pair<Integer, Integer>>();
				List<Pair<Integer, Integer>> Map2_loop = new ArrayList<Pair<Integer, Integer>>();

				for(int i=0;i<Map1_stack_strand.size()-1;i++){
					Pair<Integer,Integer> left_strand1 = Map1_stack_strand.get(i);
					Pair<Integer,Integer> right_strand1 = Map1_stack_strand.get(i+1);
					Pair<Integer,Integer> left_strand2 = Map2_stack_strand.get(i);
					Pair<Integer,Integer>  right_strand2 = Map2_stack_strand.get(i+1);

					Pair<Integer, Integer> cur_loop1 = new Pair<Integer, Integer>(-1,-1, 0);
					Pair<Integer, Integer> cur_loop2= new Pair<Integer, Integer>(-1,-1, 0);

					if(no_dangle_end == true && ( left_strand1.v == right_strand1.v || left_strand2.v == right_strand2.v
							|| left_strand1.right+1 > right_strand1.left-1 || left_strand2.right+1 > right_strand2.left-1))
						continue;

                    if(stack_strand_to_loop1.keySet().contains(new Pair<Pair<Integer, Integer>, Integer>(left_strand1, right_strand1, 0)) &&
                            stack_strand_to_loop2.keySet().contains(new Pair<Pair<Integer, Integer>, Integer>(left_strand2, right_strand2, 0)) ){ //8_6 seeds internal
						cur_loop1 = stack_strand_to_loop1.get(new Pair<Pair<Integer, Integer>, Integer>(left_strand1, right_strand1, 0));
						cur_loop2 = stack_strand_to_loop2.get(new Pair<Pair<Integer, Integer>, Integer>(left_strand2, right_strand2, 0));
					}else{
						if(!(left_strand1.right == ResID1_list.size()-1 || right_strand1.left== 0 || left_strand1.right+1 == right_strand1.left))
							cur_loop1 = new Pair<Integer, Integer>(left_strand1.right + 1, right_strand1.left - 1, 0);

						if(!(left_strand2.right == ResID2_list.size()-1 || right_strand2.left== 0 || left_strand2.right+1 == right_strand2.left))
							cur_loop2 = new Pair<Integer, Integer>(left_strand2.right + 1, right_strand2.left - 1, 0);
					}

					Map1_loop.add(cur_loop1);
					Map2_loop.add(cur_loop2);
				}

				// add the loop between last and first stackmap
				Pair<Integer,Integer> first_strand1 = Map1_stack_strand.get(0);
				Pair<Integer,Integer> last_strand1 =  Map1_stack_strand.get(Map1_stack_strand.size()-1);
				Pair<Integer,Integer> first_strand2 = Map2_stack_strand.get(0);
				Pair<Integer,Integer> last_strand2 =  Map2_stack_strand.get(Map1_stack_strand.size()-1);

				if(no_dangle_end == false || (( first_strand1.v != last_strand1.v || first_strand2.v != last_strand2.v) &&
						last_strand1.right + 1 <= first_strand1.left - 1 && last_strand2.right + 1 <= first_strand2.left - 1)) {
					Pair<Integer, Integer> cur_loop1 = new Pair<Integer, Integer>(last_strand1.right + 1, first_strand1.left - 1, 0);
					Pair<Integer, Integer> cur_loop2 = new Pair<Integer, Integer>(last_strand2.right + 1, first_strand2.left - 1, 0);
					Map1_loop.add(cur_loop1);//8_4
					Map2_loop.add(cur_loop2);//8_4
				}

				// get rotation matrix from stack
				Point stack_XC = new Point(-1, -1, -1);
				Point stack_YC = new Point(-1, -1, -1);
				DenseMatrix64F stack_R = new DenseMatrix64F(3,3);
				Boolean check_rmsd = Lib.get_rotation_matrix(Map1_stack_index, Map2_stack_index, stack_XC, stack_YC, stack_R, STAR3D.rmsd_cutoff);
				if(!check_rmsd) continue;
				/*************************** loop aln for extended stack sets**************************/
				List<Pair<Integer, Integer>> Map_loop_align = new ArrayList<Pair<Integer, Integer>>();


				double Map_loop_score = Lib.run_loop_aln(Map_loop_align, Map1_loop, Map2_loop, bp1_map, bp2_map, coord1, coord2, stack_R, stack_XC, stack_YC, false);

				ArrayList<Integer> Map1_index = new ArrayList<Integer>(Map1_stack_index);
				ArrayList<Integer> Map2_index = new ArrayList<Integer>(Map2_stack_index);
				ArrayList<Integer> Map1_loop_index = new ArrayList<Integer>();
				ArrayList<Integer> Map2_loop_index = new ArrayList<Integer>();
				Lib.add_loop_to_map_index(Map_loop_align, Map1_loop_index, Map2_loop_index, Map1_index, Map2_index);

				int aln_stack_nt_number = Map1_stack_index.size();

				double map_stack_score = match_score * aln_stack_nt_number * stack_discount;
				double comb_score = Map_loop_score + map_stack_score;

				double Map_rmsd = Geom.cal_rmsd_from_index(Map1_index, Map2_index, coord1, coord2);

				boolean contain_junction = true;
				for(Integer i : C_seed_loop) {
					if (!junction_sm.contains(i)) {
						contain_junction = false;
						break;
					}
				}

				List<List<Integer>> cur_info = seed_loop_to_info.get(C_seed_loop);
				if (Map_rmsd <= STAR3D.rmsd_cutoff){
					Aln new_aln = new Aln(Map1_index, Map2_index, Map1_loop_index, Map2_loop_index,
							Map1_stack_index, Map2_stack_index, Map_rmsd, comb_score, Map_loop_score,
							map_stack_score, s, contain_junction, cur_info, C_seed_loop
							);

					Aln overlapped_aln = Lib.overlapped_with_old_alns(new_aln, Alns);
					if(overlapped_aln==null){
						Alns.add(new_aln);
					}else{
						if(new_aln.compareTo(overlapped_aln) < 0){
							Alns.remove(overlapped_aln);
							Alns.add(new_aln);
						}
					}

					while(Alns.size()>map_num)
						Alns.poll();
				}
			}
		}

		if(STAR3D.pdb_output >0) {
			PriorityQueue<Aln> Alns_copy = new PriorityQueue<Aln>(map_num);
			Alns_copy.addAll(Alns);
			for (int i = 0; i< STAR3D.pdb_output; i++) {
				//get the transition and rotate matrix for the stacking mapping
				Point Map1_XC = new Point(-1, -1, -1);
				Point Map2_YC = new Point(-1, -1, -1);
				DenseMatrix64F Map_R = new DenseMatrix64F(3,3);

				Aln cur_aln = Alns_copy.poll();
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

		Lib.print_result(Alns, PDBID1, chainID1, PDBID2, chainID2);

		System.exit(0);
	}
}
