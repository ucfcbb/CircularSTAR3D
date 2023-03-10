import java.io.*;
import java.net.URL;
import java.net.HttpURLConnection;
import java.nio.file.*;
import java.util.*;

public class Preprocess {
	public static void main(String[] args) throws Exception{
		if(args.length!=2){
			System.err.println("usage: java -cp LocalSTAR3D.jar Preprocess PDB chain");
			System.exit(0);
		}

		String DSSR_err_msg = "Failed to run DSSR or locate DSSR annotations.\n" +
				"How to solve this issue?\n" +
				"1. Download DSSR into tools/. After that, you will have tools/DSSR/x3dna-dssr.\n" +
				"2. If you don't have DSSR, you can copy the DSSR annotation for your PDB into STAR3D_struct_info/." +
				"For example, if you want to preprocess PDB 6d90, you can manually copy your 6d90.dssr file into STAR3D_struct_info/" +
				"For users to duplicate the results in our NAR paper, the related DSSR files are provided along with the package.";

		String PDBID=args[0];
		String chainID=args[1];

		int pdb_type=-1;//0 pdb, 1 pdbx

		//get the path STAR3D.jar
		String STAR3D_PATH = new File(Preprocess.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParentFile().getPath();
		String ROOT_PATH = new File(STAR3D_PATH).getParentFile().getPath();

		File PDB_DATA_PATH=new File(ROOT_PATH, "PDB");
		File DSSR_ANNO_PATH=new File(ROOT_PATH, "DSSR_annotation");
		File SI_DATA_PATH=new File(STAR3D_PATH, "STAR3D_struct_info");

		if(!(PDB_DATA_PATH.exists() && PDB_DATA_PATH.isDirectory())) PDB_DATA_PATH.mkdir();
		if(!(DSSR_ANNO_PATH.exists() && DSSR_ANNO_PATH.isDirectory())) DSSR_ANNO_PATH.mkdir();
		if(!(SI_DATA_PATH.exists() && SI_DATA_PATH.isDirectory())) SI_DATA_PATH.mkdir();

		PDBID=PDBID.toLowerCase();
		//chain ID can be lower case
		//chainID=chainID.toUpperCase();

		File PDB_fn=new File(PDB_DATA_PATH, PDBID+".pdb");
		File PDBx_fn=new File(PDB_DATA_PATH, PDBID+".cif");
		URL link = null;
		HttpURLConnection connection = null;

		if((!PDB_fn.exists() || ! PDB_fn.isFile()) && (!PDBx_fn.exists() || ! PDBx_fn.isFile())){

			System.out.println("Downloading PDB file for "+PDBID+"...");

			List<String> pdb_url = new ArrayList<>();
			pdb_url.add("https://files.rcsb.org/download/");
			pdb_url.add("http://www.rcsb.org/pdb/files/");
			pdb_url.add("https://www.rcsb.org/structure/");

			List<String> pdb_suffix = new ArrayList<>();
			pdb_suffix.add(".pdb");
			pdb_suffix.add(".cif");

			int pdb_code = -1;
			for(String url : pdb_url){
				for(String suffix : pdb_suffix){
					link=new URL(url + PDBID + suffix);
					connection = (HttpURLConnection)link.openConnection();
					connection.setRequestMethod("GET");
					connection.connect();

					pdb_code = connection.getResponseCode();

					if(pdb_code<400) {
						if (suffix == ".pdb")
							pdb_type = 0;
						if (suffix == ".cif")
							pdb_type = 1;
						break;
					}
				}
				if (pdb_type != -1)
					break;
			}

			if (pdb_type == -1){
				System.out.println("Cannot download PDB file automatically. You can download manually into the PDB folder.");
				System.exit(0);
			}

			InputStream in =new BufferedInputStream(link.openStream());
			ByteArrayOutputStream out=new ByteArrayOutputStream();
			byte[] buf=new byte[1024];
			int n=0;
			while(-1!=(n=in.read(buf))) out.write(buf, 0, n);
			out.close();
			in.close();
			byte[] response=out.toByteArray();

			FileOutputStream fos = null;
			if(pdb_type==0)
				fos=new FileOutputStream(PDB_fn);
			if(pdb_type==1)
				fos=new FileOutputStream(PDBx_fn);
			fos.write(response);
			fos.close();
		}else{
			if(PDB_fn.exists() && PDB_fn.isFile())
				pdb_type=0;
			if(PDBx_fn.exists() && PDBx_fn.isFile())
				pdb_type=1;
		}

		File anno_file = null;
		anno_file = new File(DSSR_ANNO_PATH, PDBID + ".dssr");
		File DSSR_exe_file = new File("tools/DSSR/x3dna-dssr\"");

		if(anno_file.exists() && anno_file.length() >= 0){
			System.out.println("Using existing DSSR annotation files");
		} else {
			try {
				Process p = null;
				if (pdb_type == 1) {
					System.out.println("Processing " + PDBID + ".cif with DSSR...");
					p = Runtime.getRuntime().exec(new File(STAR3D_PATH, "tools/DSSR/x3dna-dssr").toString() + " --nested -i=" + PDBx_fn);
				}
				if (pdb_type == 0) {
					System.out.println("Processing " + PDBID + ".pdb with DSSR...");
					p = Runtime.getRuntime().exec(new File(STAR3D_PATH, "tools/DSSR/x3dna-dssr").toString() + " --nested -i=" + PDB_fn);
				}

				BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
				PrintWriter out = new PrintWriter(anno_file);
				String line;
				while ((line = in.readLine()) != null) out.println(line);
				out.close();
				in.close();
				p.waitFor();

				p = Runtime.getRuntime().exec(new File(STAR3D_PATH, "tools/DSSR/x3dna-dssr").toString() + " --cleanup");
				p.waitFor();
			} catch(Exception e)  {
				System.out.println(DSSR_err_msg);
				System.exit(1);
			}
		}

		File npk_ct_fn=new File(SI_DATA_PATH, PDBID+'_'+chainID+".ct");
//		if(!npk_ct_fn.exists() || npk_ct_fn.length()==0){
		if(true){
			PDBParser pdb_parser=null;
			if(pdb_type==0){
				pdb_parser=new PDB(PDB_fn,chainID);
			}
			if(pdb_type==1) {
				pdb_parser = new PDBx(PDBx_fn, chainID);
			}

			String chain_seq=pdb_parser.get_chain_seq(chainID);

			//get the base pairs
			HashSet<Pair<Integer, String>> chain_bp;

			HashSet<Pair<ResID, String>> basepair = new HashSet<Pair<ResID, String>>();

			DSSR dssr_parser = new DSSR(anno_file, chainID);

			basepair = dssr_parser.basepair;

			chain_bp=Lib.get_seq_bp(pdb_parser, basepair, chainID);
			Object[][] m=new Object[chain_seq.length()][6];

			for(int i=0; i<chain_seq.length(); i++){
				m[i][0]=new Integer(i+1);
				m[i][1]=new Character(chain_seq.charAt(i));
				m[i][2]=new Integer(i);
				m[i][3]=new Integer(i+2);
				m[i][4]=new Integer(0);
				m[i][5]=new Integer(i+1);
			}

			List<Integer> multi_pair=new ArrayList<Integer>();

			//number of Watson-Crick base pairs
			//xiaoli get wwc first, then get the ones with name and seanger
			int num_WC=0;
			List<Pair<Integer, String>> to_rm = new ArrayList<>();
			Set<Integer> paired_nt = new HashSet<>();

			for(Pair<Integer, String> P : chain_bp){
				if(P.v.equals("WWc")){
					if((Integer)m[P.left][4]!=0) multi_pair.add(P.left);
					m[P.left][4]=P.right+1;
					if((Integer)m[P.right][4]!=0) multi_pair.add(P.right);
					m[P.right][4]=P.left+1;
					num_WC+=1;
					paired_nt.add(P.left);
					paired_nt.add(P.right);
					to_rm.add(P);
				}
			}
			chain_bp.removeAll(to_rm);
			to_rm.clear();
			for(Pair<Integer, String> P : chain_bp){
				if(paired_nt.contains(P.left) || paired_nt.contains(P.right))
					continue;
				if(!P.v.equals("--n") && !P.v.equals("--h")){
					if((Integer)m[P.left][4]!=0)
						multi_pair.add(P.left);
					else
						m[P.left][4]=P.right+1;
					if((Integer)m[P.right][4]!=0)
						multi_pair.add(P.right);
					else
						m[P.right][4]=P.left+1;
					paired_nt.add(P.left);
					paired_nt.add(P.right);
					to_rm.add(P);
				}
			}
			chain_bp.removeAll(to_rm);
			for(Pair<Integer, String> P : chain_bp){
				if(paired_nt.contains(P.left) || paired_nt.contains(P.right))
					continue;
				if(P.v.equals("--n")){
					if((Integer)m[P.left][4]!=0)
						multi_pair.add(P.left);
					else
						m[P.left][4]=P.right+1;
					if((Integer)m[P.right][4]!=0)
						multi_pair.add(P.right);
					else
						m[P.right][4]=P.left+1;
				}
			}

			File ct_fn=new File(SI_DATA_PATH, PDBID+'_'+chainID+".ct");
			PrintWriter out=new PrintWriter(new BufferedWriter(new FileWriter(ct_fn)));
			out.printf("%d\t%s_%s\n", chain_seq.length(), PDBID, chainID);
			for(int i=0; i<chain_seq.length(); i++){
				out.printf("%d\t%s\t%d\t%d\t%d\t%d\n", m[i][0], m[i][1], m[i][2], m[i][3], m[i][4], m[i][5]);
			}
			out.close();

//			ProcessBuilder pb = new ProcessBuilder();
//			Map<String, String> env = pb.environment();
//			env.put("DATAPATH", "tools/RNAstructure/data_tables");
//			p.waitFor();
		}
	}
}
