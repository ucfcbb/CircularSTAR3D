import java.util.*;
import java.util.regex.*;
import java.io.*;

public class MCA {
	HashSet<Pair<ResID, String>> basepair;

	//return the end point of first residue in bp section of MC-Annotate
	private int parse_BP(String BP){ // A1-A93
		String pattern="(([a-zA-Z])|('\\d'))(-)*(\\d)+(\\.[a-zA-Z])*";
		Pattern r=Pattern.compile(pattern);
		Matcher m=r.matcher(BP);

		if (m.find()){
			// System.out.println("BP=" + m.group());
			return m.end(0);
		}
		else{
			System.out.println("Can not find the residue of base pair in MC-Annotate");
			return -1;
		}
	}

	private Pair<ResID, String> parse_MCA_BP_info(String MCA_BP_line){

		String[] token=MCA_BP_line.split(" : "); // A1-A93 : G-C Ww/Ww pairing antiparallel cis XIX
		String BP=token[0]; // A1-A93
		String BP_info=token[1]; // G-C Ww/Ww pairing antiparallel cis XIX

		String index1, index2;

		int r1_e=parse_BP(BP); // A1-A93 => 2
		index1=BP.substring(0, r1_e); // A1-A93 => A1

		int r2_e=parse_BP(BP.substring(r1_e+1)); // A93 => 2
		index2=BP.substring(r1_e+1, r1_e+r2_e+1); // A93
		//System.out.println(index1 + " " + index2);

		// ResID R1=new ResID(index1);
		// ResID R2=new ResID(index2);

		ResID R1 = parse_mca_index(index1); // A1 => new ResID(chainID (A), seqnum (1), icode (' '))
		ResID R2 = parse_mca_index(index2); // A93 => new ResID(chainID (A), seqnum (93), icode (' '))

		char e1, e2, o;
		if(BP_info.contains("cis")){
			o='c';
		}
		else if(BP_info.contains("trans")){
			o='t';
		}
		else{
			o='n';
		}

		token=BP_info.split("\\s+");
		String edge=token[1];

		/*for(int i = 0; i < token.length; i++) {
			System.out.print(token[i] + ", ");
		}
		System.out.print("\n");*/

		if(edge.contains("/")==false) return new Pair<ResID, String>(R1, R2, "--h");

		token=edge.split("/");

		e1=token[0].charAt(0);
		e2=token[1].charAt(0);
		if ((e1=='W' || e1=='S' || e1=='H') && (e2=='W' || e2=='S' || e2=='H') && o!='n'){
			return new Pair<ResID, String>(R1, R2, ""+e1+e2+o);
		}
		else{
			return new Pair<ResID, String>(R1, R2, "--h");
		}

	}

	//retrieve all the base pairs from the output file of MCAnnotate
	MCA(File fname) throws IOException{
		basepair=new HashSet<Pair<ResID, String>>();

		boolean in_bp=false;

		BufferedReader in=new BufferedReader(new FileReader(fname));
		String line;
		while((line=in.readLine())!=null){
			if(line.startsWith("Residue conformations")) {if(in_bp==true) break;}	//only the first base pair section is used (from the 1st MODEL)
			else if(line.startsWith("Adjacent stackings")) in_bp=false;
			else if(line.startsWith("Non-Adjacent stackings")) in_bp=false;
			else if(line.startsWith("Base-pairs")) in_bp=true;
			else if(in_bp==true){
				line=line.trim();
				Pair<ResID, String> BP=parse_MCA_BP_info(line);
				//System.out.println(BP.left + " " + BP.right + " " + BP.v);
				basepair.add(BP);
			}
		}
		in.close();
	}

	ResID parse_mca_index (String MCA_res_index){
		String[] token;
		String chainID_seqnum;

		String chainID;
		int seqnum;
		char icode;

		if(MCA_res_index.contains(".")){
			token=MCA_res_index.split("\\.");
			chainID_seqnum=token[0];
			icode=token[1].charAt(0);
		}
		else{
			chainID_seqnum=MCA_res_index;
			icode=' ';
		}

		if(chainID_seqnum.contains("\'")){
			chainID = "" + chainID_seqnum.charAt(1);
			seqnum = Integer.parseInt(chainID_seqnum.substring(3));
		}
		else{
			chainID = "" + chainID_seqnum.charAt(0);
			seqnum = Integer.parseInt(chainID_seqnum.substring(1));
		}

		return new ResID(chainID, seqnum, icode);
	}
}