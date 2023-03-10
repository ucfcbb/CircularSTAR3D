import java.io.*;
import java.util.*;

public class PDB extends PDBParser {
	public PDB(File fname, String _chainID) throws IOException{
		
		chain_res=new HashMap<String, ArrayList<Residue>>();
		chain_atom=new HashMap<String, ArrayList<Atom>>();
		HashSet<String> chain_TER=new HashSet<String>();	//the terminated chains
		
		String line;
		
		String chainID = null;
		char icode='\0';
		int seqnum=0;
		int sn=0;
		String symbol;
		double x=0., y=0., z=0.;
		float occupancy=1;
		float B_iso_or_equiv=99;
		String atom;
		ResID rid;
		ResID cur_rid=new ResID(null, -1, '\0');	//set current resid to null

		BufferedReader in=new BufferedReader(new FileReader(fname));
		
		while((line=in.readLine())!=null){
			if(line.startsWith("ENDMDL")) break;	//only consider the first model
			if(line.startsWith("TER")) { chain_TER.add(chainID); cur_rid=new ResID(null, -1, '\0');}	//add the terminated chain to chain_TER and set the current residue to null
			if(line.startsWith("ATOM") || line.startsWith("HETATM")){
				chainID = "" + line.charAt(21);
				if(chain_TER.contains(chainID)==true) continue;		//if chain is terminated, do nothing
				if(chain_res.get(chainID)==null) chain_res.put(chainID, new ArrayList<Residue>());	//create the residue list and coordinates list for the chain (first appearance)
				if(chain_atom.get(chainID)==null) chain_atom.put(chainID, new ArrayList<Atom>());	
				
				symbol=line.substring(17, 20).trim();
				seqnum=Integer.valueOf(line.substring(22, 26).trim());
				sn=Integer.valueOf(line.substring(6, 11).trim());
				icode=line.charAt(26);
				atom=line.substring(12, 16).trim();
				x=Double.valueOf(line.substring(30, 38).trim());
				y=Double.valueOf(line.substring(38, 46).trim());
				z=Double.valueOf(line.substring(46, 54).trim());

				occupancy=Float.valueOf(line.substring(56, 61).trim());
				B_iso_or_equiv=Float.valueOf(line.substring(61, 70).trim());

				rid=new ResID(chainID, seqnum, icode);
				
				if(rid.equals(cur_rid) == false) {	//new residue
					chain_res.get(chainID).add(new Residue(rid, symbol));
					cur_rid=rid;
				}
				chain_atom.get(chainID).add(new Atom(sn, new Residue(cur_rid, symbol), atom, new Point(x, y, z), occupancy, B_iso_or_equiv));	//add the new atom
			}
		}
		in.close();
	}

	/*
	public ArrayList<Point> get_chain_coord(Character chainID){
		ArrayList<Point> coord =new ArrayList<Point>();
		for(ArrayList<Point> L: chain_coord.get(chainID)){
			Point P=Geom.centroid(L);
			coord.add(P);
		}
		return coord;
	}
	*/
}
