import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;


public class PDBParser {
    ArrayList<String> Reijmers_atoms=new ArrayList<String>(Arrays.asList("C3'", "C4'", "C5'", "O3'", "O5'", "P"));	//heavy atoms used in Reijmers' paper

    public HashMap<String, ArrayList<Residue>> chain_res = new HashMap<String, ArrayList<Residue>>();;	//key: chain; value: residue
    public HashMap<String, ArrayList<Atom>> chain_atom=  new HashMap<String, ArrayList<Atom>>();	//key: chain; value: coordinates of residue;

    // Return mapping between chain and their residues.
    public HashMap<String, ArrayList<Residue>> get_chain_res() {
        return chain_res;
    }
    public HashMap<String, ArrayList<Atom>> get_chain_atom() {
        return chain_atom;
    }
    //return the RNA sequence for chain (N for unknown residue)

    //return the RNA sequence for chain (N for unknown residue)
    public String get_chain_seq(String chainID){
        String seq="";

        for(Residue R: chain_res.get(chainID)){
            if(R.symbol.equals("A") || R.symbol.equals("C")|| R.symbol.equals("G") || R.symbol.equals("U"))
                seq+=R.symbol;
            else
                seq+="N";
        }
        return seq;
    }

    public ArrayList<Point> get_chain_centroid(String chainID){
        ArrayList<ArrayList<Point>> coord = new ArrayList<ArrayList<Point>>();
        ArrayList<Point> centroid = new ArrayList<Point>();

        ResID cur_rid=new ResID(null, -1, '\0');

        for(Atom A: chain_atom.get(chainID)){
            if (A.res.rid.equals(cur_rid)==false) {coord.add(new ArrayList<Point>());  cur_rid=A.res.rid;};
            if (Reijmers_atoms.contains(A.atom)) coord.get(coord.size()-1).add(A.coord);
        }

        for(ArrayList<Point> L: coord){
            centroid.add(Geom.centroid(L));
        }

        return centroid;
    }
}