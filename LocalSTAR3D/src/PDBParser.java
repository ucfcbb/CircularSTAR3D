import java.util.ArrayList;
import java.util.HashMap;


public interface PDBParser {
    public HashMap<String, ArrayList<Residue>> chain_res = new HashMap<String, ArrayList<Residue>>();;	//key: chain; value: residue
    public HashMap<String, ArrayList<Atom>> chain_atom=  new HashMap<String, ArrayList<Atom>>();	//key: chain; value: coordinates of residue;

    // Return mapping between chain and their residues.
    public HashMap<String, ArrayList<Residue>> get_chain_res();

    // Return mapping between chain and their atoms.
    public HashMap<String, ArrayList<Atom>> get_chain_atom();

    //return the RNA sequence for chain (N for unknown residue)
    public String get_chain_seq(String chainID);

    public ArrayList<Point> get_chain_centroid(String chainID);
}