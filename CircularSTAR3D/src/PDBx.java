import java.io.*;
import java.util.*;

public class PDBx extends PDBParser{
    public PDBx(File fname, String _chainID)  throws IOException{

        chain_res = new HashMap<String, ArrayList<Residue>>();
        chain_atom = new HashMap<String, ArrayList<Atom>>();

        BufferedReader in = new BufferedReader(new FileReader(fname));

        String chain = null;
        char icode = '\0';
        int seqnum = 0, sn = 0;
        String symbol = null;
        double x = 0, y = 0, z = 0;
        float occupancy=1;
        float B_iso_or_equiv=99;
        String atom = null;
        Residue rid;

        String line;
        int symbol_index = -1, chain_index = -1, label_seqnum_index = -1, auth_seqnum_index = -1, atom_index = -1, sn_index = -1,
                x_index = -1, y_index = -1, z_index = -1, icode_index = -1, occupancy_index = -1, B_iso_or_equiv_index = -1;
        boolean parsedAtomItems = false;

        while((line = in.readLine()) != null) {
            if(line.startsWith("_atom_site.")) {
                ArrayList<String> items = new ArrayList<String>();
                while(line.startsWith("_atom_site.")) {
                    items.add(line.trim());
                    line = in.readLine();
                }

                symbol_index = items.indexOf("_atom_site.label_comp_id");
                chain_index = items.indexOf("_atom_site.auth_asym_id");
                label_seqnum_index = items.indexOf("_atom_site.label_seq_id");
                auth_seqnum_index = items.indexOf("_atom_site.auth_seq_id");
                atom_index = items.indexOf("_atom_site.label_atom_id");
                sn_index = items.indexOf("_atom_site.id");
                x_index = items.indexOf("_atom_site.Cartn_x");
                y_index = items.indexOf("_atom_site.Cartn_y");
                z_index = items.indexOf("_atom_site.Cartn_z");
                icode_index = items.indexOf("_atom_site.pdbx_PDB_ins_code");
                occupancy_index = items.indexOf("_atom_site.occupancy");
                B_iso_or_equiv_index = items.indexOf("_atom_site.B_iso_or_equiv");

                if(symbol_index == -1 || chain_index == -1 || label_seqnum_index == -1 || auth_seqnum_index == -1 || atom_index == -1 || sn_index == -1 ||
                        x_index == -1 || y_index == -1 || z_index == -1 || icode_index == -1 || occupancy_index == -1 || B_iso_or_equiv_index == -1) {
                    System.err.println("The mmCIF structure does not have all the necessary _atom_site attributes.");
                    return;
                }

                parsedAtomItems = true;
            }

            if (parsedAtomItems) {
                line = line.trim();
                if(line.equals("#")) break;

                String token[] = line.split("\\s+");

                if(chain_index != -1) chain = token[chain_index];
                if(!chain.equals(_chainID)) continue;

                if(symbol_index != -1) symbol = token[symbol_index];


                if(label_seqnum_index != -1 && auth_seqnum_index != -1 && token[label_seqnum_index].matches("^([+-]?(\\d+\\.)?\\d+)$")) {
                    seqnum = Integer.valueOf(token[auth_seqnum_index]);
                } else {
                    continue;
                }

                if(atom_index != -1) atom = token[atom_index].replace("\"", "");
                if(sn_index != -1) sn = Integer.valueOf(token[sn_index]);
                if(x_index != -1) x = Double.valueOf(token[x_index]);
                if(y_index != -1) y = Double.valueOf(token[y_index]);
                if(z_index != -1) z = Double.valueOf(token[z_index]);

                if(occupancy_index != -1) occupancy = Float.valueOf(token[occupancy_index]);
                if(B_iso_or_equiv_index != -1) B_iso_or_equiv = Float.valueOf(token[B_iso_or_equiv_index]);

                if(icode_index != -1) {
                    if(!token[icode_index].equals("?"))
                        icode = token[icode_index].charAt(0);
                    else
                        icode = ' ';
                } else {
                    icode = '\0';
                }

                if(!chain_res.containsKey(chain)) {
                    chain_res.put(chain, new ArrayList<Residue>());
                }

                if(!chain_atom.containsKey(chain)) {
                    chain_atom.put(chain, new ArrayList<Atom>());
                }

                rid = new Residue(new ResID(chain, seqnum, icode), symbol);

                if(!chain_res.get(chain).contains(rid)) {
                    chain_res.get(chain).add(rid);
                }

                chain_atom.get(chain).add(new Atom(sn, rid, atom, new Point(x, y, z), occupancy, B_iso_or_equiv));
            }
        }
        in.close();
    }
}
