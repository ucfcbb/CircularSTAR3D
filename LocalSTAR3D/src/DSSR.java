import java.util.*;
import java.util.regex.*;
import java.io.*;

public class DSSR {
    HashSet<Pair<ResID, String>> basepair;

    HashSet<String> modBase;
    HashMap<String, String> modBaseMap;

    String caret = "^";

    private void parseModifiedBase(String DSSR_MN_line) {
        String[] token = DSSR_MN_line.trim().split("\\s+");

        String[] list = token[3].split(",");
        String nt = token[1].split("-")[0];

        for(String base : list) {
            modBase.add(base);
            modBaseMap.put(base, nt);
        }
    }

    ResID parse_dssr_index (String DSSR_res_index){

        int p = DSSR_res_index.indexOf('.'),
                p1 = DSSR_res_index.indexOf(':'),
                c = DSSR_res_index.indexOf('^'),
                len = DSSR_res_index.length();

        if(!modBase.contains(DSSR_res_index)) {

            String[] index1 = DSSR_res_index.substring(p+1, ((c == -1) ? len : c))
                    .split("(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)");

            String chain = DSSR_res_index.substring(p1+1, p);
            String seqnum = index1[1];
            char icode = ((c == -1) ? ' ' : DSSR_res_index.charAt(len-1));

            return new ResID(chain, Integer.parseInt(seqnum), icode);
        } else {

            String curr_mod_base = this.modBaseMap.get(DSSR_res_index);
            int modBaseLen = curr_mod_base.length(),
                    modBaseEnd = DSSR_res_index.indexOf(curr_mod_base) + modBaseLen - 1;

            if(DSSR_res_index.charAt(modBaseEnd + 1) == '/') modBaseEnd++;

            String chain = DSSR_res_index.substring(0, p);
            String seqnum = DSSR_res_index.substring(modBaseEnd + 1, ((c == -1) ? len : c));
            char icode = ((c == -1) ? ' ' : DSSR_res_index.charAt(len-1));

            return new ResID(chain, Integer.parseInt(seqnum), icode);
        }
    }

    private Pair<ResID, String> parse_DSSR_BP_info(String DSSR_BP_line) {

//        String[] token = DSSR_BP_line.trim().split("\\s+");
//
//        ResID R1 = parse_dssr_index(token[1]);
//        ResID R2 = parse_dssr_index(token[2]);
//
//        String lw = token[6];
//        if(lw.matches("^[ct][WHS][WHS]$")) {
//            lw = token[6].substring(1) + token[6].charAt(0);
//        } else {
//            lw = "--h";
//        }
//
//        Pair<ResID, String> result = new Pair<ResID, String>(R1, R2, lw);
//
//        return result;

        String[] token = DSSR_BP_line.trim().split("\\s+");

        ResID R1 = parse_dssr_index(token[1]);
        ResID R2 = parse_dssr_index(token[2]);

        String name = token[4];
        String Saenger = token[5];
        String lw = token[6];

        if (lw.equals("--"))
            return new Pair<ResID, String>(R1, R2, "--h");

        char e1, e2, o;
        e1=lw.charAt(1);
        e2=lw.charAt(2);
        o = lw.charAt(0);
        if ((e1=='W' || e1=='S' || e1=='H') && (e2=='W' || e2=='S' || e2=='H')){
            if (!name.equals("--"))
                return new Pair<ResID, String>(R1, R2, ""+e1+e2+o);
            else
                return new Pair<ResID, String>(R1, R2, "--n");
        }
        else{
            return new Pair<ResID, String>(R1, R2, "--h");
        }

    }

    DSSR(File fname) throws IOException {
        basepair = new HashSet<Pair<ResID, String>>();
        modBase = new HashSet<String>();
        modBaseMap = new HashMap<String, String>();

        BufferedReader in = new BufferedReader(new FileReader(fname));
        String line;

        String bp_start = "nt1            nt2           bp  name        Saenger    LW  DSSR";
        String end_section = "****************************************************************************";
        boolean readModBases = false;

        while((line = in.readLine()) != null) {

            if(!readModBases && line.trim().matches("List of [\\d]+ type[s]? of [\\d]+ modified nucleotide[s]?")) {
                in.readLine(); line = in.readLine();
                while(!line.isEmpty()) {
                    parseModifiedBase(line);
                    line = in.readLine();
                }
                readModBases = true;
            }

            if(line.trim().matches("List of [\\d]+ base pairs")) {
                in.readLine(); line = in.readLine();
                while(!line.isEmpty()) {
                    Pair<ResID, String> BP = parse_DSSR_BP_info(line);
                    if(BP != null) basepair.add(BP);
                    line = in.readLine();
                }
                break;
            }
        }
        in.close();
    }


    public static void main(String args[]) {
        File fname = new File("/home/xiaoli/software/ori_star/STAR3D_source/pdbx/4xej.dssr");

        try {
            new DSSR(fname);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}