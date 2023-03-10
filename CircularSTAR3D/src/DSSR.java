import java.util.*;
import java.io.*;

public class DSSR {
    HashSet<Pair<ResID, String>> basepair;
    //List<List<Integer>> junctions;
    //List<List<Integer>> internal_loops;
    List<List<Integer>> loops_close_bp;//including junctions and internal loops
    List<List<Integer>> loops_all_nt;
    HashSet<String> modBase;
    HashMap<String, String> modBaseMap;

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

        String[] token = DSSR_BP_line.trim().split("\\s+");

        ResID R1 = parse_dssr_index(token[1]);
        ResID R2 = parse_dssr_index(token[2]);

        String name = token[4];
        String Saenger = token[5];
        String lw = token[6];

        if (lw.equals("--"))
            return new Pair<ResID, String>(R1, R2, "--h");

        char e1, e2, o;
        e1 = lw.charAt(1);
        e2 = lw.charAt(2);
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

    DSSR(File fname, String chainID) throws IOException {
        basepair = new HashSet<Pair<ResID, String>>();
        modBase = new HashSet<String>();
        modBaseMap = new HashMap<String, String>();
        loops_all_nt = new ArrayList<List<Integer>>();
        loops_close_bp = new ArrayList<List<Integer>>();

        BufferedReader in = new BufferedReader(new FileReader(fname));
        String line;

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

            if(line.trim().matches("List of [\\d]+ base pair.*")) {
                in.readLine(); line = in.readLine();
                while(!line.isEmpty()) {
                    Pair<ResID, String> BP = parse_DSSR_BP_info(line);
                    if (BP != null) basepair.add(BP);
                    line = in.readLine();
                }
            }

            if(line.trim().matches("List of [\\d]+ internal loop.*")) {
                line = in.readLine();
                while(!line.isEmpty()) {
                    String nums = line.split("; ")[1];

                    List<String> strlist = Arrays.asList(nums.substring(1,nums.length()-1).split(","));
                    List<Integer> intlist = new ArrayList<Integer>();
                    for(String s : strlist){
                        intlist.add(Integer.valueOf(s));
                    }

                    line = in.readLine().trim();

                    if(line.startsWith("summary:"))
                        line = in.readLine();

                    boolean wrong_chain = false;//7_28

                    //all nts in the internal loop
                    List<Integer> all_nts =  new ArrayList<Integer>();
                    for(String all_nt_str : Arrays.asList(line.split("nts")[1].split(" ")[2].split(","))) {
                        String cur_chain = Arrays.asList(all_nt_str.split("\\.")).get(0);
                        if (cur_chain.contains(":"))
                            cur_chain = Arrays.asList(cur_chain.split(":")).get(1);
                        if(!cur_chain.equals(chainID)) {//7_28
                            wrong_chain = true;//7_28
                        }
                        List<String> tmp = Arrays.asList(all_nt_str.replaceAll("[^0-9]"," ").split(" "));
                        all_nts.add(Integer.valueOf(tmp.get(tmp.size()-1)));
                    }

                    if(!wrong_chain)//7_29
                        loops_all_nt.add(new ArrayList<>(all_nts));

                    //remove all nts in the loop region, and get the basepairs close the loop
                    line = in.readLine();
                    while(!line.isEmpty() && line.replaceAll(" ", "").startsWith("nts") ) {
                        List<String> tmp = Arrays.asList(line.split(" "));
                        List<String> loop_strs = Arrays.asList(tmp.get(tmp.size()-1).split(","));
                        for(String s : loop_strs){
                            List<String> loop_nt_str_list = Arrays.asList(s.replaceAll("[^0-9]"," ").split(" "));

                            if(loop_nt_str_list.size()>=2) // is something like 0.G1532
                                all_nts.remove(Integer.valueOf(loop_nt_str_list.get(loop_nt_str_list.size()-1)));
                        }
                        line = in.readLine();
                    }
                    if(!wrong_chain)//7_28
                        loops_close_bp.add(new ArrayList<>(all_nts));
                }
            }

            if(line.trim().matches("List of [\\d]+ junction.*")) {
                line = in.readLine();
                while(!line.isEmpty()) {
                    String nums = line.split("; ")[1];

                    List<String> strlist = Arrays.asList(nums.substring(1,nums.length()-1).split(","));
                    List<Integer> intlist = new ArrayList<Integer>();
                    for(String s : strlist){
                        intlist.add(Integer.valueOf(s));
                    }

                    line = in.readLine().trim();
                    if(line.startsWith("summary:"))
                        line = in.readLine();

                    boolean wrong_chain = false;//7_28

                    //all nts in the junctions
                    List<Integer> all_nts =  new ArrayList<Integer>();
                    for(String all_nt_str : Arrays.asList(line.split("nts")[1].split(" ")[2].split(","))) {
                        String cur_chain = Arrays.asList(all_nt_str.split("\\.")).get(0);
                        if (cur_chain.contains(":"))
                            cur_chain = Arrays.asList(cur_chain.split(":")).get(1);
                        if(!cur_chain.equals(chainID)) {//7_28
                            wrong_chain = true;//7_28
                        }
                        List<String> tmp = Arrays.asList(all_nt_str.replaceAll("[^0-9]"," ").split(" "));//only keep the digit, index of nt
                        all_nts.add(Integer.valueOf(tmp.get(tmp.size()-1)));
                    }

                    if(!wrong_chain)//7_29
                        loops_all_nt.add(new ArrayList<>(all_nts));

                    //remove all nts in the loop region, and get the basepairs close the loop
                    line = in.readLine();
                    while(!line.isEmpty() && line.replaceAll(" ", "").startsWith("nts") ) {
                        List<String> tmp = Arrays.asList(line.split(" "));
                        List<String> loop_strs = Arrays.asList(tmp.get(tmp.size()-1).split(","));
                        for(String s : loop_strs){
                            List<String> loop_nt_str_list = Arrays.asList(s.replaceAll("[^0-9]"," ").split(" "));

                            if(loop_nt_str_list.size()>=2) // is something like 0.G1532
                                all_nts.remove(Integer.valueOf(loop_nt_str_list.get(loop_nt_str_list.size()-1)));
                        }
                        line = in.readLine();
                    }
                    if(!wrong_chain)//7_28
                        loops_close_bp.add(new ArrayList<>(all_nts));
                }
                break;
            }
        }
        in.close();
    }
}