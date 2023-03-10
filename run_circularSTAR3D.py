import argparse
import os

option_map = {
    "output_prefix": "o",
    "rmsd": "r",
    "minimum_stack_size": "s",
    "gap_open_penalty": "g",
    "gap_extension_penalty": "e",
    "number_of_threads": "t",
    "match_score": "m",
    "mismatch_score": "i",
    "include_dangling_end": "d",
    "stack_distance_cutoff": "c",
    "number_of_alignments": "n",
    "fix_stacks": "f",
    "clique_search_timeout": "l",
    "print_PDB": "p",
    "help": "h"
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-id", help="only for preprocess")
    parser.add_argument("--chain-id", help="only for preprocess")
    parser.add_argument("--pdb-id1", help="only for alignment")
    parser.add_argument("--chain-id1", help="only for alignment")
    parser.add_argument("--pdb-id2", help="only for alignment")
    parser.add_argument("--chain-id2", help="only for alignment")
    parser.add_argument("--non-rotated", action="store_true",
                        help="when using non-rotated, LocalSTAR3D is used for alignment and the rotated match will not be output")
    parser.add_argument("--preprocess", action="store_true",
                        help="run preprocess automatically, including downloading PDB and base-pair annotation")
    parser.add_argument("--rmsd", help="rmsd cutoff for output alignments, default is 4 (Armstrong)")
    parser.add_argument("--minimum-stack-size", help="the minimum conserved stack size, default is 3 (base-pair)")
    parser.add_argument("--gap-open-penalty")
    parser.add_argument("--gap-extension-penalty")
    parser.add_argument("--number-of-threads")
    parser.add_argument("--match-score")
    parser.add_argument("--mismatch-score")
    parser.add_argument("--include-dangling-end",
                        help="whether to include the dangling end in alignment, default is true")
    parser.add_argument("--stack-distance-cutoff",
                        help="the adjacency constrain between conserved stack pairs, default is 15 (nt)")
    parser.add_argument("--number-of-alignments", help="the maximum number of alignment in output, default is 5")
    parser.add_argument("--fix-stacks", help="tolerant 1nt gap in stack, default is true")
    parser.add_argument("--clique-search-timeout")
    parser.add_argument("--print-PDB", help="output PDB for the number of top alignments, default is 0")
    parser.add_argument("--output-prefix", help="the prefix of the output files, default is \"output\"")

    args = parser.parse_args()
    jar = "LocalSTAR3D/LocalSTAR3D.jar" if args.non_rotated else "CircularSTAR3D/CircularSTAR3D.jar"

    if args.preprocess:
        run_preprocess(jar, args)
    else:
        run_align(jar, args)


def run_preprocess(jar, args):
    cmd = f"java -cp {jar} Preprocess {args.pdb_id} {args.chain_id}"
    run_cmd(cmd)


def run_align(jar, args):
    other_args = dict(args.__dict__)
    for a in ["non_rotated", "preprocess", "chain_id", "pdb_id", "chain_id1", "pdb_id1", "chain_id2", "pdb_id2"]:
        other_args.pop(a)

    options = []

    if "output_prefix" not in other_args or other_args["output_prefix"] is None:
        other_args["output_prefix"] = "output"

    for arg_name, arg_val in other_args.items():
        if arg_val:
            options.append("-" + option_map[arg_name])
            options.append(arg_val)

    option_str = " ".join(options)

    cmd = f"java -jar {jar} {option_str} {args.pdb_id1} {args.chain_id1} {args.pdb_id2} {args.chain_id2}"
    run_cmd(cmd)


def run_cmd(cmd):
    os.system(cmd)


if __name__ == "__main__":
    main()
