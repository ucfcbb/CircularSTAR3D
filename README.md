## CircularSTAR3D: a stack-based RNA 3D structural alignment tool for circular matching

### Terms

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)  

Where appropriate, please cite the following circularSTAR3D paper:
Chen et al. "CircularSTAR3D: a stack-based RNA 3D structural alignment tool for circular matching.

### Requirements
CircularSTAR3D is implemented by using java 1.8 and can be executed in 64-bit 
Linux. Two java packages, "commons-cli-1.2.jar" and "EJML-core-0.26.jar", 
are used in the program to support argument parsing and efficient 
matrix computation. The two jar files have been already included into the package.

 "Apache ant" is required to compile the source code.

Apache ant: http://ant.apache.org/bindownload.cgi 
For Linux system, ant can be installed directly by using command line.  
Debian/Ubuntu: "apt-get install ant"  
Fedora/CentOS: "yum install ant"   

### Installation
CircularSTAR3D package provides compiled jar files that you can run directly.  
If you want to build from the source code, go to the root directory of this package and run "chmod +x ./setup.sh && ./setup.sh". 

### Preprocessing
CircularSTAR3D downloads PDB files and preprocesses them to retrieve the 
secondary structural information.

Preprocess for circular matching, go to the root directory of CircularSTAR3D and execute 
"python3 run_circularSTAR3D.py --preprocess --pdb-id [PDB ID] --chain-id [CHAIN ID]".

Preprocess for regular local alignment, go to the root directory of CircularSTAR3D and execute 
"python3 run_circularSTAR3D.py --non-rotate --preprocess --pdb-id [PDB ID] --chain-id [CHAIN ID]".

If you want to use your custom PDB files, copy them into PDB/ before running the preprocess command. CircularSTAR3D will skip downloading from PDB website if the PDB files present in the PDB folder. Please use suffix ".pdb" for PDB format and suffix ".cif" for PDBx/mmcif format files.    

CircularSTAR3D extracts base-pair information from DSSR's annotations. CircularSTAR3D has been tested with DSSR basic version 1.5.3 and the version in http://skmatic.x3dna.org/ as of 01/31/2023. The DSSR annotation files that were used to generate the results in our manuscript are provided along with CircularSTAR3D package. There are two ways for users to use DSSR for their PDB files.
1. If you have DSSR software, copy DSSR into CircularSTAR3D/tools/. After that, you will have CircularSTAR3D/tools/DSSR/x3dna-dssr.
2. If you don't have DSSR software, you can copy DSSR annotation files for your PDBs into CircularSTAR3D/STAR3D_struct_info/.
   For example, if you want to preprocess PDB 6d90, you can copy your 6d90.dssr file into CircularSTAR3D/STAR3D_struct_info/. The DSSR annotation files can be downloaded from http://skmatic.x3dna.org/, where you can specifiy the PDB ID if it is in PDB database or upload the PDB file to run DSSR webserver. Then copy the content in "DSSR-derived features in text format" into a file in CircularSTAR3D/STAR3D_struct_info with file name "PDB_ID.dssr", replace PDB_ID with the PDB ID that you use to run CircularSTAR3D.

Notice: Make sure the fold "tools" and file "CircularSTAR3D.jar" in the same directory. If you see "Exception in thread "main" java.lang.UnsupportedClassVersionError", try to delete the files in class/ folder and run "ant jar".

### Structural Alignment
Go to the CircularSTAR3D home directory and execute 
"python3 run_circularSTAR3D.py [options]"  
More options for the program can be seen by executing 
"python3 run_circularSTAR3D.py -h"

Use "--non-rotate" to generate local alignment from LocalSTAR3D.  
Use "--output-prefix" to set the prefix of the output files, The default prefix is "output".  
Use "--print-PDB n" to generate PDB for top n alignments. The output PDB files will be named by using the same prefix of the output alingment file. For example, if the prefix is "4frn_B_3f2x_X", the output alignment file will be 4frn_B_3f2x_X.aln and the PDB file of the first alignment will be 4frn_B_3f2x_X.aln1.pdb. The default prefix is "output". The users can add path string to the prefix. The output files will be in the CircularSTAR3D root directory by default.  

Notice: Make sure the folds "PDB", "STAR3D_struct_info" and the file "CircularSTAR3D.jar" 
in the same directory. 

### Output files  
The alignment output files have the suffix "aln". Each alignment output file includes the key parameters that are used to run CircularSTAR3D and the sorted local alignments. Each local alignment includes the rank of the alignment, the alignment score, the alignment nucleotide number, RMSD, and the nucleotide mapping. The format of the nucleotide mapping is "RNA1_chain_id:RNA1_nucleotide_id<->RNA2_chain_id:RNA2_nucleotide_id".  

The PDB output files have the suffix "alnx.pdb", where "x" is the rank of the local alignment. The PDB files store the coordinate data for the local alignment that can be used for visualization and futher analysis. These files are in PDB format ([an introduction of PDB format](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html)). Each PDB ouput file includes 2 models. Model 1 and Model 2 are the aligned substructures in input RNA 1 and RNA 2 respectively. The coodinates of the aligned substructure in input RNA 1 are rotated and translated to superimpose with RNA 2. The coodinates of the aligned substructure in input RNA 2 are the same with the original input files. To view the alignment in PyMOL, click "Movie" in the menu and select "Show All States".  

### Examples
```
python3 run_circularSTAR3D.py --preprocess --pdb-id 3f2x --chain-id X
python3 run_circularSTAR3D.py --preprocess --pdb-id 4frn --chain-id B
python3 run_circularSTAR3D.py --print-PDB 1 --minimum-stack-size 2 --pdb-id1 4frn --chain-id1 B --pdb-id2 3f2x --chain-id2 X

python3 run_circularSTAR3D.py --non-rotate --preprocess --pdb-id 3f2x --chain-id X  
python3 run_circularSTAR3D.py --non-rotate --preprocess --pdb-id 4frn --chain-id B  
python3 run_circularSTAR3D.py --non-rotate --minimum-stack-size 2 --pdb-id1 4frn --chain-id1 B --pdb-id2 3f2x --chain-id2 X
```

### Acknowledgments
CircularSTAR3D is developed for an NIH funded project (R01GM102515).

CircularSTAR3D uses java packages Commons CLI and EJML.  
Commons CLI: http://commons.apache.org/proper/commons-cli/  
EJML: http://ejml.org/wiki/index.php?title=Main_Page  

CircularSTAR3D modified and reused code from [STAR3D](http://genome.ucf.edu/STAR3D/) for stack detection and loop alignment.  
  
### Contact
Xiaoli Chen is the author of CircularSTAR3D. For bug reports or comments please contact xiaoli.chen@knights.ucf.edu.
