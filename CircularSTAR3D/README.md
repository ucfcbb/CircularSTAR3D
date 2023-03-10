## CircularSTAR3D: a stack-based RNA 3D structural alignment tool for circular matching

### Terms

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)  

Where appropriate, please cite the following circularSTAR3D paper:
Chen et al. "CircularSTAR3D: a stack-based RNA 3D structural alignment tool for circular matching.

### Requirement
CircularSTAR3D is implemented by using java 1.8 and can be executed in 64-bit 
Linux. Two java packages, "commons-cli-1.2.jar" and "EJML-core-0.26.jar", 
are used in the program to support argument parsing and efficient 
matrix computation. The two jar files have been already included into the package.

 "Apache ant" is required to compile the source code.

Apache ant: http://ant.apache.org/bindownload.cgi 
For Linux system, ant can be installed directly by using command line.
Debian/Ubuntu: "apt-get install ant"
Fedora/CentOS: "yum install ant"

CircularSTAR3D extracts base-pair information from DSSR's annotations. CircularSTAR3D has been tested with DSSR basic version 1.5.3 and 2.3.2. There are two ways for users to use DSSR.
1. If you have DSSR software, copy DSSR into CircularSTAR3D/tools/. After that, you will have CircularSTAR3D/tools/DSSR/x3dna-dssr.
2. If you don't have DSSR software, you can copy DSSR annotation files for your PDBs into CircularSTAR3D/STAR3D_struct_info/.
   For example, if you want to preprocess PDB 6d90, you can copy your 6d90.dssr file into CircularSTAR3D/STAR3D_struct_info/.
   The DSSR annotation files that were used to generate the results in our manuscript are provided along with CircularSTAR3D package.

### Installation
Go to the CircularSTAR3D home directory and execute "ant jar".
(The users can also clean the compiled files by executing "ant clean")

### Preprocessing
CircularSTAR3D downloads PDB files and preprocesses them to retrieve the 
secondary structural information.

Go to the CircularSTAR3D home directory and execute 
"java -cp CircularSTAR3D.jar Preprocess [PDB ID] [Chain ID]".

Notice: Make sure the fold "tools" and file "CircularSTAR3D.jar" in the same directory. If you see "Exception in thread "main" java.lang.UnsupportedClassVersionError", try to delete the files in class/ folder and run "ant jar".

### Structural Alignment
Go to the CircularSTAR3D home directory and execute 
"java -jar CircularSTAR3D.jar [PDB1 ID] [Chain1 ID] [PDB2 ID] [Chain2 ID]"
More options for the program can be seen by executing 
"java -jar CircularSTAR3D.jar -h"

Notice: Make sure the folds "PDB", "STAR3D_struct_info" and the file "CircularSTAR3D.jar" 
in the same directory. 

### Example
```
java -cp CircularSTAR3D.jar Preprocess 3f2x X
java -cp CircularSTAR3D.jar Preprocess 4frn B
java -jar CircularSTAR3D.jar -s 2 4frn B 3f2x X
```

### ACKNOWLEDGEMENTS
CircularSTAR3D is developed for an NIH funded project (R01GM102515).

CircularSTAR3D uses java packages Commons CLI and EJML.
Commons CLI: http://commons.apache.org/proper/commons-cli/
EJML: http://ejml.org/wiki/index.php?title=Main_Page
  
### CONTACT
Xiaoli Chen is the author of CircularSTAR3D. For bug reports or comments please contact xiaoli.chen@knights.ucf.edu.
