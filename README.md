LSP
===
LSP is a project for Computer Science 75/175 at Dartmouth College which implements an algorithm for loop structure 
prediction. 

Our code and supporting materials are contained in the project LSP (short for loop structure prediction). We 
acquired our loop data from Dr. Yoonjoo Choi; the provided dataset included protein sequence and structure 
information for over 32,000 proteins in PDB file format. For our prediction algorithm, we first extract loop 
sequence and structure data from the PDB files, along with information about each loop's anchor points 
and flanking SSEs. We use this loop data to cluster the loops hierarchically into groups of structurally 
similar loops, and then generate a representative probabilistic sequence for each of these loop families. 
Once our source loop data has been clustered, we have the infrastructure necessary to classify new loops. 
Given a new loop sequence, we first match it to the available loop clusters based on sequence similarity 
to each cluster's probabilistic sequence. To predict the loop's structure, we then transfer the structural 
characteristics of the most closely matched loop cluster to the new loop.

Usage:
===
In order to run this code you will need to make sure that:
- (!) you have functioning DSSP executable in the root of the project
- you have a pdb directory in the root of the project with PDB files
- by executing main.py, the cross-validation will kick-off; we have yet to
implement the user-interface for structure prediction as we've been focusing
on testing our methods. 

Data:
===
- PDB files complements of Yoonjoo Choi. [1]

Dependencies:
===
- Python 2.7
- Biopython, see https://github.com/biopython/biopython
- DSSP, see http://swift.cmbi.ru.nl/gv/dssp/ (If you are running Mac OSX and are having trouble getting DSSP installed, I suggest this resource: http://proteinz.blogspot.com/2013/02/compiling-dssp-on-osx-lion-redux.html).

References:
===
[1] Choi, Y. and C.M. Deane. (2009) “FREAD revisited: Accurate loop structure prediction using a 
database search algorithm.” Proteins: Structure, Function, and Bioinformatics, 78 (6).

[2] Verschueren, E. et al. (2011) “Protein design with fragment databases.” Current Opinion in 
Structural Biology, 21 (4). 

License:
===
Licensed under the Apache License, Version 2.0 (the "License"); you may not use this except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 either express or implied. See the License for the specific language governing permissions and limitations under the License.

Authors:
===
Ryan Amos <<Ryan.B.Amos.16@dartmouth.edu>>

Travis Peters <<Travis.W.Peters.GR@dartmouth.edu>>

Serena Liu <<serena.x.liu.14@dartmouth.edu>>
