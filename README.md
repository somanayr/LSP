LSP
===
LSP is a project for Computer Science 75/175 at Dartmouth College which implements an algorithm for loop structure 
prediction. 

Loop structure is highly variable and, thus, is difficult to model accurately. However, loop configuration plays 
an integral role in determining protein tertiary structure and can greatly impact protein function and binding. 
Hence, developing reliable and accurate methods for loop structure prediction is essential in order to better 
model protein structure and function.

Previous research in this area has generally taken one of two approaches: (1) ab initio methods that use a statistical 
model to impute loop structure from the primary amino acid sequence, and (2) database search methods that predict 
a new loop structure based on similarity to previously known loop structures. The effectiveness of database-search 
methods was initially constrained by the limited availability of determined protein structures, but more recent 
implementations of database-search techniques, like FREAD (Choi 2009), have been more successful and have even been known
to outperform ab initio methods for some loop subsets.

TODO
===
We can use this as a main source for keeping track of our own progress as well as the progress of others. We can also use
this space for taking general notes about the project. 

A link to where you can learn about all this stuff can be found [here](https://help.github.com/articles/github-flavored-markdown).

**Serena:**
- [ ] K-means clustering (noise - later)

**Ryan:**
- [ ] Scoring of loops (1) RMSD
  - [ ] align sequences (NW)
  - [ ] Use RMSD to compare aligned sequences
  - [ ] Gap scoring?

**Travis:**
- [ ] Code to extract loop sequences from PDB file(s)
- [ ] Email Yoonjoo

**All:**
- [ ] Backbone framework
  - [ ] clustering
  - [ ] classification
- [ ] Model for representative sequence (probabalistic?)
- [ ] Given sequence, predict structure...
  - [ ] find best cluster
  - [ ] then find best sequence within the cluster

Outstanding Questions:
===
1. Classification: how to find a representative structure? --> RANSAC

Usage:
===
(How to use this applications/what people need to do in order to use this stuff...)

Data:
===
- 

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
