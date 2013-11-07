LSP
===
We can include a general description of our project here...

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


Dependencies:
===
- Python 2.7
- Biopython, see https://github.com/biopython/biopython
- DSSP, see http://swift.cmbi.ru.nl/gv/dssp/ (If you are running Mac OSX and are having trouble getting DSSP installed, I suggest this resource: http://proteinz.blogspot.com/2013/02/compiling-dssp-on-osx-lion-redux.html).
