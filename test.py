'''
Created on Nov 2, 2013

@author: ramos
'''
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

if __name__ == "__main__":
    dssp = dssp_dict_from_pdb_file("pdb/membrane/1a0sQ.pdb")
    for key in dssp[1]:
        print(str(key) + ":" + str(dssp[0][key]))
    print("(")
    for dict in dssp:
        print("\t" + "{")
        try:
            for key in dict.keys():
                print("\t\t" + str(key) + ":" + str(dict[key]))
        except:
            for key in dict:
                print("\t\t" + str(key))
        print("\t" + "}")
    print(")")    