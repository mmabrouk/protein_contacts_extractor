import os, sys
from optparse import OptionParser
from Bio.PDB.PDBParser import PDBParser
from math import sqrt
#Calculate Euclidean distance between two points in 3D space
def dist(point1, point2):
   return sqrt( (point1[0] - point2[0])*(point1[0] - point2[0]) + (point1[1] - point2[1])*(point1[1] - point2[1]) + (point1[2] - point2[2])*(point1[2] - point2[2]) )
   
#Check if a two residues are in long-range or medium-range
def isFromType(res2, res1, c_type):
    if c_type == "long":
        return (res2 - res1 > 23)
    if c_type == "medium":
        return ( (res2-res1 >=12) and (res2-res1<=23) )
    else:
        print "Wrong parameter contactTypes. Choose between 'long' or 'medium'"
        sys.exit()
        return False

#Load pdb file and return object of the structure
def loadPDB(pdbFile):
    p=PDBParser(PERMISSIVE=1)
    structure_id="native"
    return p.get_structure(structure_id, pdbFile)

def main():
    parser = OptionParser()
    parser.add_option('-p', '--pdb', dest='pdb', help="The pdb file of the native structure")
    parser.add_option('-c', '--contactTypes', dest='contactType', default ="long", help = 'The type of the contacts (long, medium)')
    parser.add_option('-t', '--threshold', dest='threshold', default=8, help="Distance threshold for contacts (default 8)")
    (option, args) = parser.parse_args()
    pdbFile = option.pdb
    c_type = option.contactType
    thresh = option.threshold
    
    struct_copy = {}
    structure = loadPDB(pdbFile)
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                if residue.has_id("CA"):
                    ca=residue["CA"]
                    struct_copy[residue.get_id()[1]] = ca.get_coord()
                    #print residue.get_id()[1], ca.get_coord()
                    
    for res1 in struct_copy:
        for res2 in struct_copy:
            if (dist(struct_copy[res1], struct_copy[res2]) < thresh) and (isFromType(res2, res1, c_type)):
                print str(res1)+" "+str(res2)
                
if __name__ == "__main__":
    main()