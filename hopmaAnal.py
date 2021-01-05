#!/usr/local/bin/python

import Bio.PDB
import numpy
import sys
import subprocess
from Bio.PDB.PDBIO import Select

############# FIX 1

copy=Bio.PDB.Atom.copy
def myCopy(self):
    shallow = copy.copy(self)
    for child in self.child_dict.values():
        shallow.disordered_add(child.copy())
    return shallow
Bio.PDB.Atom.DisorderedAtom.copy=myCopy

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_residue_dist_min(resid1, resid2) :
    """Returns the minimum distance between two residues"""
    distmin = 10000
    for atom1 in resid1:
    	for atom2 in resid2:
		    vec  = atom1.coord - atom2.coord
		    dist =  numpy.sqrt(numpy.sum(vec * vec))
		    if dist < distmin:
		    	distmin = dist
    return distmin 

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :
    	#print row, residue_one
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist_min(residue_one, residue_two)
    return answer

class caSelect(Select):
    def accept_atom(self, atom):
        if atom.get_name() == 'CA':
            return True
        else:
            return False

def computeDist(pdb_code, chain):
    """Computte a calpha-calpha distance matrix for a protein chain"""
    pdb_filename = pdb_code+'.pdb'
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_code, pdb_filename)
    first_model = structure[0]
    dist_matrix = calc_dist_matrix(first_model[chain], first_model[chain])
    numpy.savetxt(pdb_code+'.dist',dist_matrix)

def writeCA(pdb_code, chain):
    pdb_filename = pdb_code+'.pdb'
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_code, pdb_filename)
    first_model = structure[0]
    io = Bio.PDB.PDBIO()
    io.set_structure(first_model[chain])
    io.save(pdb_code+chain+'_CA.pdb', caSelect())
    subprocess.call('$HOPMA_PATH/PDBConverterIN.pl '+pdb_code+chain+'_CA',shell=True)

def preprocessInput(pdb_code, chain):
    writeCA(pdb_code, chain)
    computeDist(pdb_code+chain+"_CA", chain)

def createColoredMaps(pdb_code,chain,k,cutoff,dmin,dmax):
    subprocess.call('Rscript --vanilla $HOPMA_PATH/doCreateColoredMatSingle.R '+pdb_code+chain+' '+k+' '+cutoff+' '+dmin+' '+dmax,shell=True)

def extractExcludedContacts(pdb_code,chain,k,size,cutoff):
    subprocess.call('Rscript --vanilla $HOPMA_PATH/doExtractFromColoredMatSingle.R '+pdb_code+chain+' '+k+' '+size+' '+cutoff,shell=True)
