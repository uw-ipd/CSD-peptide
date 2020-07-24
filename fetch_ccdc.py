#!/usr/bin/env python

from ccdc import io
from ccdc.search import SubstructureSearch, MoleculeSubstructure, SMARTSSubstructure

# main
if __name__ == "__main__":
	strings = [
		"NCC(=O)NCC(=O)",
		"N(C)CC(=O)NCC(=O)",
		"NCC(=O)N(C)CC(=O)",
		"N(C)CC(=O)N(C)CC(=O)"
	]

	for pepsmile in strings:
		pep = SMARTSSubstructure(pepsmile)
		substructure_search_smiles = SubstructureSearch()
		sub_id = substructure_search_smiles.add_substructure(pep)
		print ("Searching...")
		hits_smiles = substructure_search_smiles.search(max_hits_per_structure=1)
		print (len(hits_smiles),"hits found")
		for hit in hits_smiles:
			#print (hit.identifier)
			with io.MoleculeWriter("from_CSD/nmeth/"+hit.identifier+".pdb") as pdb_writer:
				pdb_writer.write(hit.molecule)