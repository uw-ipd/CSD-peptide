#!/usr/bin/env python

import sys,os
import argparse
import numpy
import json
import itertools
import glob

import matplotlib.pyplot as plt
from scipy import ndimage

REFDIR = "refpep/canonic/"

atom_dtype = numpy.dtype( [
    ("atomname", numpy.unicode_, 4),
    ("resname", numpy.unicode_, 3),
    ("chnid", numpy.unicode_, 1),
    ("resid", numpy.int32),
    ("X", numpy.float32, 3),
    ("element", numpy.unicode_, 2),
] )

# pdb reader
def parse_pdb(pdbfile):
    allatoms = []
    for line in pdbfile:
        if line[:4] == 'ATOM' or line[:6] == "HETATM":
            resid = line[22:26]
            if (resid == '    '):
                resid = 0
            else:
                resid = int(resid)
                
            elt = line[76:78].strip()
            if (elt != 'H'):
                split_line = (
                    line[12:16], line[17:20], line[21], resid, 
                    (line[30:38], line[38:46], line[46:54]), elt
                )
                allatoms.append(split_line)
    return (numpy.array(allatoms, dtype=atom_dtype))

# pdb writer
def write_pdb(atoms, f=None):
    if (f is None):
        f = sys.stdout

    counter = 1
    for atom in atoms:
        f.write ("%-6s%5s %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(
            "ATOM", counter, atom['atomname'], atom['resname'], 
            atom['chnid'], atom['resid'], atom['X'][0], atom['X'][1], atom['X'][2],
            1.0, 0.0 ) )
        counter += 1


def get_torsion(X,Y,Z,W):
    def normalize(v):
        norm = numpy.linalg.norm(v)
        if norm == 0: 
           return v
        return v / norm

    A = normalize(Y-X)
    B = normalize(Z-Y)
    C = normalize(W-Z)

    angle = 0 # "undefined" value
    if ( numpy.linalg.norm(B) != 0 ):
        xval = -numpy.dot(A,C) + numpy.dot(A,B)*numpy.dot(B,C)
        yval = numpy.dot(A, numpy.cross(B,C))
        if (xval != 0 or yval != 0):
            angle = numpy.degrees( numpy.arctan2(yval,xval) )

    return angle


class PeptideParser:
    def __init__( self, minlen=0 ):
        self.minlen = minlen
        self.sc_names = []
        self.sc_codes = []
        self.sc_bfs_element = []
        self.sc_bfs_atomname = []
        self.sc_bfs_index = []
        self.sc_bonds = []
        self.sc_chiralities = []

        # parse reference sidechains
        for reffile in glob.glob(REFDIR+'*.pdb'):
            base=os.path.basename(reffile)
            stem=os.path.splitext(base)[0]
            #print (reffile)

            with open(reffile,'r') as ref:
                refatoms = parse_pdb( ref )
            refbonds = self.get_bondgraph( refatoms )
            refpeptide = self.get_peptides( refatoms, refbonds )
            assert (len(refpeptide) == 1)
            sc_atms, sc_bonds = self.get_sidechain( refpeptide[0], refatoms, refbonds, exhaustive=True )

            # AIB special case
            if (len(sc_atms[0]) <= 1 or stem == 'AIB'):
                chirality = 0
            else:
                chirality = int (numpy.sign( get_torsion(
                    refatoms['X'][refpeptide[0][0]],
                    refatoms['X'][refpeptide[0][2]],
                    refatoms['X'][refpeptide[0][1]],
                    refatoms['X'][sc_atms[0][1]]
                ) ) )

            for i in range(len(sc_atms)):
                self.sc_names.append( stem )
                self.sc_codes.append( refatoms[0]['resname'] )
                self.sc_bfs_element.append( [refatoms[j]['element'] for j in sc_atms[i]] )
                self.sc_bfs_atomname.append( [refatoms[j]['atomname'] for j in sc_atms[i]] )
                self.sc_bfs_index.append( refpeptide[0] + sc_atms[i][1:] )
                self.sc_bonds.append( sc_bonds[i] )
                self.sc_chiralities.append( chirality )

    def parse( self, atoms ):
        bonds = self.get_bondgraph( atoms )
        atoms,bonds = self.remove_unbound( atoms, bonds )
        peptides = self.get_peptides( atoms, bonds )
        
        if (len(peptides) < self.minlen):
            return None
        
        pep_indices = [atm for pep in peptides for atm in pep]
        
        mappedatoms = atoms.copy()
        atomorder = []

        for resid,pep in enumerate(peptides):
            # sidechain
            scatms, scbonds = self.get_sidechain( pep, atoms, bonds )
            scelements = [atoms[i]['element'] for i in scatms]

            if (len(scatms) <= 1):
                chirality = 0
            else:
                chirality = int (numpy.sign( get_torsion(
                    atoms['X'][pep[0]],
                    atoms['X'][pep[2]],
                    atoms['X'][pep[1]],
                    atoms['X'][scatms[1]]
                ) ) )

            nmethyl = self.get_nmethyl( pep, pep_indices+scatms, atoms, bonds )

            # 1 identify SC
            scindex = None
            for i in range(len(self.sc_names)):
                if (
                    self.sc_bfs_element[i] == scelements 
                    and self.sc_bonds[i] == scbonds
                    and (
                        self.sc_chiralities[i] == chirality 
                        or self.sc_chiralities[i] == 0 # AIB special case
                    )
                ):
                    scindex = i
                    break
                    
            if (scindex == None):
                return None

            mappedatoms['atomname'][pep] = [" N  "," CA ", " C  ", " O  "]
            mappedatoms['atomname'][scatms] = self.sc_bfs_atomname[scindex]

            if (nmethyl is not None):
                mappedatoms['atomname'][nmethyl] = " CN "
                atom_indices = pep + [nmethyl] + scatms[1:]
            else:
                atom_indices = pep + scatms[1:]
            #print (i,self.sc_codes[scindex],self.sc_bfs_index[scindex], scelements, scatms )
            atom_indices_ordered = [ atom_indices[x] for x in (0,1,2,3) ]  ## bb
            if (nmethyl is not None):
                atom_indices_ordered += [ atom_indices[4] ]   ## nmeth (CN)
                atom_indices_ordered += [ 
                    atom_indices[x+1] for x in self.sc_bfs_index[scindex][4:] ] ## sc
            else:
                atom_indices_ordered += [ 
                    atom_indices[x] for x in self.sc_bfs_index[scindex][4:] ] ## sc
            atomorder += atom_indices_ordered

            mappedatoms['resname'][atom_indices] = self.sc_codes[scindex]
            mappedatoms['resid'][atom_indices] = resid+1

        # ensure that nothing we've mapped is bonded to something we haven't mapped
        if (len(numpy.unique(atomorder)) != len(atoms)):
            return None

        return mappedatoms[atomorder]
        
    # bondgraph
    def get_bondgraph( self, atoms ):
        dists = numpy.linalg.norm( atoms['X'][:,None] - atoms['X'][None,:], axis=-1)
        bonds = ((dists>0.0) * (dists<2.0)).astype(numpy.bool)
        return bonds

    # remove disconnected components (water & ions)
    # maxconn is max connected components to remove
    def remove_unbound( self, atoms, bonds, maxconn=5 ):
        # recursive edge traversal
        def expand(toexpand, edges, visited):
            new_nodes = []

            for idx in range(len(edges[0])):
                i,j = edges[0][idx], edges[1][idx]
                if (i in toexpand and j not in visited):
                    new_nodes.append(j)
                    
            if (len(new_nodes) == 0):
                return visited
            else:
                return expand(new_nodes, edges, [*new_nodes,*visited])

        largest_subgraph = numpy.zeros( atoms.shape[0], dtype=numpy.bool )
        visited = numpy.zeros( atoms.shape[0], dtype=numpy.bool )
        maxsize=0
        
        # bfs
        for i,atm in enumerate (atoms):
            if (not visited[i]):
                graph_i = expand([i], bonds.nonzero(), [i])
                visited[graph_i] = True
                if (len(graph_i)>maxsize):
                    maxsize = len(graph_i)
                    largest_subgraph[:] = False
                    largest_subgraph[graph_i] = True
            
        mi = (largest_subgraph[:,None]*largest_subgraph[None,:]).nonzero()

        return atoms[largest_subgraph],bonds[mi].reshape(sum(largest_subgraph),-1)

    # peptide bonds
    def get_peptides( self, atoms, bonds ):
        # recursive path expansion
        def grow(direction, path, edges):
            longest_path, unassigned = path, edges
            for [i,j] in edges:
                new_path = None
                if (direction<0 and path[0]==j and i not in path):
                    new_path = [i,*path]
                elif (direction>0 and path[-1]==i and j not in path):
                    new_path = [*path,j]

                if (new_path is not None):
                    path_i,edges_i = grow(
                        direction, new_path, [e for e in edges if e != [i,j]])
                    if (len(path_i)>len(longest_path)):
                        longest_path = path_i
                        unassigned = edges_i
            return longest_path, unassigned

        peptides = []
        allbonds = bonds.nonzero()

        for idx in range(len(allbonds[0])):
            i,j = allbonds[0][idx], allbonds[1][idx]
            if (atoms[i]['element'] == 'N' and atoms[j]['element'] == 'C'):
                for k in bonds[j].nonzero()[0]:
                    if (atoms[k]['element'] == 'C'):
                        for l in bonds[k].nonzero()[0]:
                            if (atoms[l]['element'] == 'O'):
                                peptides.append([i,j,k,l])

        if (len(peptides) == 0):
            return []

        # find peptide bonds
        nlinks = numpy.zeros(len(peptides))
        for i,pep1 in enumerate(peptides):
            for j,pep2 in enumerate(peptides):
                if ( bonds[pep1[2],pep2[0]] ):
                    nlinks[i] += 1
                    nlinks[j] += 1

        # if there is overlap, favor the peptide with more connections
        peptide_dedup = []
        for i,pep1 in enumerate(peptides):
            keep_i = True
            for j,pep2 in enumerate(peptides):
                if (i==j):
                    continue
                if (len (set(pep1).intersection(pep2))>0):
                    if (nlinks[j] > nlinks[i] or
                        (nlinks[j] == nlinks[i] and i>j)
                    ):
                        keep_i = False
                        break
            if (keep_i):
                peptide_dedup.append(pep1)

        #for pep in peptide_dedup:
        #    print ([atoms[k]['atomname'] for k in pep])

        links = []
        for i,pep1 in enumerate(peptide_dedup):
            for j,pep2 in enumerate(peptide_dedup):
                if ( bonds[pep1[2],pep2[0]] ):
                    links.append((i,j))

        longest_path = [0]
        for (i,j) in links:
            path,unassigned = grow(
                1, [i,j], 
                [e for e in links if e != [i,j]])
            path,unassigned = grow(-1, path,unassigned)
            if (len(path) > len(longest_path)):
                longest_path = path
        linked_peptides = [peptide_dedup[i] for i in longest_path]

        return linked_peptides

    # get nmethylation
    def get_nmethyl( self, peptide, pepatoms, atoms, bonds ):
        link_atom = peptide[0]
        bonds_i = bonds[link_atom].nonzero()[0]
        methyls = numpy.setdiff1d(bonds_i,pepatoms)
        if (len(methyls) == 1 and atoms[methyls[0]]['element'] == 'C'):
            return (methyls[0])
        return None


    # get SC
    def get_sidechain( self, peptide, atoms, bonds, exhaustive=False ):
        def generate_perturbations(sc_atms, sc_atom_depth, sc_bonds):
            sc_atms_enum, sc_bonds_enum = [sc_atms],[sc_bonds]
            ndepths = max(sc_atom_depth)
            for depth in range(1,ndepths+1):
                atms_at_depth = [a for i,a in enumerate(sc_atms) if sc_atom_depth[i]==depth]
                elts_at_depth = [atoms[j]['element'] for j in atms_at_depth]

                toadd_sc_atms, toadd_sc_bonds = [],[]
                for scramble in itertools.permutations(atms_at_depth):
                    if (scramble == tuple(atms_at_depth)):
                        continue

                    scramble_elts = [atoms[j]['element'] for j in scramble]
                    if (elts_at_depth == scramble_elts):
                        for i in range(len(sc_atms_enum)):
                            new_sc_bonds = []
                            for [a,b] in sc_bonds_enum[i]:
                                mapa = a
                                if (sc_atms[a] in atms_at_depth):
                                    mapa = sc_atms.index( scramble[atms_at_depth.index(sc_atms[a])] )
                                mapb = b
                                if (sc_atms[b] in atms_at_depth):
                                    mapb = sc_atms.index( scramble[atms_at_depth.index(sc_atms[b])] )
                                new_sc_bonds.append([mapa,mapb])
                            new_sc_bonds.sort()
                            if (new_sc_bonds != sc_bonds):
                                new_sc_atms = [
                                    a if a not in atms_at_depth else scramble[atms_at_depth.index(a)]
                                    for a in sc_atms_enum[i]]
                                sc_atms_enum.append(new_sc_atms)
                                sc_bonds_enum.append(new_sc_bonds)

            return sc_atms_enum, sc_bonds_enum

        link_atom = peptide[1]
        sc_atms, sc_atom_depth, sc_bonds = [link_atom], [0], []

        # breadth-first traversal
        to_expand = [link_atom]
        round = 0
        while (len(to_expand) > 0):
            expanding = to_expand
            to_expand = []
            round += 1
            for i in expanding:
                bonds_i = bonds[i].nonzero()[0]
                for j in bonds_i:
                    if (not j in peptide and not j in sc_atms and not j in to_expand):
                        to_expand.append(j)
            reindex = numpy.argsort(atoms[to_expand]['element'])
            to_expand = [to_expand[x] for x in reindex]

            for j in to_expand:
                sc_atms.append(j)
                sc_atom_depth.append(round)

        Blist = bonds.nonzero()
        for i in range(len(Blist[0])):
            ai,aj = Blist[0][i],Blist[1][i]
            if (ai in sc_atms and aj in sc_atms):
                idxi = sc_atms.index(ai)
                idxj = sc_atms.index(aj)
                if (idxi<idxj):
                    sc_bonds.append( [idxi, idxj] )
            elif (ai in sc_atms and aj in peptide):
                idxi = sc_atms.index(ai)
                idxj = peptide.index(aj)
                sc_bonds.append( [idxi, -idxj] )

        sc_bonds.sort()
        
        if (exhaustive):
            sc_atms, sc_bonds = generate_perturbations(
                sc_atms, sc_atom_depth, sc_bonds)

        return sc_atms, sc_bonds
        
# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='ccdc converter')
    parser.add_argument('pdbs', nargs='+', help='input pdbs')
    parser.add_argument('--minlen', help='min peptide length', type=int, default=3)
    args = parser.parse_args()

    parser = PeptideParser( args.minlen )

    # 1) read pdb
    for pdb in (args.pdbs):
        print (">",pdb)
        base=os.path.basename(pdb)
        stem=os.path.splitext(base)[0]

        with open(pdb) as pdbfile:
            atoms = parse_pdb( pdbfile )
        fixedatoms = parser.parse(atoms)
        if (fixedatoms is not None):
            filename = stem+"_fixed.pdb"
            with open(filename,'w') as outfile:
                write_pdb( fixedatoms, outfile )
        else:
            print("no good")

            
