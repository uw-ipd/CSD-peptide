# CSD-peptide

This repository contains several utilities for finding and processing peptide structures deposited in the CDC small molecule database.  Two main scripts are included:
 * _fetch_ccdc.py_ : A script for finding all possible peptide-containing structures in the CSD.  Requires CSD and the python bindings to be installed.
 * _convert.py_ : A script that a) identifies all CSD structures that contain peptides, and b) writes them as "PDB-style"
 
### Notes on _convert.py_
 * Usage: _convert.py_ [--minlen MINLEN] pdbs
 * the "MINLEN" argument only generates peptides of at least the given length
 * the folder _refpep/canonic_ contains all the PDBs of which the script is aware.  Changing this to _refpep/all_ (or adding PDBs to this folder) allow understanding of additional PDBs.
 * outputs are found in _converted_
 * currently only alpha AAs and N-methylated AAs are understood.  Additional variants (peptoids) will be added in the future.
