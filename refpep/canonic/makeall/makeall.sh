#!/bin/bash

for aa in \
  ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL \
  DALA DARG DASN DASP DCYS DGLN DGLU DHIS DILE DLEU DLYS DMET DPHE DPRO DSER DTHR DTRP DTYR DVAL \
  AIB
do
	echo "> dummy" > dummy.fasta
	echo "X["$aa"]" >> dummy.fasta
	~/Rosetta_master/main/source/bin/rosetta_scripts.default.macosclangrelease \
		-in:file:fasta dummy.fasta \
		-parser::protocol generator.xml \
		-out:pdb -in:file:fullatom
	mv S_0001.pdb ../$aa.pdb
done