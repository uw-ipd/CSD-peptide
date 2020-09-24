#!/bin/bash

for aa in \
    001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 \
    101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 \
    122 123 124 125 126 127 128 129 130 131 132 202 203 204 205 206 207 208 209 210 211 \
    302 303 304 305 306 307 308 309 313 314 315 316 317 318 320 332 333 401 402 403 404 \
    405 406 407 408 409 410 411 412 413 501 502 503 504 505 506 507 508 509 601 621 623 \
    631 633 701 702 703 704
do
	echo "> dummy" > dummy.fasta
	echo "X["$aa":NtermPeptoidFull:CtermPeptoidFull]" >> dummy.fasta
	~/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease \
		-in:file:fasta dummy.fasta \
		-parser::protocol generator.xml \
		-out:pdb -in:file:fullatom
	if [[ -f S_0001.pdb ]]; then
    	cat S_0001.pdb | grep -v OXT > ../$aa.pdb
	    rm S_0001.pdb
	fi
done