# Vcontacts

Commandline tool to compute surface areas in contact using the constrained voronoi procedure from the paper: Quantification of protein surfaces, volumes and atom-atom contacts using a constrained Voronoi procedure. (doi: 10.1093/bioinformatics/18.10.1365)

To compile:
gcc -c Vcontacts-v1-2.c
gcc Vcontacts-v1-2.o -o vcontacts

To run:
./vcontacts yourpdbfile.pdb
