#!/bin/bash

I=1
E=10 # How many folder/how many conformers 
i=$I

# I=1
# i=$I
path=`pwd`
files=$(ls $path/*.pdb)
for filename in $files
do
     outname="${filename%.pdb}.json"
     /home/ainan/x3dna-dssr -i=$filename -o=$outname --json --more
done
  
