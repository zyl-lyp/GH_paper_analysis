#!/bin/bash

####################################################################
# Author: Yupeng Liang
# Date: 2024-06-12
# Script Purpose: Convert SDF files to PDB, prepare receptors and ligands, run docking, and process results.
####################################################################

# Define the path to obabel
obabel=/pub/anaconda3/bin/obabel

# Convert SDF file to PDB file
$obabel -isdf pNPR.sdf -opdb -O ligand.pdb

# Initialize loop variables
i=1
index=1
end=20

# Create directory to store PDB files
mkdir -p pdb

# Get all subfolders in the parent directory
subfolders=($(ls -d ../*/))

# Main loop to process each subfolder
for folder in "${subfolders[@]}"; do
    # Strip the path and trailing slash from the folder name
    folder_name=$(basename "${folder%/}")

    # Prepare receptor and ligand files using the folder name as a prefix
    prepare_receptor -r ${folder_name}.pdb -o protein.pdbqt -A hydrogens
    prepare_ligand -l ligand.pdb -o ligand.pdbqt

    while (( i <= 50 )); do
        # Perform docking and save output
        vina --config config.txt --out dock_$i.pdbqt

        # Create subdirectory and move docking results
        mkdir -p dock_$i
        cp dock_$i.pdbqt dock_$i/

        # Enter subdirectory and split docking results
        cd dock_$i
        vina_split --input dock_$i.pdbqt
        cd ..

        # Inner loop to process each docking result
        while (( index <= end )); do
            t=$(printf "%02d" $((index + 20 - end)))
            echo $obabel -ipdbqt dock_$i/dock_${i}_ligand_$t.pdbqt -opdb -O pdb/${index}.pdb
            $obabel -ipdbqt dock_$i/dock_${i}_ligand_$t.pdbqt -opdb -O pdb/${index}.pdb
            index=$((index + 1))
        done

        # Update end index and main loop counter
        end=$((end + 20))
        i=$((i + 1))
    done
done

