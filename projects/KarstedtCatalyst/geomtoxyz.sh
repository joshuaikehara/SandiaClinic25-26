#!/bin/bash
#SBATCH --job-name=geom_xyz
#SBATCH --output=geom_xyz.out
#SBATCH --error=geom_xyz.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

# Luis Lorenzana
# Harvey Mudd College
# Sandia National Laboratories Clinic 2025

# how to run: sbatch convert_geom2xyz.sh input.geom [output.xyz]
# if no output file is specified, the script will create one using the input filename with a .xyz extension

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 input.geom [output.xyz]"
    exit 1
fi

INPUT="$1"
if [ "$#" -eq 2 ]; then
    OUTPUT="$2"
else
    OUTPUT="${INPUT%.*}.xyz"
fi

FACTOR=0.529177  # conversion from bohr to angstrom
# count non-empty lines excluding header lines that contain "atom, type, position;"
NUM_ATOMS=$(grep -v 'atom, type, position;' "$INPUT" | grep -cv '^\s*$')

{
    # write the header for the .xyz file.
    echo "$NUM_ATOMS"
    echo ".geom to .xyz and scaled by factor of $FACTOR"

    # process the .geom file:
    # - Skip any line containing the header text
    # - If the line has 5 fields (index, symbol, x, y, z), drop the index
    # - Otherwise, assume the line contains at least 4 fields: symbol, x, y, z
    # - Multiply the coordinates by the factor
awk -v factor="$FACTOR" '
    {
        if ($0 ~ /atom, type, position;/) next; # Skip header

        if (NF==5 && $1 ~ /^[0-9]+$/) {          # Case 1: 5 fields, 1st is number
            symbol = $2;
            x = $3 * factor;
            y = $4 * factor;
            z = $5 * factor;
            printf "%s %.10f %.10f %.10f\n", symbol, x, y, z;
        } else if (NF>=5) {                     # Case 2: 5 or more fields, 1st is NOT number (like "AT1 Pt X Y Z")
            symbol = $2;
            x = $3 * factor;
            y = $4 * factor;
            z = $5 * factor;
            printf "%s %.10f %.10f %.10f\n", symbol, x, y, z;
        } else {                                # Case 3: Anything else (skip the line)
            next;
        }
    }' "$INPUT"
} > "$OUTPUT"  # <-- ADD THIS LINE to close the block and redirect output

echo "Conversion complete. Output written to $OUTPUT"
