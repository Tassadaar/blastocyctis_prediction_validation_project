#!/bin/bash

gff_file=$1

# Extract the ID attribute values for genes from the GFF3 file
id_values=$(awk -F'\t' '!/^#/ && $3 == "gene" {split($9,a,";"); for(i in a) if(a[i] ~ /^ID=/) print a[i]}' "$gff_file")

# Check for duplicate ID values
duplicate_ids=$(echo "$id_values" | sort | uniq -d)

# Print the duplicate ID values
if [ -n "$duplicate_ids" ]; then
  echo "Duplicate ID values found for genes:"
  echo "$duplicate_ids"
else
  echo "No duplicate ID values found for genes."
fi

