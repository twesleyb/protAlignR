#!/bin/bash

protA="Q3UXZ6" # Fam81a
protB="Q9JK97" # Kcnq4

# Download protein sequences.
echo "Downloading fasta sequences for $protA and $protB..."
wget "https://www.uniprot.org/uniprot/$protA.fasta" --quiet
wget "https://www.uniprot.org/uniprot/$protB.fasta" --quiet

# Run water to align the sequences.
water $protA.fasta $protB.fasta -auto -stdout

# Remove temporary files.
rm $protA.fasta
rm $protB.fasta
