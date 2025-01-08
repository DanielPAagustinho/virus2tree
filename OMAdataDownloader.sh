#!/bin/bash
#SBATCH -A proj-gm0001

# Input file containing accession numbers
input_file=$1

# Check if the input file exists
if [ ! -f "$input_file" ]; then
	    echo "Input file not found: $input_file"
	        exit 1
	fi

mkdir -p db	
# Loop through each accession number in the input file
while IFS= read -r accession || [[ -n "$accession" ]]; do
    # Fetch CDS data in FASTA format
    efetch -db nucleotide -id "$accession" -format fasta_cds_na > db/"${accession}_cds_from_genomic.fna"

    #Convert to protein sequences.
    python /stornext/snfs4/next-gen/scratch/daniel/projects/RSV/p1442/r2t/OMAdb/protein_converter.py db/"${accession}_cds_from_genomic.fna" db/"${accession}_translated_cds.faa" 



    # Add a delay between requests (optional, to avoid overloading the server)
    sleep 1
done < "$input_file"

# Script downloaded from the R2T github: https://github.com/DessimozLab/read2tree/blob/main/archive/scripts/clean_fasta_cdna_cds.py#L12:~:text=adjust_mapping_names.py-,clean_fasta_cdna_cds,-.py
# usage is like: python clean_fasta_cdna_cds.py input_folder_fasta output_folder output_cdna.fa
python clean_fasta_cdna_cds.py db DB output_cdna.fa

