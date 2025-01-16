#!/bin/bash
#SBATCH -A proj-gm0001
set -e 
set -u
set -o pipefail

MAIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${MAIN_DIR}/../scripts"
OMA="${MAIN_DIR}/../oma/bin"
PARAMETERS_FILE="parameters.drw"
FIVE_LETTER_FILE="clean_five_letter_species.tsv"
read_type="long-hifi"  # Default read type
threads=12
outgroup_file=""
reads_dir=$(pwd)  # Default to current working directory


function show_help() {
    echo "Usage: $0 -i <input_file> [-o <outgroup_file>] [-t <read_type>] [-n <threads>] [-r <reads_dir>] [-h]"
    echo
    echo "Required arguments:"
    echo "  -i, --input        Path to the accession file (required)."
    echo
    echo "Optional arguments:"
    echo "  -o, --outgroup     Path to the outgroup accession file."
    echo "  -t, --read-type    Type of sequencing reads ('long-ont', 'long-hifi', or 'short'). Default: 'long-hifi'."
    echo "  -n, --threads      Number of threads to use for parallel processing. Default: 12."
    echo "  -r, --reads-dir    Directory containing the read files. Default: current working directory."
    echo "  -h, --help         Display this help message."
}

# Parse flags
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) input_file="$2"; shift ;;
        -o|--outgroup) outgroup_file="$2"; shift ;;
        -t|--read-type) read_type="$2"; shift ;;
        -n|--threads) threads="$2"; shift ;;
        -r|--reads-dir) reads_dir="$2"; shift ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
    esac
    shift
done
# Check required parameters
if [[ -z "${input_file:-}" ]]; then
    echo "Error: --input (-i) is required."
    show_help
    exit 1
fi

# Validate input file existence
if [[ ! -f "$input_file" ]] || [[ ! -s "$input_file" ]]; then
    echo "Error: The input file '$input_file' is missing or empty."
    exit 1
fi

# Validate outgroup file if provided
if [[ -n "$outgroup_file" ]] && ([[ ! -f "$outgroup_file" ]] || [[ ! -s "$outgroup_file" ]]); then
    echo "Error: The outgroup file '$outgroup_file' is missing or empty."
    exit 1
fi

if [[ ! -d "$reads_dir" || -z "$(ls -A "$reads_dir" 2>/dev/null)" ]]; then
    echo "Error: The reads directory '$reads_dir' does not exist or is empty."
    exit 1
fi

mkdir -p db

echo "Fetching data from accession numbers"	
while IFS= read -r accession || [[ -n "$accession" ]]; do
    efetch -db nucleotide -id "$accession" -format fasta_cds_na > db/"${accession}_cds_from_genomic.fna"
    
    #datasets download virus genome accession "${accession}.1" --include cds --filename "${accession}.zip" && \
    #unzip -p "${accession}.zip" "ncbi_dataset/data/cds.fna" > db/"${accession}_cds_from_genomic.fna" && \
    #rm "${accession}.zip"
    
    #Convert to protein sequences.
    #python "${SCRIPTS_DIR}"/protein_converter.py db/"${accession}_cds_from_genomic.fna" db/"${accession}_translated_cds.faa" 
    sleep 1
done < "$input_file"

echo "Finished retrieval and conversion of nucleotide sequences"

#Adding 5 letter code and cleaning for r2t
python "${SCRIPTS_DIR}"/clean_fasta_cdna_cds.py db
cp -r clean_aa DB 

echo "Running OMA standalone"
"${OMA}/oma" -p  

outgroup_codes=()
if [ -n "$outgroup_file" ]; then
    echo "Reading outgroup accessions from $outgroup_file..."
    while IFS= read -r accession || [[ -n "$accession" ]]; do
        if grep -q "^${accession}_" "$FIVE_LETTER_FILE"; then
            code=$(grep "^${accession}_" "$FIVE_LETTER_FILE" | awk '{print $2}')
            outgroup_codes+=("$code")
        else
            echo "Warning: The accession $accession was not found in $FIVE_LETTER_FILE."
        fi
    done < "$outgroup_file"
fi

# Verificar los códigos obtenidos y actualizar el parámetro OutgroupSpecies
if [ ${#outgroup_codes[@]} -gt 0 ]; then
    outgroup_list=$(IFS=,; echo "${outgroup_codes[*]}")
    OUTGROUPS="OutgroupSpecies := [${outgroup_list}];"
    echo "Using the following outgroup codes: ${outgroup_list}."
else
    echo "No valid outgroup accessions provided. Using default value 'none'."
    OUTGROUPS="OutgroupSpecies := 'none';"
fi
#map initial filen with the new 5 letter code given by the clean...py
grep -q "^OutgroupSpecies" "$PARAMETERS_FILE" && \
    sed -i "s/^OutgroupSpecies.*/$OUTGROUPS/" "$PARAMETERS_FILE"

echo "Updated parameters file with outgroups: $OUTGROUPS."

"${OMA}/oma"

if ls Output/OrthologousGroupsFasta/*.fa >/dev/null 2>&1; then
    mkdir -p marker_genes
    #cat Output/OrthologousGroupsFasta/*.fa > dna_ref.fa
    mv Output/OrthologousGroupsFasta/*.fa marker_genes
else
    error "No files found in Output/OrthologousGroupsFasta."
    exit 1
fi

echo "Running Read2Tree with "${threads}" threads"

#read2tree --standalone_path ./marker_genes --output_path read2tree_output --dna_reference dna_ref.fa
read2tree --step 1marker --standalone_path marker_genes --dna_reference dna_ref.fa --output_path read2tree_output --debug 


# Procesar archivos según el tipo de lectura
#if [[ "$read_type" == "long-ont" || "$read_type" == "long-hifi" ]]; then
#   reads=($(find "$reads_dir" -name "*.fastq.gz" -exec basename {} .fastq.gz \;))  # Archivos de lecturas largas
#   for read in "${reads[@]}"; do
        #read2tree --step 2map --standalone_path marker_genes --dna_reference dna_ref.fa --reads "${read}.fastq.gz" --read_type "$read_type" --threads "$threads" --output_path read2tree_output --debug
        #echo "Processing long read file: ${read}.fastq.gz"
    #done
#elif [[ "$read_type" == "short" ]]; then
    # Archivos de lecturas cortas emparejadas (R1 y R2)
    #reads=($(find "$reads_dir" -name "*_1.fastq.gz" -exec basename {} _1.fastq.gz \;))  # Detectar muestras basadas en R1
    #for read in "${reads[@]}"; do
        #R1="${read}_1.fastq.gz"
        #R2="${read}_2.fastq.gz"
        #if [[ -f "$R1" && -f "$R2" ]]; then
            #read2tree --step 2map --standalone_path marker_genes --dna_reference dna_ref.fa --reads "$R1" "$R2" --read_type "$read_type" --threads "$threads" --output_path read2tree_output --debug
            #echo "Processing paired-end reads: $R1 and $R2"
        #else
            #echo "Warning: Missing pair for $read. Skipping..."
        #fi
    #done
#else
    #echo "Error: Unsupported read type '$read_type'. Use 'long-ont', 'long-hifi', or 'short'."
    #exit 1
#fi



#echo "Merging..."
#read2tree --step 3combine --standalone_path marker_genes --dna_reference dna_ref.fa --output_path read2tree_output --tree --debug
#iqtree -T ${threads} -s read2tree_output/concat_*_aa.phy -bb 1000
#iqtree -T ${threads} -s read2tree_output/concat_*_dna.phy 