#!/bin/bash
set -e 
set -u
set -o pipefail

if [ -t 1 ]; then
  RED="\033[1;31m"
  GREEN="\033[1;32m"
  YELLOW="\033[1;33m"
  NC="\033[0m"
else
  RED=""
  GREEN=""
  YELLOW=""
  NC=""
fi

log_info() {
  echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]${NC} $*"
}

log_warn() {
  echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]${NC} $*" >&2
}

log_error() {
  echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR]${NC} $*" >&2
}

MAIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${MAIN_DIR}/../scripts"
OMA="${MAIN_DIR}/../oma/bin"
PARAMETERS_FILE="parameters.drw"
FIVE_LETTER_FILE="clean_five_letter_species.tsv"

show_help() {
  cat << EOF
Usage: $0 [options]

Required:
  -i, --input        Path to the accession file (comma-separetad file with the name of species/strain as its first column)

Optional:
  -o, --outgroup     Path to the outgroup species/strain file
  -t, --read-type    Type of reads: long-ont, long-hifi, or short [default: long-hifi]
  -n, --threads      Number of threads [default: 12]
  -r, --reads-dir    Directory with read files [default: current directory]
  -h, --help         Show this help message

Example:
  $0 -i accessions.txt -o outgroups.txt -t short -r /path/to/reads
EOF
}

fetch_data() {
  local line="$1"

  # Extract strain name (first column) and accession(s) (remaining columns)
  local strain
  local accessions
  strain="$(echo "$line" | cut -d ',' -f1 | tr -d ' ')"
  accessions="$(echo "$line" | cut -d ',' -f2- | tr -d ' ')"

  # Skip empty lines
  [[ -z "$strain" || -z "$accessions" ]] && exit 0

  log_info "Fetching data for strain: $strain (Accessions: $accessions)"
  efetch -db nucleotide -id "$accessions" -format fasta_cds_na > "db/${strain}_cds_from_genomic.fna" || {
    log_error "Failed to fetch accession(s) for ${strain}: ${accessions}"
    exit 1
  }

  log_info "Done fetching data for: $strain"
}

# Function to process files and generate the output table
generate_og_gene_tsv() {
    local fna_dir=$1     # Directory containing .fna files downloaded from NCBI using efetch
    local fa_dir=$2      # Directory containing .fa files from OMA output (Output/OrthologousGroupsFasta)
    local output_file=$3 # Output file to store results

    # Example usage
    # process_genes ~/oma/test3_illumina/db ~/oma/test3_illumina/marker_genes output_table.tsv
    #temp file
    local tmp_file
    tmp_file="$(mktemp)"
    while IFS=: read -r file line; do
        local gene
        gene="$(echo "$line" | grep -oP '\[gene=\K[^\]]+')"
        local protein
        protein="$(echo "$line" | grep -oP '\[protein=\K[^\]]+')"
        local protein_id
        protein_id="$(echo "$line" | grep -oP '\[protein_id=\K[^\]]+')"
        local location
        location="$(echo "$line" | grep -oP '\[location=\K[^\]]+')"

        local compressed_id
        compressed_id="$(echo "$protein_id" | sed 's/[^a-zA-Z0-9]//g')"
        local match_found=false
        while IFS=: read -r og_file og_line; do
            local og
            og="$(basename "$og_file" .fa)"
            echo -e "${og}\t${gene}\t${protein}\t${protein_id}\t${location}" >> "$tmp_file"
            match_found=true
        done < <(grep '^>' "$fa_dir"/*.fa | grep "$compressed_id")
        if ! $match_found; then
            log_warn "OG-Gene TSV generation: No match found for Protein ID: $protein_id"
        fi
    done < <(grep '^>' "$fna_dir"/*.fna)
    #-V option in sort from GNU coreutils
    sort -k1,1 -V -k2,2 "$tmp_file" > "$output_file"
    sed -i '1i OG\tGene\tProtein\tProtein_ID\tLocation' "$output_file"
    rm -f "$tmp_file"

    log_info "OG-Gene TSV generation complete: $output_file"
}
READ_TYPE="long-hifi"  # Default read type
THREADS=12
OUTGROUP_FILE=""
READS_DIR=$(pwd)  # Default to current working directory

####################################################

#MAIN

###################################################

# Parse flags
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_FILE="$2"; shift ;;
        -o|--outgroup) OUTGROUP_FILE="$2"; shift ;;
        -t|--read-type) READ_TYPE="$2"; shift ;;
        -n|--THREADS) THREADS="$2"; shift ;;
        -r|--reads-dir) READS_DIR="$2"; shift ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
    esac
    shift
done
# Check required parameters
if [[ -z "${INPUT_FILE:-}" ]]; then
    echo "Error: --input (-i) is required."
    show_help
    exit 1
fi

# Validate input file existence
if [[ ! -f "$INPUT_FILE" ]] || [[ ! -s "$INPUT_FILE" ]]; then
    echo "Error: The input file '$INPUT_FILE' does not exist or is empty."
    exit 1
fi

# Validate outgroup file if provided
if [[ -n "$OUTGROUP_FILE" ]] && ([[ ! -f "$OUTGROUP_FILE" ]] || [[ ! -s "$OUTGROUP_FILE" ]]); then
    echo "Error: The outgroup file '$OUTGROUP_FILE' is missing or empty."
    exit 1
fi

if [[ ! -d "$READS_DIR" || -z "$(ls -A "$READS_DIR" 2>/dev/null)" ]]; then
    echo "Error: The reads directory '$READS_DIR' does not exist or is empty."
    exit 1
fi

mkdir -p db
log_info "Starting parallel retrieval of sequences from $INPUT_FILE..."
# Export the function so GNU Parallel can see it
export -f fetch_data log_info log_warn log_error
export RED YELLOW BLUE GREEN NC

grep -v '^$' "$INPUT_FILE" | parallel -j 3 fetch_data {}
log_info "Finished retrieval of nucleotide sequences."

#Adding 5 letter code and cleaning for r2t
python "${SCRIPTS_DIR}"/clean_fasta_cdna_cds.py db
cp -r clean_aa DB 
#remove then clean_aa?
log_info "Editing parameters file..."
"${OMA}/oma" -p 
#Editing parameters file
outgroup_codes=()
if [ -n "$OUTGROUP_FILE" ]; then
    log_info "Reading outgroup species/strain from $OUTGROUP_FILE..."
    while IFS= read -r species || [[ -n "$species" ]]; do
        if grep -q "^${species}_cds_" "$FIVE_LETTER_FILE"; then
            code=$(grep "^${species}_" "$FIVE_LETTER_FILE" | awk '{print $2}')
            outgroup_codes+=("$code")
        else
            log_warn "Outgroup species '$species' not found in $FIVE_LETTER_FILE."
        fi
    done < "$OUTGROUP_FILE"
fi

#Verify codes and edit the parameters file
if [ ${#outgroup_codes[@]} -gt 0 ]; then
    outgroup_list=$(IFS=,; echo "${outgroup_codes[*]}")
    OUTGROUPS="OutgroupSpecies := [${outgroup_list}];"
    log_info "Using the following 5 letter code outgroup(s): ${outgroup_list}."
else
    log_warn "No valid outgroup species provided. Using default 'none'."
    OUTGROUPS="OutgroupSpecies := 'none';"
fi
#map initial filen with the new 5 letter code given by the clean...py
grep -q "^OutgroupSpecies" "$PARAMETERS_FILE" && \
    sed -i "s/^OutgroupSpecies.*/$OUTGROUPS/" "$PARAMETERS_FILE"

log_info "Running OMA..."
"${OMA}/oma" -n ${THREADS}

if ls Output/OrthologousGroupsFasta/*.fa >/dev/null 2>&1; then
    #Create tsv
    generate_og_gene_tsv db Output/OrthologousGroupsFasta OG_genes.tsv
    mkdir -p marker_genes
    #cat Output/OrthologousGroupsFasta/*.fa > dna_ref.fa
    mv Output/OrthologousGroupsFasta/*.fa marker_genes
else
    log_error "No files found in Output/OrthologousGroupsFasta."
    exit 1
fi


log_info "Running Read2Tree (step 1marker) with ${THREADS} threads..."

#read2tree --standalone_path ./marker_genes --output_path read2tree_output --dna_reference dna_ref.fa
read2tree --step 1marker --standalone_path marker_genes --dna_reference dna_ref.fa --output_path read2tree_output --debug 


# Procesar archivos según el tipo de lectura
#if [[ "$READ_TYPE" == "long-ont" || "$READ_TYPE" == "long-hifi" ]]; then
#   reads=($(find "$READS_DIR" -name "*.fastq.gz" -exec basename {} .fastq.gz \;))  # Archivos de lecturas largas
#   for read in "${reads[@]}"; do
        #read2tree --step 2map --standalone_path marker_genes --dna_reference dna_ref.fa --reads "${read}.fastq.gz" --READ_TYPE "$READ_TYPE" --THREADS "$THREADS" --output_path read2tree_output --debug
        #echo "Processing long read file: ${read}.fastq.gz"
    #done
#elif [[ "$READ_TYPE" == "short" ]]; then
    # Archivos de lecturas cortas emparejadas (R1 y R2)
    #reads=($(find "$READS_DIR" -name "*_1.fastq.gz" -exec basename {} _1.fastq.gz \;))  # Detectar muestras basadas en R1
    #for read in "${reads[@]}"; do
        #R1="${read}_1.fastq.gz"
        #R2="${read}_2.fastq.gz"
        #if [[ -f "$R1" && -f "$R2" ]]; then
            #read2tree --step 2map --standalone_path marker_genes --dna_reference dna_ref.fa --reads "$R1" "$R2" --READ_TYPE "$READ_TYPE" --THREADS "$THREADS" --output_path read2tree_output --debug
            #echo "Processing paired-end reads: $R1 and $R2"
        #else
            #echo "Warning: Missing pair for $read. Skipping..."
        #fi
    #done
#else
    #echo "Error: Unsupported read type '$READ_TYPE'. Use 'long-ont', 'long-hifi', or 'short'."
    #exit 1
#fi



#echo "Merging..."
#read2tree --step 3combine --standalone_path marker_genes --dna_reference dna_ref.fa --output_path read2tree_output --tree --debug
#iqtree -T ${THREADS} -s read2tree_output/concat_*_aa.phy -bb 1000
#iqtree -T ${THREADS} -s read2tree_output/concat_*_dna.phy 