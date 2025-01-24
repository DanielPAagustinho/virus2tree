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

show_help() {
  cat << EOF
Usage: $0 [options]

Required:
  -i, --input        Path to the accession file (comma-separetad file with the structure: taxon(s)/species/strain(s),code(optional),accession(s)). Note: If provided, the code should have exactly 5 alphanumeric characters.

Optional:
  -o, --outgroup     Path to the outgroup taxon/species/strain file
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
  local clean_line=$(echo "$line" | tr -d '[:space:]')
  IFS=',' read -ra columns <<< "$clean_line"
  local strain="${columns[0]}"
  if [[ -z "$strain" ]]; then
    log_warn "Skipping line with empty strain: $line"
    return 0
  fi
  # Decide if we have a code column or not
  local code=""
  local accessions_list=""
  if $HAS_CODE_COLUMN; then
    # The second column must be CODE
    if [[ "${#columns[@]}" -lt 3 ]]; then
      log_error "Line has fewer than 3 columns but we expected CODE. Line: $line"
      exit 1
    fi
    code="${columns[1]}"
    # The rest (from columns[2] onward) are accession(s)
    accessions_list=$(IFS=','; echo "${columns[@]:2}")
  else
    # No code column: second column onward => accessions
    if [[ "${#columns[@]}" -lt 2 ]]; then
      log_error "Line has fewer than 2 columns but we expected STRAIN,ACCESSIONS. Line: $line"
      exit 1
    fi
    accessions_list=$(IFS=','; echo "${columns[@]:1}")
  fi

  # Write the code map
  if $HAS_CODE_COLUMN; then
    echo "WRITING TO ${FIVE_LETTER_FILE}"
    echo -e "${strain}\t${code}" >> "${FIVE_LETTER_FILE}"
  fi

  if [[ "$accessions_list" == *GCF_* || "$accessions_list" == *GCA_* ]]; then
    local assembly_accessions=()
    local regular_accessions=()
    #From list to array
    IFS=',' read -ra accessions <<< "$accessions_list"

    for acc in "${accessions[@]}"; do
      if [[ $acc == GCF_* || $acc == GCA_* ]]; then
        assembly_accessions+=("$acc")
      else
        regular_accessions+=("$acc")
      fi
    done

    local nc_accessions=()

    for assembly in "${assembly_accessions[@]}"; do
      echo "Processing assembly: ${assembly} from ${strain}"
      nc_accessions+=($(esearch -db assembly -query "$assembly" \
        | elink -target nuccore \
        | efetch -format acc \
        | grep "^NC_")) #allow NP_ or NM?
      if [[ ${#nc_accessions[@]} -eq 0 ]]; then
        log_error "No NC_ accessions found for assembly: $assembly"
        exit 1 
      fi
    done
    #To avoid problem with commas if the regular array is void?
    if [[ -z "${regular_accessions[*]}" ]]; then
      accessions_list=$(IFS=','; echo "${nc_accessions[*]}")
    else
      accessions_list=$(IFS=','; echo "${nc_accessions[*]},${regular_accessions[*]}")
    fi
    echo "Final accesion list has ${#regular_accessions[@]} given accessions and ${#nc_accessions[@]} assembly accessions. The assembly accessions are:${nc_accessions[@]}"
  fi
  # Fetch data
  log_info "Fetching data for strain: $strain (Accessions: $accessions_list)"
  efetch -db nucleotide -id "$accessions_list" -format fasta_cds_na \
    > "db/${strain}_cds_from_genomic.fna" \
    || {
      log_error "Failed to fetch accession(s) for ${strain}: ${accessions_list}"
      exit 1
    }
  sleep 1
  log_info "Done fetching data for: $strain"
}
#declare -A SPECIES_TO_CODE
# Function to process files and generate the output table
generate_og_gene_tsv() {
    local fna_dir=$1     # Directory containing .fna files downloaded from NCBI using efetch
    local fa_dir=$2      # Directory containing .fa files from OMA output (Output/OrthologousGroupsFasta)
    local output_file=$3 # Output file to store results
    local unique_output_file="${output_file%.tsv}-unique.tsv"
    local tmp_file  
    local -A SPECIES_TO_CODE
    while IFS=$'\t' read -r species_val code_val; do
          echo "HEREEE, ${FIVE_LETTER_FILE}, ${species_val} and ${code_val}"
          SPECIES_TO_CODE["$species_val"]="$code_val"
    done < "$FIVE_LETTER_FILE"
    # Example usage
    # process_genes ~/oma/test3_illumina/db ~/oma/test3_illumina/marker_genes output_table.tsv
    #temp file

    tmp_file="$(mktemp)"
    while IFS=: read -r file line; do
        #Example of the complete input line: db/rsv_11_cds_from_genomic.fna:>lcl|MG813984.1_cds_AZQ19553.1_6 [gene=SH] [protein=small hydrophobic protein] [protein_id=AZQ19553.1] [location=4251..4445] [gbkey=CDS]
        local species=$(basename "$file" | awk -F '_cds_' '{print $1}')
        local code="${SPECIES_TO_CODE[$species]}"
        local accession="$(awk -F '_cds_' '{print $1}' <<< "${line#*|}")"
        #local protein_id2=$(awk -F '_cds_' '{print $2}' <<< "${line#*|}")

        local gene="NA"
        local protein="NA"
        local protein_id="NA"
        local location="NA"
        #local accession="NA"; local code="NA"
        #local gene="$(echo "$line" | grep -oP '\[gene=\K[^\]]+')"
        #local protein="$(echo "$line" | grep -oP '\[protein=\K[^\]]+')"
        #local protein_id="$(echo "$line" | grep -oP '\[protein_id=\K[^\]]+')"
        #local location="$(echo "$line" | grep -oP '\[location=\K[^\]]+')"

        # gene
        if [[ "$line" =~ \[gene=([^]]+)\] ]]; then
          gene="${BASH_REMATCH[1]}"
        fi

        # protein
        if [[ "$line" =~ \[protein=([^]]+)\] ]]; then
          protein="${BASH_REMATCH[1]}"
        fi

        # protein_id
        if [[ "$line" =~ \[protein_id=([^]]+)\] ]]; then
          protein_id="${BASH_REMATCH[1]}"
        fi

        # location
        if [[ "$line" =~ \[location=([^]]+)\] ]]; then
          location="${BASH_REMATCH[1]}"
        fi

        #Compressed_id uses the protein_id without alphanumeric characters to match the header in OrthologousGroupsFasta:
        #Output/OrthologousGroupsFasta/OG11.fa:>3355X||rsv_11||lcl|MG8139841cdsAZQ1955316 cleaned for r2t [rsv_11_3355X]
        local compressed_id="$(echo "$protein_id" | sed 's/[^a-zA-Z0-9]//g')"
        #echo "Got correct features from *fna for ${protein_id}"

        [[ "$gene" == "NA" ]]       && log_warn "Missing gene for line in $file"
        [[ "$protein" == "NA" ]]    && log_warn "Missing protein for line in $file"
        [[ "$protein_id" == "NA" ]] && log_warn "Missing protein_id for line in $file"
        [[ "$location" == "NA" ]]  && log_warn "Missing location for line in $file"
        # Buscar OG directamente con grep
        local og_matches
        if og_matches=$(grep -l "$compressed_id" "$fa_dir"/*.fa 2>/dev/null | xargs -I {} basename {} .fa); then
              echo "Matches found for $compressed_id: $og_matches"
              # Iterar sobre los OGs encontrados y escribir al archivo temporal
              for og in $og_matches; do
                  echo -e "${og}\t${gene}\t${protein}\t${protein_id}\t${location}\t${accession}\t${species}\t${code}" >> "$tmp_file"
              done
        else
            log_warn "Error: grep command failed while searching for compressed ID: $protein_id in directory: $fa_dir"
        fi
    done < <(grep '^>' "$fna_dir"/*.fna)
    #-V option in sort from GNU coreutils
    sort -k1,1 -V -k2,2 "$tmp_file" > "$output_file"
    sed -i '1i OG\tGene\tProtein\tProtein_ID\tLocation\tAccession\tTaxon\tCode' "$output_file"
        {
      echo -e "OG\tGene\tProtein\ttaxon"
      # Skip the header of the main TSV
      tail -n +2 "$output_file" \
      | awk -F'\t' '{
          # Key is combination of columns OG,Gene,Protein
          key=$1"\t"$2"\t"$3
          # We append species/strain from col 7
          if(!seen[key]) {
             order[++cnt]=key
             taxa[key]=$7
             seen[key]=$7
          } else {
             # Only add if not already present
             split(taxa[key], check, ";")
             found=0
             for(i in check) {
               if(check[i] == $7) {found=1; break}
             }
             if(!found) { taxa[key]=taxa[key]";"$7 }
          }
      }
      END {
          for(i=1;i<=cnt;i++){
             k=order[i]
             print k"\t"taxa[k]
          }
      }'
    } > "$unique_output_file"

    rm -f "$tmp_file"

    log_info "OG-Gene TSV generation complete: $output_file"
}
READ_TYPE="long-hifi"  # Default read type
THREADS=12
FIVE_LETTER_FILE="five_letter_taxon.tsv"
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
        *) log_error "Unknown parameter passed: $1"; show_help; exit 1 ;;
    esac
    shift
done
# Check required parameters
if [[ -z "${INPUT_FILE:-}" ]]; then
    log_error "Error: --input (-i) is required."
    show_help
    exit 1
fi

# Validate input file existence
if [[ ! -f "$INPUT_FILE" ]] || [[ ! -s "$INPUT_FILE" ]]; then
    log_error "Error: The input file '$INPUT_FILE' does not exist or is empty."
    exit 1
fi

# Validate outgroup file if provided
if [[ -n "$OUTGROUP_FILE" ]] && ([[ ! -f "$OUTGROUP_FILE" ]] || [[ ! -s "$OUTGROUP_FILE" ]]); then
    log_error "Error: The outgroup file '$OUTGROUP_FILE' is missing or empty."
    exit 1
fi

if [[ ! -d "$READS_DIR" || -z "$(ls -A "$READS_DIR" 2>/dev/null)" ]]; then
    log_error "Error: The reads directory '$READS_DIR' does not exist or is empty."
    exit 1
fi

# Extract the header line from the input file
HEADER_LINE="$(head -n1 "$INPUT_FILE" | tr -d '[:space:]')"
IFS=',' read -ra HEADER_COLS <<< "$HEADER_LINE"
LOWER_HEADER=("${HEADER_COLS[@],,}")
HAS_CODE_COLUMN=false
# Validate the number of columns first
if [[ "${#HEADER_COLS[@]}" -lt 2 ]]; then
  log_error "Error: The header must have at least 2 columns: TAXON/SPECIES/STRAIN and ACCESSIONS (or equivalent)."
  exit 1
elif [[ "${#HEADER_COLS[@]}" -gt 3 ]]; then
  log_error "Error: The header has more than 3 columns. This is not supported."
  exit 1
fi
# Check specific column structure
if [[ "${#HEADER_COLS[@]}" -eq 3 ]]; then
  # If has exactly 3 columns, check if the second column is CODE
  if [[ "${LOWER_HEADER[1]}" =~ ^code(s)?$ ]]; then
    HAS_CODE_COLUMN=true
    line_number=1
    while IFS= read -r line || [[ -n $line ]]; do
      ((line_number+=1)) #How do I sum from
      if ! [[ "$line" =~ ^[[:alnum:]]{5}$ ]]; then
        log_error "Invalid CODE on line $line_number:'$line'. CODE must have exactly 5 alphanumeric characters."
        exit 1
      fi
    done <<< "$(cat "orthomyxoviridae_accessions.txt" | tail -n +2 |grep -v '^$'|cut -d',' -f2|tr -d '[:blank:]')"
    #first I deleted complete void lines, so if -z works, it is detecting a row that has its second field void
    #echo -e "species\tcode" >> $FIVE_LETTER_FILE
    log_info "CODE column has been detected and validated. (3 columns)."
  else
    log_error "The header has 3 columns, but the second is not a CODE column (Found: '${HEADER_COLS[1]}')."
    exit 1
  fi
elif [[ "${#HEADER_COLS[@]}" -eq 2 ]]; then
  # If has exactly 2 columns, check if they are valid (STRAIN/SPECIES and ACCESSIONS)
  if [[ "${LOWER_HEADER[0]}" =~ ^(taxon(s)?|strain(s)?|species)$ ]] && [[ "${LOWER_HEADER[1]}" =~ ^accession(s)?$ ]]; then
    log_info "Only 2 columns found in header: TAXON/SPECIES/STRAIN, ACCESSIONS (or equivalent). No CODE column."
  else
    log_error "The header columns are not correct. Expected TAXON(S)/STRAIN/SPECIE(S) and ACCESSION(S) (or equivalent). Found: '${HEADER_COLS[0]}' and '${HEADER_COLS[1]}'."
    exit 1
  fi
fi
mkdir -p db

log_info "Starting retrieval of sequences from $INPUT_FILE..."
# Export the function so GNU Parallel can see it
count=0
export -f fetch_data log_info log_warn log_error
export RED YELLOW BLUE GREEN NC FIVE_LETTER_FILE HAS_CODE_COLUMN count
 # <-- Aquí inicializas la variable
tail -n +2 "$INPUT_FILE" | grep -v '^$' | parallel -j 1 fetch_data {}
#while IFS= read -r line || [[ -n "$line" ]]; do
#  if ! fetch_data "$line"; then
#    echo "Error processing line: $line"
#  fi
#done < <(tail -n +2 "$INPUT_FILE" | grep -v '^$')
#tail -n +2 "$INPUT_FILE" | grep -v '^$' | xargs -n 1 -I {} bash -c 'fetch_data "{}"'

#echo "Total lines processed: $count"

log_info "Finished retrieval of nucleotide sequences."

#Adding 5 letter code and cleaning for r2t
if $HAS_CODE_COLUMN; then
  # If the user gave a CODE column, we pass db_codes.tsv as second argument
  python "${SCRIPTS_DIR}/clean_fasta_cdna_cds2.py" db "$FIVE_LETTER_FILE"
  echo "HEREEEEE"
else
  # If no code was provided, the Python script generates codes on its own
  python "${SCRIPTS_DIR}/clean_fasta_cdna_cds2.py" db
  echo "PIPIPI"
fi

log_info "Editing parameters file..."
"${OMA}/oma" -p 

#Editing parameters file
log_info "Uncommenting last four steps of the parameters file..."
sed -i '/#WriteOutput_\(Phy\|Par\|H\)/ s/^#//' parameters.drw

outgroup_codes=()
if [ -n "$OUTGROUP_FILE" ]; then
    log_info "Reading outgroup taxon from $OUTGROUP_FILE..."
    while IFS= read -r species || [[ -n "$species" ]]; do
        species=$(echo "$species" | tr -d '[:space:]')
        if grep -q "^${species}[[:space:]]" "$FIVE_LETTER_FILE"; then
            code=$(awk -v sp="$species" -vOFS="_" '$1 == sp {print $1}' "$FIVE_LETTER_FILE")
            outgroup_codes+=("$code")
        else
            log_warn "Outgroup taxon '$species' not found in $FIVE_LETTER_FILE."
        fi
    done < "$OUTGROUP_FILE"
fi

#Verify codes and edit the parameters file
if [ ${#outgroup_codes[@]} -gt 0 ]; then
    outgroup_list=$(IFS=,; echo "${outgroup_codes[*]}")
    OUTGROUPS="OutgroupSpecies := [${outgroup_list}];"
    log_info "Using the following 5 letter code outgroup(s) to edit the parameters file: ${outgroup_list}."
else
    log_warn "No valid outgroup taxon provided. Using default 'none'."
    OUTGROUPS="OutgroupSpecies := 'none';"
fi
#map initial file with the new 5 letter code given by the clean...py
grep -q "^OutgroupSpecies" "$PARAMETERS_FILE" && \
    sed -i "s/^OutgroupSpecies.*/$OUTGROUPS/" "$PARAMETERS_FILE"

log_info "Running OMA..."
"${OMA}/oma" -n ${THREADS}
echo "Status"
#oma-status
#echo "End status"
if ls Output/OrthologousGroupsFasta/*.fa >/dev/null 2>&1; then
    ###Create tsv
    generate_og_gene_tsv db Output/OrthologousGroupsFasta OG_genes.tsv
    mkdir -p marker_genes
    #####cat Output/OrthologousGroupsFasta/*.fa > dna_ref.fa
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