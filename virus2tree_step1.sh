#!/bin/bash
set -euo pipefail
#set -x
PROGNAME="$(basename "$0")"
MAIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${MAIN_DIR}"
WORK_DIR=""
#OMA="${MAIN_DIR}/../oma/bin"
PARAMETERS_FILE="parameters.drw"
FIVE_LETTER_FILE="five_letter_taxon.tsv"
INPUT_FILE=""
OUTGROUP_FILE=""
THREADS=12
TEMP_DIR=""
OUT_DIR=""
DEBUG=false
MAT_PEPTIDES=false
ONLY_MAT_PEPTIDES=false
RES_DOWN=false
RES_DOWN_VOID=false
NCBI_DOWNLOAD_COUNT=0
P_FLAG=false
Q_FLAG=false
SKIP_STEP4=false

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
# rest
log_info() {
  echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]${NC} $*"
}

log_warn() {
  echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]${NC} $*"
}

log_error() {
  echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR]${NC} $*" >&2
}

check_dependencies() {
  #Initialize array to map tools to messages
  local missing=0
  declare -A tools=(
    ["read2tree"]="read2tree"
    ["oma"]="OMA standalone"
    ["efetch"]="Entrez Direct utilities"
  )

  #log_info "Checking system dependencies..."
  #Looping through the array to test for commands and detect possible missing tools  
  for cmd in "${!tools[@]}"; do
    if ! command -v "$cmd" &>/dev/null; then
      log_error "Missing requirement: ${tools[$cmd]}"
      ((missing++))
    #else
      #log_info "Found: ${tools[$cmd]}"
    fi
  done
}

usage() {
  log_info "Usage: ${PROGNAME} -i <accession_file> [-g <outgroup_file>] [options]\n"
  #echo "Try '$0 --help' for more information."
}

show_help() {
  cat << EOF
$(usage)

Required:
  -i, --input <file>                           Path to the accession file. Should be a comma-separated file with the structure: taxon(s)/species/strain(s),code(optional),accession(s). If provided, the code should have exactly 5 alphanumeric characters.

Optional:
  -g, --outgroup <file>                        Path to the outgroup taxon/species/strain file
  -p, --use_mat_peptides                       Downloads the gbk file for each taxon's accession(s). If at least one mature peptide feature is detected, these features are used as the coding sequences; otherwise, the standard CDS features are downloaded.
  -q, --use_only_mat_peptides                  Downloads the gbk file for each taxon's accession(s). If at least one mature peptide feature is detected, these features are used as the coding sequences; if none are detected, that taxon is skipped.
  -T, --threads <int>                          Number of threads [default: 12]
  --root_dir <dir>                             Specify root directory where all the outputs will be saved [default: current directory]
  --temp_dir <dir>                             Specify temporary directory (otherwise mktemp -d is used). If relative, it will be relative to the root_dir.
  --out_dir <dir>                              Specify output directory for read2tree step1 [default: read2tree_output]. If relative, it will be relative to the root_dir.
  --resume_download                            Skips taxa whose coding sequences have already been downloaded from NCBI to the db folder                                
  --debug                                      Keeps temporary directory
  -h, --help                                   Show this help message

Example:
  $PROGNAME -i accessions.txt -g outgroups.txt

EOF
exit 0
}

clean_line() {
  local input="$1"
  echo "$input" | tr -cd '[:alnum:]'
}

skip_taxa() {
    local CLEAN_FILE="$1"
    local filtered_input_file="${CLEAN_FILE%.*}_filtered.csv"
    local taxon_list=""
    #echo "This is the filtered_input_file: ${filtered_input_file}"
    # Create an associative array of downloaded taxa already in db dir
    declare -A downloaded_taxa_map
    while IFS= read -r file; do
        taxon=$(basename "$file" | awk -F '_cds_' '{print $1}' | tr -cd '[:alnum:]')
        downloaded_taxa_map["$taxon"]=1
    done < <(find db -maxdepth 1 -type f -name "*_cds_from_genomic.fna")
    {
        head -n1 "$CLEAN_FILE" #iterate over each taxon to see if it was already downloaded (if it matches at least one name from db folder)
        while IFS= read -r line || [[ -n "$line" ]]; do
            # Clean taxa
            raw_taxon=$(cut -d',' -f1 <<<"$line" | tr -cd '[:alnum:]')
            if [[ -n "${downloaded_taxa_map[$raw_taxon]:-}" ]]; then
              if [ -s "db/${raw_taxon}_cds_from_genomic.fna" ]; then
                taxon_list+="${raw_taxon}_cds_from_genomic.fna, "
              else
                rm -f "db/${raw_taxon}_cds_from_genomic.fna"
                log_warn "Removed invalid file: ${raw_taxon}_cds_from_genomic.fna. Will re-download" >&2
                echo "$line"
              fi
              unset 'downloaded_taxa_map["$raw_taxon"]'
            else
              #if it wasnt found, then the line is conserved in the filtered file
              echo "$line"
            fi
        done < <(tail -n +2 "$CLEAN_FILE")
    } > "$filtered_input_file"

    if [ ${#downloaded_taxa_map[@]} -gt 0 ]; then
        local extra_files=()
        for taxon in "${!downloaded_taxa_map[@]}"; do
            extra_files+=( "${taxon}_cds_from_genomic.fna" )
        done
        log_error "Found file(s) in the 'db/' folder that are not present in the accession file: ${extra_files[*]}. Please delete them"
        return 1
    fi

    if awk 'NR >=2 && NF {found=1; exit} END {exit !found}' "$filtered_input_file"; then #The NF verifies if there are fields, if so, found=1, and the exit will give 0, which is success
        log_info "Skipping the download of the following detected input files: ${taxon_list%, }"
    else
        log_info "All taxa had a corresponding file, skipping downloading from NCBI"
        RES_DOWN_VOID=true
        #log_info "Look at the value of RESDOWN VOID = ${RES_DOWN_VOID}"
        #Here somehow skip all the downloading section, use a flag?
    fi

    cp "$filtered_input_file" "$CLEAN_FILE"
    rm -f "$filtered_input_file"
}

fetch_data() {
  #error or warn and skip? better error!!
  local line="$1"
  # if [[ "$line" =~ ^# ]]; then
  #   log_warn "Skipping commented line: $line"
  #   return 0
  # fi
  #delete all spaces
  IFS=',' read -ra columns <<< "$(echo "$line" | tr -d '[:space:]')"
  local strain=$(clean_line "${columns[0]}") #Here I'm revoing everything that is not alnum, so it is not necessary to do this in the funcion clean..py
  if [[ -z "$strain" ]]; then
    log_error "Line with empty taxon: ${line}"
    return 1
  #else
    #log_info "Processing line for taxon ${strain}"
  fi
  # Decide if we have a code column or not
  local code=""
  local accessions_list=""
  if $HAS_CODE_COLUMN; then
    # The second column must be CODE
    if [[ "${#columns[@]}" -lt 3 ]]; then
      log_error "Line has fewer than 3 columns but we expected accessions. Line: $line"
      #this exit only gets me out of the subprocess, I don't think so, only of one instance of parallel?. The point is that the rest continue to run
      #Pending see what happens when parallel is remove and we use while instead
      return 1
    else
      code="${columns[1]}"
      accessions_list=$(IFS=','; echo "${columns[@]:2}")
    fi
  else
    # No code column: second column onward => accessions
    if [[ "${#columns[@]}" -lt 2 ]]; then
      log_error "Line has fewer than 2 columns but we expected taxa and accessions. Line: $line"
      return 1
    fi
    accessions_list=$(IFS=','; echo "${columns[@]:1}")
  fi

  if [[ -z "$accessions_list" ]]; then
    log_error "Taxon with empty accessions: ${strain}"
    return 1
  fi

  if [[ "$accessions_list" == *GCF_* || "$accessions_list" == *GCA_* ]]; then
    local assembly_accessions=()
    local regular_accessions=()
    #From list to array
    IFS=',' read -ra accessions <<< "$accessions_list"
    for acc in "${accessions[@]}"; do
      [[ $acc == GCF_* || $acc == GCA_* ]] && assembly_accessions+=("$acc") || regular_accessions+=("$acc")
    done
    log_info "Analyzing taxon ${strain}. Detected the following GCF/GCA assemblies: ${assembly_accessions[*]}"
    local nc_accessions=()
    for assembly in "${assembly_accessions[@]}"; do
	    log_info "Retrieving nucleotide accession(s) associated with assembly ${assembly}"
      # Input all the accessions to the array all_accs
      mapfile -t all_accs < <(esearch -db assembly -query "$assembly" \
        | elink -target nuccore \
        | efetch -format acc ) 

      # Output error if no accessions were found
      if [[ ${#all_accs[@]} -eq 0 ]]; then
        log_error "No nucleotide accessions found for assembly: $assembly"
        return 1
      else
        # Let' see if at least one accession begins with NC_ (RefSeq)
        mapfile -t nc_only < <(printf '%s\n' "${all_accs[@]}" | grep '^NC_')
        if [[ ${#nc_only[@]} -gt 0 ]]; then
		log_info "Retrieved ${#nc_only[@]} RefSeq accession(s) for assembly $assembly: ${nc_only[*]}."
          nc_accessions+=("${nc_only[@]}")
        else
		log_info "No RefSeq (NC_*) accessions found for assembly: $assembly. Falling back to ${#all_accs[@]} GenBank accession(s): ${all_accs[*]}."
          # Add all accessions (GenBank)
          nc_accessions+=("${all_accs[@]}")
        fi
      fi
    done
    accessions=("${nc_accessions[@]}" "${regular_accessions[@]}")
    # #To avoid problem with commas if the regular array is void?
    # if [[ -z "${regular_accessions[*]}" ]]; then
    #   #echo "There are no regular accessions: ${regular_accessions[@]}"
    #   accessions_list=$(IFS=','; echo "${nc_accessions[*]}")
    # else
    #   #echo "There are regular accessions: ${regular_accessions[@]}"
    #   accessions_list=$(IFS=','; echo "${nc_accessions[*]},${regular_accessions[*]}")
    # fi
    #echo "Final accesion list has ${#regular_accessions[@]} accession(s) and ${#nc_accessions[@]} assembly accession(s). The assembly accession(s) are: ${nc_accessions[@]}"
    local accessions_list=$(IFS=','; echo "${accessions[*]}")
  fi
  # Fetch data
  
  log_info "Fetching data for taxon: ${strain}. Final nucleotide accession(s) used: ${accessions_list//,/ }."
  
  if $MAT_PEPTIDES; then
    #downloads all the ids in a single gbk file
    log_info "--use_mat_peptides parameter specified, searching for gbk file and evaluating..."
    
    efetch -db nucleotide -id "$accessions_list" -format gbwithparts > "${TEMP_DIR}/gbk_dir/${strain}.gbk"
    if python "${SCRIPTS_DIR}/write_mat_peptides.py" "${TEMP_DIR}/gbk_dir/${strain}.gbk" "db/${strain}_cds_from_genomic.fna"; then
      if [[ -f "db/${strain}_cds_from_genomic.fna" ]]; then
        if [[ -s "db/${strain}_cds_from_genomic.fna" ]]; then
          if $HAS_CODE_COLUMN; then
            log_info "Writing 5-letter code for taxon ${strain} to ${FIVE_LETTER_FILE}"
            echo -e "${strain}\t${code}" >> "${FIVE_LETTER_FILE}"
          fi
          ((NCBI_DOWNLOAD_COUNT++))
          return 0
        else
          log_error "Writing of mat_peptide features to db/${strain}_cds_from_genomic.fna failed: file is empty for taxon ${strain}: ${accessions_list}"
          return 1
        fi

      elif [[ $ONLY_MAT_PEPTIDES == true ]]; then
	      log_info "Skipping taxon ${strain}..."
	      return 0
      else
        log_info "Fetching CDS features..."
      fi
      
    else
      log_error "Processing of '$( realpath "${TEMP_DIR}/gbk_dir/${strain}.gbk")' file failed"
      return 1
    fi
  fi

  efetch -db nucleotide -id "$accessions_list" -format fasta_cds_na \
    > "db/${strain}_cds_from_genomic.fna" \
    || {
      log_warn "Command efetch failed while fetching accession(s) for taxon ${strain}: ${accessions_list}"
      rm -f "db/${strain}_cds_from_genomic.fna" #Ensures this file won't exist
      return 0
    }

  if [ ! -s "db/${strain}_cds_from_genomic.fna" ]; then
	  log_warn "Command efetch failed to fetch accession(s) for taxon ${strain}: ${accessions_list}"
	  rm -f "db/${strain}_cds_from_genomic.fna"
    return 0
  fi


  ((NCBI_DOWNLOAD_COUNT++))

  #write to the file only if fetching was successful
  if $HAS_CODE_COLUMN; then
    log_info "Writing 5-letter code for taxon ${strain} to ${FIVE_LETTER_FILE}"
    echo -e "${strain}\t${code}" >> "${FIVE_LETTER_FILE}"
  fi
  sleep 1
  log_info "Successfully fetched data for taxon ${strain}\n"
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
    log_info "Loading 5-letter code for each taxon from file ${FIVE_LETTER_FILE}"
    while IFS=$'\t' read -r species_val code_val; do
          SPECIES_TO_CODE["$species_val"]="$code_val"
    done < "$FIVE_LETTER_FILE"
    # Example usage
    # process_genes ~/oma/test3_illumina/db ~/oma/test3_illumina/marker_genes output_table.tsv
    # process_genes ~/oma/test3_illumina/db ~/oma/test3_illumina/marker_genes output_table.tsv
    tmp_file="$TEMP_DIR/OG_genes_unsorted.tsv"
    output_file2="$TEMP_DIR/${output_file}.tmp"
    unique_output_file2="$TEMP_DIR/${unique_output_file}.tmp"
    log_info "Processing gene features from coding sequences files in directory: $fna_dir"
    while IFS=: read -r file line; do
        #Example of the complete input line: db/rsv_11_cds_from_genomic.fna:>lcl|MG813984.1_cds_AZQ19553.1_6 [gene=SH] [protein=small hydrophobic protein] [protein_id=AZQ19553.1] [location=4251..4445] [gbkey=CDS]
        local protein_id="NA"
        # protein_id
        if [[ "$line" =~ \[protein_id=([^]]+)\] ]]; then
          protein_id="${BASH_REMATCH[1]}"
        fi
        if [[ -z "$protein_id" ]] || [[ "$protein_id" == "NA" ]]; then
          log_error "Error: Missing or invalid protein_id in file: $file, line: $line"
          return 1
        fi

        local compressed_id="$(clean_line "$protein_id")"
        # Search OG using grep
        local og_matches
        #Here it can match partial substrings
        if og_matches=$(grep -l "$compressed_id" "$fa_dir"/*.fa 2>/dev/null | xargs -I {} basename {} .fa); then
              log_info "Success: Match found for protein ID $protein_id in OG $og_matches"
              # Iterate over the OGs found and write to temporary file
              #This loop assumes a gene could be present in more than one OG, is that actually the case? don't think so
              local species=$(basename "$file" | awk -F '_cds_' '{print $1}')
              local code="${SPECIES_TO_CODE[$species]}"
              local accession="$(awk -F '_cds_' '{print $1}' <<< "${line#*|}")"
              local gene="NA"
              local protein="NA"
              local location="NA"
              local locus_tag="NA"
              local db_xref="NA"
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

              # location
              if [[ "$line" =~ \[location=([^]]+)\] ]]; then
                location="${BASH_REMATCH[1]}"
              fi
              #Locus tag
              if [[ "$line" =~ \[locus_tag=([^]]+)\] ]]; then
                locus_tag="${BASH_REMATCH[1]}"
              fi
              #Db_xref
              if [[ "$line" =~ \[db_xref=([^]]+)\] ]]; then
                db_xref="${BASH_REMATCH[1]}"
              fi
              #Compressed_id uses the protein_id without alphanumeric characters to match the header in OrthologousGroupsFasta:
              #Output/OrthologousGroupsFasta/OG11.fa:>3355X||rsv_11||lcl|MG8139841cdsAZQ1955316 cleaned for r2t [rsv_11_3355X]
              #local compressed_id="$(echo "$protein_id" | sed 's/[^a-zA-Z0-9]//g')"
              #echo "Got correct features from *fna for ${protein_id}"
              [[ "$gene" == "NA" ]]       && log_info "Missing gene for protein ID $protein_id in $file"
              [[ "$protein" == "NA" ]]    && log_info "Missing product(protein) for protein ID $protein_id in $file"
              [[ "$location" == "NA" ]]  && log_info "Missing location for protein ID $protein_id in $file"
              [[ "$locus_tag" == "NA" ]]  && log_info "Missing locus_tag for protein ID $protein_id in $file"
              [[ "$db_xref" == "NA" ]]  && log_info "Missing db_xref for protein ID $protein_id in $file"

              for og in $og_matches; do
                  echo -e "${og}\t${gene}\t${protein}\t${protein_id}\t${location}\t${accession}\t${species}\t${code}\t${locus_tag}\t${db_xref}" >> "$tmp_file"
              done
        else #Should this be really a warn?, it just means that that gene is not in the OG from OMA
            log_info "No match found for protein ID '$protein_id' in OGs directory. This sequence is not part of any OG from OMA."
        fi
    done < <(grep '^>' "$fna_dir"/*.fna)
    log_info "Finished processing all coding sequences files in directory: $fna_dir"
    #-V option in sort from GNU coreutils
    sort -k1,1 -V -k2,2 "$tmp_file" > "$output_file2"
    sed -i '1iOG\tGene\tProtein\tProtein_ID\tLocation\tAccession\tTaxon\tCode\tLocus_tag\tDb_xref' "$output_file2"
    {
      head -n +1 "$output_file2"
      # Skip the header of the main TSV
      tail -n +2 "$output_file2" \
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
    } > "$unique_output_file2"

    for file in "$output_file2" "$unique_output_file2"; do
      final_file_base="${file##*/}"; final_file="${final_file_base%.*}"
      awk 'BEGIN { FS = OFS = "\t" }
      NR==1 {
          n = NF
          for (i = 1; i <= NF; i++) {
             header[i] = $i
          }
          next
      }
      {
          for (i = 1; i <= NF; i++) {
             if ($i != "NA") {
                nonNA[i] = 1
             }
          }
          data[NR] = $0
      }
      END {
      #Now remades the header suing directly the nonNA array
          out = ""
          for (i = 1; i <= n; i++) {
             if (nonNA[i])
                out = out header[i] OFS
          }
          sub(OFS "$", "", out)
          print out
	  #Now writes each row of data using nonNA array
          for (j = 2; j <= NR; j++) {
             split(data[j], fields, FS)
             out_line = ""
             for (i = 1; i <= n; i++) {
                if (nonNA[i])
                   out_line = out_line fields[i] OFS
             }
             sub(OFS "$", "", out_line)
             print out_line
          }
      }' "$file" > "$final_file"
    done
    log_info "OG-Gene TSV generation complete: $output_file"
}

####################################################

#MAIN

###################################################

if [[ $# -eq 0 ]]; then
  usage
  log_info "Try '$PROGNAME --help' for more information."
  exit 1
fi

log_info "Script invoked with: $PROGNAME $*\n"

log_info "========== Step 1.1: Validating parameters =========="

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_FILE="$2"; shift ;;
        -g|--outgroup) OUTGROUP_FILE="$2"; shift ;;
        -p|--use_mat_peptides) MAT_PEPTIDES=true; P_FLAG=true;;
        -q|--use_only_mat_peptides) MAT_PEPTIDES=true; ONLY_MAT_PEPTIDES=true; Q_FLAG=true;;
        -T|--threads) THREADS="$2"; shift ;;
        --temp_dir) TEMP_DIR="${2%/}"; shift ;;
        --root_dir) WORK_DIR="${2%/}"; shift ;;   # NUEVO
        --out_dir) OUT_DIR="${2%/}"; shift ;;
        --debug) DEBUG=true;;
        --resume_download) RES_DOWN=true;;
        -h|--help) show_help;;
        *) log_error "Unknown parameter passed: $1"; usage; log_info "Try '$PROGNAME --help' for more information."; exit 1 ;;
    esac
    shift
done

# Validate dependencies at the very beggining
check_dependencies
log_info "Checked system dependencies"

if [[ -z "${INPUT_FILE}" ]]; then
    log_error "Error: --input (-i) is required."
    usage
    exit 1
fi

# Validate input file existence
if [[ ! -f "$INPUT_FILE" ]] || [[ ! -s "$INPUT_FILE" ]]; then
    log_error "The input file '$INPUT_FILE' does not exist or is empty."
    usage
    exit 1
fi

# Validate outgroup file if provided, yes -f tolerates sym links
if [[ -n "$OUTGROUP_FILE" ]] && { [[ ! -f "$OUTGROUP_FILE" ]] || [[ ! -s "$OUTGROUP_FILE" ]]; }; then
    log_error "Error: The outgroup file '$OUTGROUP_FILE' is missing or empty."
    exit 1
fi

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
  log_error "Error: THREADS must be a positive integer."
  exit 1
fi

if [ "$P_FLAG" = true ] && [ "$Q_FLAG" = true ]; then
    log_error "Options -p and -q cannot be used together"
    exit 1
fi

if [[ -n "$WORK_DIR" ]]; then
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"
    log_info "Using work directory: '$(pwd)'"
else
    # por defecto: cwd del usuario al ejecutar el script
    WORK_DIR="$(pwd)"
    log_info "No --work_dir specified. Outputs will be written in '$(pwd)'"
fi

if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR=read2tree_output
  log_info "No name for the read2tree directory was specified, using '$(realpath "$OUT_DIR")'"
  if [[ -d "$OUT_DIR" ]]; then
	  log_error "The read2tree output directory '$(realpath "$OUT_DIR")' already exists. Please provide a novel read2tree directory"
    exit 1
  fi
else
  if [[ -d "$OUT_DIR" ]]; then
	  log_error "The read2tree output directory '$(realpath "$OUT_DIR")' already exists. Please provide a novel read2tree directory"
    exit 1
  fi
  log_info "Using output directory: '$(realpath "$OUT_DIR")'"
fi

if [[ "$RES_DOWN" == true ]]; then
    #This should be done at the very beggining
    if [[ ! -d "db" ]]; then
        log_error "Directory 'db/' does not exist. Are you sure about resuming a previous download?"
        exit 1
    fi
    if [[ -d "$TEMP_DIR" ]]; then
        if [[ $(ls "$TEMP_DIR" 2>/dev/null) ]]; then
            log_error "Temporary directory '$(realpath "$TEMP_DIR")' already exists and is not void. Please specify a novel directory to avoid conflicts when resuming the download."
            exit 1
        fi
    fi
else
  if [[ -d "db" ]]; then
        log_error "Directory $(realpath "db/")' already exists, please remove it"
        exit 1
  fi
  if [[ -d "DB" ]]; then
        log_error "Directory $(realpath "DB/")' already exists, please remove it"
        exit 1
  fi
  if [[ -d "$TEMP_DIR" ]]; then
        if [[ $(ls "$TEMP_DIR" 2>/dev/null) ]]; then
            log_error "Temporary directory '$(realpath "$TEMP_DIR")' already exists and is not void. Please specify a novel directory to avoid conflicts with the new run."
            exit 1
        fi
  fi
  mkdir -p db
fi

if [[ -z "$TEMP_DIR" ]]; then
  TEMP_DIR="$(mktemp -d)"
  log_info "Created temp directory at '$TEMP_DIR'"
else
  # Validate if it is a directory
  if [[ ! -d "$TEMP_DIR" ]]; then
    mkdir -p "$TEMP_DIR"
    log_info "Using temp directory: '$(realpath "$TEMP_DIR")'"
  else 
    if [[ $(ls "$TEMP_DIR" 2>/dev/null) ]]; then
        log_error "Temporary directory '$(realpath "$TEMP_DIR")' already exists and is not void. Please specify a novel directory to avoid conflicts with the new run."
        exit 1
    fi
  fi
fi

if [[ "$MAT_PEPTIDES" == true ]]; then
  mkdir -p "$TEMP_DIR"/gbk_dir
fi

if [[ "$DEBUG" == false ]]; then
  trap '[[ -n "$TEMP_DIR" && -d "$TEMP_DIR" ]] && rm -rf "./$TEMP_DIR"' EXIT
else
  log_info "Debug mode enabled, keeping temporary directory: '$(realpath "$TEMP_DIR")'"
fi

log_info "========== Step 1.2: Validating input file =========="

# Export the functions and variables
export -f skip_taxa fetch_data log_info log_warn log_error clean_line 
export RED YELLOW GREEN NC FIVE_LETTER_FILE HAS_CODE_COLUMN NCBI_DOWNLOAD_COUNT RES_DOWN_VOID RES_DOWN ONLY_MAT_PEPTIDES MAT_PEPTIDES

#CLEAN_FILE="$(mktemp -p "$TEMP_DIR" file.XXXXXX)"  #Temporary file for cleaned input
CLEAN_FILE="$TEMP_DIR/input_clean_file.txt"
# Extract the header line from the input file
log_info "Cleaning accession file and saving to '$(realpath "$CLEAN_FILE")'"
grep -v '^#' "$INPUT_FILE" | grep -v '^$' > "$CLEAN_FILE"

HEADER_LINE="$(head -n1 "$CLEAN_FILE" | tr -d '[:space:]')"
IFS=',' read -ra HEADER_COLS <<< "$HEADER_LINE"
LOWER_HEADER=("${HEADER_COLS[@],,}")
HAS_CODE_COLUMN=false
# Validate the number of columns first
if [[ "${#HEADER_COLS[@]}" -lt 2 ]]; then
  log_error "Error: The header must have at least 2 columns: 'taxon/species/strain' and 'accession(s)' (or equivalent)."
  exit 1
elif [[ "${#HEADER_COLS[@]}" -gt 3 ]]; then
  log_error "Error: The header has more than 3 columns. This is not supported."
  exit 1
fi
# Check specific column structure
if [[ "${#HEADER_COLS[@]}" -eq 3 ]]; then
  if [[ "${LOWER_HEADER[0]}" =~ ^(taxon|taxa|strain(s)?|species)$ ]] && [[ "${LOWER_HEADER[2]}" =~ ^accession(s)?$ ]]; then
    log_info "Detected taxon and accession columns in the header. Now checking for the code column..." 
    # If has exactly 3 columns, check if the second column is CODE
    if [[ "${LOWER_HEADER[1]}" =~ ^code(s)?$ ]]; then
      HAS_CODE_COLUMN=true
      set +e
      awk_output=$(awk -F ',' '
      NR == 1 { next } 
      {
          codigo_col1 = $1
          original_col1 = codigo_col1  
          gsub(/[^[:alnum:]]/, "", codigo_col1)  

          codigo_col2 = $2
          gsub(/[[:blank:]]/, "", codigo_col2)

          if (codigo_col1 in visto_col1) {
              printf "Duplicated taxon in column 1: \047%s\047 (line %d). First occurrence: line %d.\nNote: Comparison is based on alphanumeric characters (ignoring symbols and spaces).\n",
                    original_col1, NR, visto_col1[codigo_col1]
              exit 1 
          }
          visto_col1[codigo_col1] = NR

          if (length(codigo_col2) != 5 || codigo_col2 !~ /^[[:alnum:]]+$/) {
              printf "Invalid 5-letter code on line %d: \047%s\047. The code must have 5 alphanumeric characters.\n",
                    NR, codigo_col2
              exit 1  
          }

          if (codigo_col2 in visto_col2) {
              printf "Duplicated 5-letter code on line %d: \047%s\047. First occurrence: line %d.\nNote: Comparison is based on alphanumeric characters (ignoring symbols and spaces).\n",
                    NR, codigo_col2, visto_col2[codigo_col2]
              exit 1  
          }
          visto_col2[codigo_col2] = NR
      }' "$CLEAN_FILE" 2>&1)
      if [[ -n "$awk_output" ]]; then
          log_error "$awk_output"
          exit 1
      fi
      set -e
      #first I deleted complete void lines, so if -z works, it is detecting a row that has its second field void
      log_info "The code column has been detected and validated."
      if [[ "$RES_DOWN" == true ]]; then
        skip_taxa "$CLEAN_FILE"
      fi
    else
      log_error "The header has 3 columns, but the second is not a 'code(s)' column (Found: '${HEADER_COLS[1]}')."
      exit 1
    fi
  else
    log_error "The header columns are not correct. Expected 'taxon/species/strain' and 'accession(s)' (or equivalent). Found: '${HEADER_COLS[0]}' and '${HEADER_COLS[2]}'."
    exit 1
  fi
elif [[ "${#HEADER_COLS[@]}" -eq 2 ]]; then
  # If has exactly 2 columns, check if they are valid (STRAIN/SPECIES and ACCESSIONS)
  if [[ "${LOWER_HEADER[0]}" =~ ^(taxon|taxa|strain(s)?|species)$ ]] && [[ "${LOWER_HEADER[1]}" =~ ^accession(s)?$ ]]; then
    log_info "Detected only taxon and accession columns in the header. No code column found. Unique codes for each taxon in the format 'sXXXX' will be automatically generated."
    # Validate duplicates in column 1
    set +e
    awk_output=$(awk -F ',' '
    NR == 1 { next }  
    {
        taxon = $1
        original_taxon = taxon  # keep original for the log message
        gsub(/[^[:alnum:]]/, "", taxon)  # only alnum characters are compared

        # check for duplicatres
        if (taxon in visto) {
            printf "Duplicated taxon: \047%s\047 (line %d). First occurrence: line %d.\nNote: Comparison is based on alphanumeric characters (ignoring symbols and spaces).\n",
                   original_taxon, NR, visto[taxon]
            exit 1  
        }
        visto[taxon] = NR
    }' "$CLEAN_FILE" 2>&1)
    # Manejar errores
    if [[ -n "$awk_output" ]]; then
        log_error "$awk_output"
        exit 1
    fi
    set -e
    if [[ "$RES_DOWN" == true ]]; then
      skip_taxa "$CLEAN_FILE"
    fi

  else
    log_error "The header columns are not correct. Expected 'taxon/species/strain' and 'accession(s)' (or equivalent). Found: '${HEADER_COLS[0]}' and '${HEADER_COLS[1]}'."
    exit 1
  fi
fi
if [[ "$RES_DOWN_VOID" == false ]]; then
  log_info "========== Step 1.3: Retrieving coding sequences from NCBI =========="

  log_info "Starting retrieval of sequences from $INPUT_FILE..."

  count=0
  #count has nothing to do, exit status are always cero, the IFS is end of line and the $ are where they should be,
  #the only clue is that this always happens in GCF/GCA lines of the clean file and never with parallel!! this did not happen with parallel
  #it is not sth of the while loop, as I ran it apart and ran well, its sth of the fetch_data function!!, but I dont know what 
  #because the exit code is always cero!!!, and when it stops it simply gets out the while loop WITHOUT passing by my error message!!!
  #And after 3 hours of trials, the solution was to add <dev/null to prevent the function (or the commands it invokes) from consuming the standard script input.
  while IFS= read -r line || [[ -n $line ]]; do
    log_info "Processing line: $line"
    if ! fetch_data "$line" < /dev/null; then
      log_error "Failed to fetch data for line: $line"
      exit 1
    fi
    ((count+=1))
  done <<< "$(tail -n +2 "$CLEAN_FILE")"

  log_info "Processed $count lines from the accession file"
  log_info "Downloaded accession(s) from $NCBI_DOWNLOAD_COUNT taxa" 
  log_info "Nucleotide sequences retrieval completed successfully.\n"
fi

if [[ "$RES_DOWN_VOID" == false && "$NCBI_DOWNLOAD_COUNT" -eq 0 && "$RES_DOWN" == true ]]; then
  RES_DOWN_VOID=true
fi

if [ "$RES_DOWN_VOID" == true ] && [ -d "DB" ] && [ -s "dna_ref.fa" ]; then
  SKIP_STEP4=true
  for DB_file in DB/*.fa; do
      base_name=$(basename "$DB_file" .fa)
      target_file="db/${base_name}_cds_from_genomic.fna"
      if [ ! -s "$target_file" ]; then
          log_error "The file ${DB_file} doesn't have a corresponding ${target_file} in the db/ directory Please check your process."
          exit 1
      fi
  done
  # Now check sequence count of dna ref to be the same in all the seq on db
  if $SKIP_STEP4; then
      log_info "All the files in the DB folder have a corresponding one in the db folder"
      dna_ref_count=$(grep -c '^>' "dna_ref.fa" 2>/dev/null || echo 0)
      db_total_count=$(grep -ch '^>' db/*_cds_from_genomic.fna | awk '{sum += $1} END {print sum}')
      if [ "$dna_ref_count" -ne "$db_total_count" ] || [ "$dna_ref_count" -eq 0 ]; then
          SKIP_STEP4=false
          log_info "Sequences in dna_ref.fa ($dna_ref_count) do not match with those in the db folder ($db_total_count)"
      else
        log_info "The dna_ref.fa file has the same number of sequences as the db folder"
      fi
  fi
fi

if $SKIP_STEP4; then
  log_info "========== Skipping Step 1.4 : All files were already downloaded from NCBI and 'DB' folder and 'dna_ref.fa' file already exist =========="
else
  log_info "========== Step 1.4: Preparing format of coding sequences for OMA and read2tree =========="
  #Adding 5 letter code and cleaning for r2t
  if $HAS_CODE_COLUMN; then
    # If the user gave a CODE column, we pass db_codes.tsv as an argument
    python "${SCRIPTS_DIR}/clean_fasta_cdna_cds.py" db "$RES_DOWN" "$FIVE_LETTER_FILE" 
  else
    # If no code was provided, the Python script generates codes on its own
    python "${SCRIPTS_DIR}/clean_fasta_cdna_cds.py" db "$RES_DOWN"
  fi
fi

log_info "========== Step 1.5: Editing parameters.drw file for running OMA =========="

outgroup_codes=()
if [ -n "$OUTGROUP_FILE" ]; then
    log_info "Reading outgroup taxon(s) from $OUTGROUP_FILE..."
    while IFS= read -r species || [[ -n "$species" ]]; do
        species=$(clean_line "${species}")
        if grep -q "^${species}[[:space:]]" "$FIVE_LETTER_FILE"; then
            code=$(awk -v sp="$species" '$1 == sp {print $1}' "$FIVE_LETTER_FILE")
            outgroup_codes+=("$code")
        else
            log_error "Outgroup taxon '$species' not found in $FIVE_LETTER_FILE."
            exit 1
        fi
    done < "$OUTGROUP_FILE"
fi

if [ -f "parameters.drw" ]; then
  rm -f parameters.drw
  log_info "Existing parameters.drw file removed"
fi
log_info "Creating the parameters.drw file for OMA"
oma -p 
sed -i '/#WriteOutput_\(Phy\|Par\|H\)/ s/^#//' parameters.drw

#Verify codes and edit the parameters file
if [ ${#outgroup_codes[@]} -gt 0 ]; then
    outgroup_list=$(IFS=,; echo "${outgroup_codes[*]}")
    OUTGROUPS="OutgroupSpecies := [${outgroup_list}];"
    log_info "Using the following 5-letter code outgroup(s) to edit the parameters file: ${outgroup_list//,/ }"
else
    log_warn "No valid outgroup taxon provided. Using default 'none'."
    OUTGROUPS="OutgroupSpecies := 'none';"
fi
#map initial file with the new 5 letter code given by the clean...py
grep -q "^OutgroupSpecies" "$PARAMETERS_FILE" && \
    sed -i "s/^OutgroupSpecies.*/$OUTGROUPS/" "$PARAMETERS_FILE"

log_info "========== Step 1.6: Running OMA =========="

oma -n "${THREADS}"
#oma-status
#echo "End status"
if oma-status && ls Output/OrthologousGroupsFasta/*.fa >/dev/null 2>&1; then
    log_info "OMA finished successfully! OMA did great!"
else
    log_error "OMA has failed."
    exit 1
fi

log_info "========== Step 1.7: Gathering OG-Gene statistics =========="
#Now it does not make sense to check if Output/Orth...*.fa exists
log_info "Generating summary OG-gene TSV file"
###Create tsv
generate_og_gene_tsv db Output/OrthologousGroupsFasta OG_genes.tsv
mkdir -p marker_genes
#####cat Output/OrthologousGroupsFasta/*.fa > dna_ref.fa
mv Output/OrthologousGroupsFasta/*.fa marker_genes
log_info "========== Step 1.8: Running Read2tree (step 1 marker) =========="
log_info "Using ${THREADS} threads..."
#read2tree --standalone_path ./marker_genes --output_path read2tree_output --dna_reference dna_ref.fa
read2tree --step 1marker --standalone_path marker_genes --dna_reference dna_ref.fa --output_path "$OUT_DIR" --debug 
# echo "Starting read processing..."
# echo "Searching for FASTQ files in: $READS_DIR"
# all_fastq=( $(find "$READS_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) ) )
# echo "Found ${#all_fastq[@]} FASTQ files."

# # Declare associative arrays for both reads
# declare -A sampleR1
# declare -A sampleR2

# # First we define sampleR1 / sampleR2
# for f in "${all_fastq[@]}"; do
#     base="$(basename "$f")"
#     noext="${base%.fastq.gz}"
#     if [[ "$noext" == "$base" ]]; then
#         noext="${base%.fastq}"
#     fi

#     # If _1 is detected => Paired R1
#     if [[ "$noext" =~ (.*)_1$ ]]; then
#         sample="${BASH_REMATCH[1]}"
#         sampleR1["$sample"]="$f"
#         echo "Detected paired-end R1: $f (Sample: $sample)"
#         echo "Detected paired-end R2: $f (Sample: $sample)"

# If _2 is detected => Paired R2
#     elif [[ "$noext" =~ (.*)_2$ ]]; then
#         sample="${BASH_REMATCH[1]}"
#         sampleR2["$sample"]="$f"
#         echo "Detected paired-end R2: $f (Sample: $sample)"

#     else

#     else
#         # Otherwise => single-end
#         sampleR1["$noext"]="$f"
#         echo "Detected single-end read: $f (Sample: $noext)"
#     fi
# done

# # Get the union of all detected samples
# all_samples=()
# all_samples+=( "${!sampleR1[@]}" )
# all_samples+=( "${!sampleR2[@]}" )
# # We create a unique list by removing duplicates
# readarray -t unique_samples < <(printf '%s\n' "${all_samples[@]}" | sort -u)

# #main loop for processing each sample
# for sample in "${unique_samples[@]}"; do
#     R1="${sampleR1[$sample]:-}"
#     R2="${sampleR2[$sample]:-}"

#     if [[ -n "$R1" && -n "$R2" ]]; then
#         # => Paired-end => force short read type
#         echo "Processing short paired-end reads for sample '$sample':"
#         echo "  R1: $R1"
#         echo "  R2: $R2"
#         read2tree --step 2map \
#                   --standalone_path marker_genes \
#                   --dna_reference dna_ref.fa \
#                   --reads "$R1" "$R2" \
#                   --read_type "short" \
#                   --threads "$THREADS" \
#                   --output_path read2tree_output \
#                   --debug
#     else
#         # => Single-end => use READ_TYPE ('short', 'long-ont', 'long-hifi')
#         #    (if only R2 or only R1 exists, treat it as single-end)
#         SE_FILE="$R1"
#         if [[ -z "$SE_FILE" ]]; then
#             SE_FILE="$R2"
#         fi
#         echo "Processing single-end reads for sample '$sample' (READ_TYPE=$READ_TYPE):"
#         echo "  Read file: $SE_FILE"
#         read2tree --step 2map \
#                   --standalone_path marker_genes \
#                   --dna_reference dna_ref.fa \
#                   --reads "$SE_FILE" \
#                   --read_type "$READ_TYPE" \
#                   --threads "$THREADS" \
#                   --output_path read2tree_output \
#                   --debug
#     fi
# done

# echo "Read processing completed successfully."

#echo "Merging..."
#read2tree --step 3combine --standalone_path marker_genes --dna_reference dna_ref.fa --output_path read2tree_output --tree --debug
#

#iqtree -T ${THREADS} -s read2tree_output/concat_*_aa.phy -bb 1000
#iqtree -T ${THREADS} -s read2tree_output/concat_*_dna.phy 
