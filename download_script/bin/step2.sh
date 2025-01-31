#!/usr/bin/env bash
set -euo pipefail
set -x #bug fixes
echo -e "Script invoked with: $0 $*\n"

############################################
# Default values
############################################
READ_TYPE=""
READS=()         # Array to store multiple input reads
TEMP_DIR=""      # If empty, we will create a new one with mktemp
DEDUP=false
DEDUP_L=""       # optional argument for czid-dedup
DOWNSAMPLE=false
# Downsampling method can be coverage & genome_size, or number_of_bases, or number_of_reads
COVERAGE=0
GENOME_SIZE=0
NUM_BASES=0
NUM_READS=0
THREADS=4        # you can adjust a default or parse it
DEBUG=false

############################################
# Functions
############################################
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
  echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]${NC} $*"
}

log_error() {
  echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR]${NC} $*"
}

usage() {
  echo "Usage: $0 -r <file1> [file2...] -t <read_type> [options]"
  echo "Try '$0 --help' for more information."
}

show_help() {
  cat << EOF
$(usage)

  Required:
    -r, --reads <file1> [file2 file3 ...]            Input fastq/fastq.gz
    -t, --read_type <se_short|pe_short|pacbio|ont>    Type of reads


  Optional:
    -T, --threads <int>      Threads to use (default 4)
    --temp_dir <dir>         Specify existing temp directory (otherwise mktemp -d is used)
    --dedup                  Run czid-dedup
    --dedup_l <int>          Provide '-l <int>' to czid-dedup
    --downsample             Run rasusa
    --coverage <int>         Coverage to target (requires --genome_size)
    --genome_size <int>      Genome size in bases
    --num_bases <int>        Target number of bases
    --num_reads <int>        Target number of reads
    --debug                  Keeps temporary directory  
    -h, --help               Show this help
Examples:
  $0 -t se_short -r sample.fastq.gz --dedup --downsample --num_reads 10000
  $0 --read_type pe_short --reads R1.fastq R2.fastq --coverage 30 --genome_size 5000000
EOF
  exit 0
}

############################################
# Parse arguments
############################################

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

#Not positional arguments different from the options below
while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--read_type)
      READ_TYPE="$2"
      shift 2
      ;;
    -r|--reads)
      shift
      while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
        READS+=("$1")
        shift
      done
      ;;
    -T|--threads)
      THREADS="$2"
      shift 2
      ;;
    --temp_dir)
      TEMP_DIR="$2"
      shift 2
      ;;
    --dedup)
      DEDUP=true
      shift
      ;;
    --dedup_l)
      DEDUP_L="$2"
      shift 2
      ;;
    --downsample)
      DOWNSAMPLE=true
      shift
      ;;
    --coverage)
      COVERAGE="$2"
      shift 2
      ;;
    --genome_size)
      GENOME_SIZE="$2"
      shift 2
      ;;
    --num_bases)
      NUM_BASES="$2"
      shift 2
      ;;
    --num_reads)
      NUM_READS="$2"
      shift 2
      ;;
    --debug)
      DEBUG=true
      shift
      ;;
    -h|--help)
      show_help
      ;;
    *)
      log_error "Unknown parameter passed: $1"
      usage
      exit 1
      ;;
  esac
done
############################################
# Validate required parameters
############################################

if [[ -z "$READ_TYPE" ]]; then
  log_error "Please specify --read_type"
  usage
  exit 1
fi
if [[ ${#READS[@]} -eq 0 ]]; then
  log_error "Please specify at least one input file with --reads"
  usage
  exit 1
fi

case "$READ_TYPE" in
  se_short|pe_short|pacbio|ont) ;;
  *)
    log_error "Invalid read_type: $READ_TYPE. Must be one of se_short, pe_short, pacbio, ont."
    usage
    exit 1
    ;;
esac

if [[ "$READ_TYPE" == "pe_short" && ${#READS[@]} -ne 2 ]]; then
  log_error "Paired-end short reads requires exactly 2 files."
  usage
  exit 1
fi

#Verify the use of just one of the 3 rasusa options
if [[ "$DOWNSAMPLE" == true ]]; then
  method_count=0
  if [[ $COVERAGE -gt 0 && $GENOME_SIZE -gt 0 ]]; then
    ((method_count++))
  fi
  if [[ $NUM_BASES -gt 0 ]]; then
    ((method_count++))
  fi
  if [[ $NUM_READS -gt 0 ]]; then
    ((method_count++))
  fi
  if [[ $method_count -eq 0 ]]; then
    log_error "Downsampling requires --coverage + --genome_size OR --num_bases OR --num_reads."
    exit 1
  fi
  if [[ $method_count -gt 1 ]]; then
    log_error "Please specify only ONE downsampling method (coverage+genome_size) OR num_bases OR num_reads."
    exit 1
  fi
fi

############################################
# Create or verify temp directory
############################################
if [[ -z "$TEMP_DIR" ]]; then
  TEMP_DIR="$(mktemp -d)"
  log_info "Created temp directory at $TEMP_DIR"
else
  # Validate if it is a directory
  if [[ ! -d "$TEMP_DIR" ]]; then
    mkdir -p "$TEMP_DIR"
  fi
  log_info "Using existing temp directory: $TEMP_DIR"
fi

if [[ "$DEBUG" == false ]]; then
  trap 'rm -rf "$TEMP_DIR"' EXIT
else
  log_info "Debug mode enabled, keeping temporary directory: $TEMP_DIR"
fi

############################################
# Decompression and Concatenation if necessary
############################################
FINAL_READS=()
if [[ ${#READS[@]} -gt 1 ]]; then
  if [[ ${READ_TYPE} == "pe_short" ]]; then
    for READ_SAMPLE in "${READS[@]}"; do
      output_file=$TEMP_DIR/${READ_SAMPLE%.gz}
      if gzip -t "$READ_SAMPLE" 2>/dev/null; then
        zcat "$READ_SAMPLE" > "$output_file"
      else
        cat "$READ_SAMPLE" > "$output_file"
      fi
      FINAL_READS+=("$output_file")
    done
  else
    log_info "Concatenating multiple files..."
    output_file=$TEMP_DIR/${READS[0]%.gz}
    for READ_SAMPLE in "${READS[@]}"; do
      if gzip -t "$READ_SAMPLE" 2>/dev/null; then
        zcat "$READ_SAMPLE" >> "$output_file"
      else
        cat "$READ_SAMPLE" >> "$output_file"
      fi
    done
    FINAL_READS=("$output_file")
  fi
else
  output_file=$TEMP_DIR/${READS[0]%.gz}
  if gzip -t "${READS[0]}" 2>/dev/null; then
    zcat "${READS[0]}" > "$output_file"
  else
    cat "${READS[0]}" > "$output_file"
  fi
  FINAL_READS=("$output_file")
fi

   
############################################
# Deduplication with czid-dedup (optional) 
############################################
if [[ "$DEDUP" == true ]]; then
  log_info "Running czid-dedup on the reads..."
  DEDUP_READS=()
  if [[ "$READ_TYPE" == "pe_short" ]]; then
    input_file_1="${FINAL_READS[0]}"
    input_file_2="${FINAL_READS[1]}"
    output_file_1="$TEMP_DIR/dedup_1.fastq"
    output_file_2="$TEMP_DIR/dedup_2.fastq"
    if [[ -n "$DEDUP_L" ]]; then
      czid-dedup-Linux -i "$input_file_1" -i "$input_file_2" -o "$output_file_1" -o "$output_file_2" -l "$DEDUP_L"
    else
      czid-dedup-Linux -i "$input_file_1" -i "$input_file_2" -o "$output_file_1" -o "$output_file_2"
    fi
    FINAL_READS=("$output_file_1" "$output_file_2")
  else
    # single file
    input_file="${FINAL_READS[0]}"
    output_file="$TEMP_DIR/dedup_single.fastq"
    if [[ -n "$DEDUP_L" ]]; then
      czid-dedup-Linux -i "$input_file" -o "$output_file" -l "$DEDUP_L"
    else
      czid-dedup-Linux -i "$input_file" -o "$output_file"
    fi
    FINAL_READS=("$output_file")
  fi
fi

############################################
# (Optional) Downsampling with rasusa
############################################
# The user may specify coverage+genome_size, or num_bases, or num_reads
if [[ "$DOWNSAMPLE" == true ]]; then
  log_info "Running rasusa for downsampling..."

  # We'll generate new file(s):
  if [[ "$READ_TYPE" == "pe_short" ]]; then
    # we have two files in FINAL_READS
    input_file_1="${FINAL_READS[0]}"
    input_file_2="${FINAL_READS[1]}"
    output_file_1="$TEMP_DIR/ds_1.fastq"
    output_file_2="$TEMP_DIR/ds_2.fastq"
    # Decide which rasusa flags to use:
    RASUSA_CMD=(rasusa reads)
    if [[ $COVERAGE -gt 0 && $GENOME_SIZE -gt 0 ]]; then
      RASUSA_CMD+=(--coverage "$COVERAGE" --genome_size "$GENOME_SIZE")
    elif [[ $NUM_BASES -gt 0 ]]; then
      RASUSA_CMD+=(--bp "$NUM_BASES")
    elif [[ $NUM_READS -gt 0 ]]; then
      RASUSA_CMD+=(--limit "$NUM_READS")
    fi
    RASUSA_CMD+=(-o "$output_file_1" -o "$output_file_2" "$input_file_1" "$input_file_2")
    #Here show commands before execute
    log_info "Executing: ${RASUSA_CMD[*]}"
    # run it
    "${RASUSA_CMD[@]}"
    FINAL_READS=("$output_file_1" "$output_file_2")
  else
    # single file
    input_file="${FINAL_READS[0]}"
    output_file="$TEMP_DIR/ds_${FINAL_READS[0]}"
    RASUSA_CMD=(rasusa reads)
    if [[ $COVERAGE -gt 0 && $GENOME_SIZE -gt 0 ]]; then
      RASUSA_CMD+=(--coverage "$COVERAGE" --genome_size "$GENOME_SIZE")
    elif [[ $NUM_BASES -gt 0 ]]; then
      RASUSA_CMD+=(--bp "$NUM_BASES")
    elif [[ $NUM_READS -gt 0 ]]; then
      RASUSA_CMD+=(--limit "$NUM_READS")
    fi
    RASUSA_CMD+=(-o "$output_file" "$input_file")
    log_info "Executing: ${RASUSA_CMD[*]}"
    "${RASUSA_CMD[@]}"
    FINAL_READS=("$output_file")
  fi
fi

############################################
# Finally, run read2tree
############################################
log_info "Running read2tree --step 2map with read_type=$READ_TYPE, threads=$THREADS"

# If we have pe_short => pass two files as arguments
# Otherwise pass single file
READ2TREE_CMD=( read2tree --step 2map
                --standalone_path marker_genes
                --dna_reference dna_ref.fa
                --threads "$THREADS"
                --output_path read2tree_output
                --debug )

# For the read_type we might do:
case "$READ_TYPE" in
  pe_short)
    READ2TREE_CMD+=(--read_type "short")  # forcing 'short' 
    # add the two read arguments
    READ2TREE_CMD+=(--reads "${FINAL_READS[0]}" "${FINAL_READS[1]}")
    ;;
  se_short)
    READ2TREE_CMD+=(--read_type "short")
    # just one file
    READ2TREE_CMD+=(--reads "${FINAL_READS[0]}")
    ;;
  pacbio)
    READ2TREE_CMD+=(--read_type "long-hifi") 
    READ2TREE_CMD+=(--reads "${FINAL_READS[0]}")
    ;;
  ont)
    READ2TREE_CMD+=(--read_type "long-ont")
    READ2TREE_CMD+=(--reads "${FINAL_READS[0]}")
    ;;
esac

log_info "Executing: ${READ2TREE_CMD[*]}"
"${READ2TREE_CMD[@]}"

log_info "read2tree step 2map completed successfully."

# Optionally do step 3combine
log_info "Running read2tree --step 3combine..."
read2tree --step 3combine \
          --standalone_path marker_genes \
          --dna_reference dna_ref.fa \
          --output_path read2tree_output \
          --tree \
          --debug

log_info "Done. If you want to run iqtree or similar, you can do so now."
