
# Welcome to virus2tree!!

This new version of read2tree enables the creation of a reference database via OMA Standalone using as input coding sequences from NCBI assemblies. The final tree combines the presence of assemblies (that is, the reference) and the read samples. It also supports read deduplication with czid-dedup and downsampling with rasusa, among other functionalities.

## Installation of dependencies

This software relies on four external tools: [OMA Standalone](https://omabrowser.org/standalone/), [Rasusa](https://github.com/mbhall88/rasusa?tab=readme-ov-file#install), [czid-dedup](https://github.com/chanzuckerberg/czid-dedup?tab=readme-ov-file#installation), and [Read2Tree](https://github.com/DessimozLab/read2tree/tree/minimap2?tab=readme-ov-file#installation). It assumes all programs are in your Conda environment or `PATH`. 
Below are two general ways to install all the required dependencies. For more details, please visit the respective web pages.

### 1. Installation with Conda (not possible for all dependencies)

[Conda](https://docs.anaconda.com/miniconda/) is a package manager that allows you to install all dependencies quickly and easily.

```bash
conda create -n my_env python=3.10.8 -y 
conda activate my_env 
conda install -c bioconda rasusa read2tree sra-tools entrez-direct -y
```

**Notes:** 
* OMA standalone and czid-dedup are not available via Conda. Please, follow the "Installation from source" instructions below.
* The Conda version of read2tree does not include the minimap2 branch. If you need this branch, follow the "Installation from source" instructions.

### 2. Installation from source

If you prefer to install the tools manually from their source code, use the following commands:

**OMA Standalone**

```bash
## Download the last version, in this example is 2.6.0
wget -O oma.tgz https://omabrowser.org/standalone/OMA.2.6.0.tgz 
tar xvzf oma.tgz 
cd OMA.2.6.0

## Below choose your install path, if not OMA will be installed in /usr/local/OMA (you might need to use sudo in this case)
./install.sh /your/install/path

## After installation, make sure the bin folder of OMA is in your PATH variable. For that, edit your shell configuration file (`~/.bashrc`, `~/.zshrc`, etc.)
echo 'export PATH=$PATH:/your/install/path/OMA/bin' >> ~/.bashrc 
source ~/.bashrc
```

**Rasusa**

```bash
## When rasusa is downloaded it is automatically added to your PATH
curl -sSL rasusa.mbh.sh | sh
```

**czid-dedup**

Take into account that `czid-dedup` requires [rust/cargo](https://www.rust-lang.org/tools/install) for compilation

```bash
git clone https://github.com/chanzuckerberg/czid-dedup.git 
cd czid-dedup 
cargo build --release 

#Make sure that the release directory is in your PATH variable
echo 'export PATH=$PATH:your/install/path/czid-dedup/target/release' >> ~/.bashrc 
source ~/.bashrc
```

**Read2Tree**

```bash
## Create conda env
conda create -n r2t python=3.10.8 -y 
conda activate r2t

## Get required python packages
conda install -c conda-forge biopython numpy Cython ete3 lxml tqdm scipy pyparsing requests natsort pyyaml filelock -y
conda install -c bioconda dendropy pysam -y

## Install required softwares
conda install -c bioconda mafft iqtree minimap2 samtools -y

## Clone minimap2 branch of read2tree
git clone --branch minimap2 https://github.com/DessimozLab/read2tree.git 
cd read2tree 
python setup.py install
## read2tree will be placed in the default bin folder of your Conda installation
```

**SRA Toolkit**

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz 
tar -xvzf sratoolkit.current-ubuntu64.tar.gz

## Add executable to your path (using your own version, in this case is 3.2.0)
echo 'export PATH="$PATH:/your/install/path/sratoolkit.3.2.0-ubuntu64/bin"' >> ~/.bashrc 
source ~/.bashrc
```

**Entrez direct**

```bash
## Get the scripts and download them in an "edirect folder" in the user's home directory
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
source ~/.bashrc
```

### 3. Verify Installation

To verify that all tools are correctly installed and available in your Conda environment or `PATH`, run the following command:

```bash
oma -h
rasusa --help
czid-dedup --help
read2tree --help
fasterq-dump --help
esearch -h
```

## Installation

You can set up virus2tree cloning the repo and running the installer:

```bash
git clone https://github.com/DanielPAagustinho/virus2tree.git
cd virus2tree
./install.sh /your/install/path
```

The installation script creates symlinks to the shell entry points.
* If you omit the install path, symlinks go to /usr/local/bin (may require sudo).
* If you use a custom path, make sure it’s in your PATH.

Finally, check your installation with:

```bash
which v2t-step1 && v2t-step1 --help
which v2t-step2 && v2t-step2 --help
which v2t-sra   && v2t-sra   --help
```

## Quick start

Minimal example (adjust paths to your data):

```bash
# Step 1: Create reference OMA database for r2t with NCBI accessions from RSV
v2t-step1 -i rsv_accessions.csv -g rsv_outgroups.txt -T 25 --root_dir virus2tree_rsv --out_dir read2tree &> rsv_long_step1.log

# Step 2: Map long nanopore RSV reads to the reference
parallel -j 4 v2t-step2 \
  -r {1} -t ont --dedup --downsample --coverage 250 --genome_size 15kb --root_dir virus2tree_rsv --out_dir read2tree -T 20 ::: \
  $(ls reads/*fastq* | sort) &>> "rsv_long_step2.log" &

# Step 3: Create the tree
read2tree --step 3combine --standalone_path marker_genes --dna_reference dna_ref.fa --output_path virus2tree_rsv/read2tree --tree --debug
```

## Running step 1: Creating the reference database

```bash
v2t-step1 -i rsv_accessions.csv -g rsv_outgroups.txt -T 25 --out_dir read2tree --temp_dir temp --debug &> def_rsv_long.log
```
To create the reference database, two key input files are required:

`-i`, `--input` (Required): A file containing the NCBI accessions to be used for reference (`rsv_accessions.csv`).

`-g`, `--outgroup` (optional but recommended): A file containing taxon(s) to be used as outgroups by OMA Standalone during orthologous group prediction (`rsv_outgroups`).  
If this file is not specified, OMA Standalone will use midpoint rooting, which is likely incorrect and will significantly affect hierarchical orthologous groups (HOGs) inferred by OMA.

### **Command Parameters**

| **Parameter**       | **Description** |
|--------------------|-----------------------------|
| `-i`, `--input`    | **Required.** CSV file with NCBI accessions. |
| `-g`, `--outgroup` | **Optional (recommended)** File with outgroup taxa used by OMA. |
| `--root_dir`       |Root directory where all outputs are written. **Default:** current directory|
| `--out_dir`        | read2tree step 1 output directory (relative to `--root_dir` or absolute). **Default:** `read2tree_output`. |
| `--temp_dir`       | Temporary directory. If relative, it’s resolved under `--root_dir`. **Default:** `/tmp`|
| `--resume_download`       |Skip taxa already downloaded from NCBI into the `db` folder. If all taxa were already downloaded, it resumes from Step 1.4. Additionally, if the required files are already present, Step 1.4 is bypassed and the script practically resumes from the OMA Standalone run (Step 1.6). |
|`--og_min_fraction`| Keep only OGs present in at least this fraction of species (0–1). If omitted, all OGs are kept. |
| `-p, --use_mat_peptides`       | Download GBK files for each taxon's accession(s) and uses the mat_peptide features instead of CDS features if at least one mat_peptide is found. |
| `-q, --use_mat_peptides_only`       | Same as --use_mat_peptides, except that if no mat_peptide feature is found, it does not download CDS features and simply skips that taxon. |
| `-T, --threads`   | Number of threads to use for Oma Standalone and the first step of read2tree. |
| `--debug`         | Keeps temp directory with intermediate files. |
| `-h, --help`         | Show help. |

### **Accession File Format**

The accession file must be a comma-separated values (CSV) text file, with the first line as the header. Each line represents a taxon/species/strain with associated accessions. The format varies depending on whether a five-letter code is included.

#### **Columns:**
1. **First column (required):** Taxon/species/strain name. Header: taxon (or taxa),species or strain(s).
2. **Second column (optional):** Five-letter code. Must be exactly 5 alphanumeric characters. Header: code(s). If not provided, a random five-letter code for each taxon will be generated and saved in the file five_letter_taxon.tsv
3. **Third and onward (required):** One or more accession numbers (comma-separated) to obtain coding sequences. Accepts NCBI Nucleotide database accessions and assembly identifiers (GCF_/GCA_). Header: accession(s).

Commented lines starting with # are ignored.

#### **Example Input Files**

##### **With a five-letter code:**
```plaintext
STRAINS,CODE,accessions
influenza A virus California,INCFA,GCF_001343785.1
Influenza A Hong Kong, INHKA,GCF_000851145.1
ebola virus,EBOLA,GCA_034098425.1
Measles morbillivirus,MEAMO,GCF_000854845.1
Lyssavirus rabies,RABIE,GCF_000859625.1
Mammarenavirus lassaense,MAMMA,GCF_000851705.1
```

##### **Without a five-letter code:**
```plaintext
taxon,accession
chikungunya virus S-27,GCA_000854045.1
SARS-COV 2,NC_045512.2
Norovirus GI,NC_044853.1
Norovirus GIV,NC_044855.1
Norovirus GII,NC_044932.1
Norovirus GIII,NC_029645.1
Norovirus GV, NC_008311.1
```
> ⚠️ **Important**: Only the alphanumeric characters in the taxon column are considered for downstream processing. Taxon names and codes must be unique; duplicates are not allowed.

### **Outgroup File Format**

The outgroup file should contain the taxa fo the accession file to be used as outgroups. It can include one or more taxa.

#### **Example Outgroup File**

```plaintext
influenza A virus California
Influenza A Hong Kong
```

### **Output Files**

| **File**                      | **Description** |
|--------------------------------|------------------------------------------------------------------|
| `db/{taxon}_cds_from_genomic.fna` | Nucleotide FASTA files for each taxon with CDS (or mature peptides if selected) retrieved from NCBI. |
| `DB/{taxon}.fa`               | Amino acid FASTA files for each taxon, prepared for running with OMA Standalone and read2tree. |
| `dna_ref.fa`                  | Reference FASTA file with all nucleotide CDS from all taxa, prepared to be used as input for read2tree. |
| `five_letter_taxon.tsv`        | Table linking taxa with five-letter codes. |
| `parameters.drw`              | Parameter file for the OMA run, modified according to the outgroup file and with the last 4 steps of OMA Standalone deactivated. |
| `Output/`                     | Folder containing the output from OMA Standalone. |
| `marker_genes/`               | Folder required by read2tree with the orthologous groups (OGs) generated by OMA Standalone (contents of `Output/OrthologousGroupsFasta`). |
| `stats/cds_count_per_accession*`                | Per-assembly CDS counts and their distribution across all downloaded assemblies: `cds_count_per_accession.tsv`, `cds_count_per_accession_frequency.tsv` and `cds_count_per_accession.png` |
| `stats/OG_genes.tsv`                | Table with all features for each CDS from the OGs identified by OMA. |
| `stats/OG_genes-unique.tsv`         | Summary table listing the OGs alongside its associated gene, protein, and the taxa in which it is found. |
|`stats/taxon_OG.tsv`| Table containing per-taxon summary: total CDS, missing protein_id, no-OG matches, and matched counts. |
|`stats/OG_taxa.tsv	`|Summary of species coverage per OG and whether it is kept (only when --og_min_fraction is used)|
| `read2tree_output`     | Named according to the `--out_dir` parameter, this folder contains the output of step 1 of read2tree. |

## Running step 2: Processing sample reads and adding them to the read2tree folder

After generating the reference database of orthologous groups, we proceed to add the sample reads.

```bash
#For long nanopore reads
parallel -j 4 v2t-step2 \
  -r {1} -t ont --dedup --downsample --coverage 250 --genome_size 15kb --out_dir read2tree -T 20 ::: \
  $(ls reads/*fastq* | sort) &>> "rsv_long_step2.log" &

#For paired end illumina reads
parallel -j 4 v2t-step2 \
  -r {1} {2} -t pe_short --dedup --downsample --coverage 250 --genome_size 15kb --out_dir read2tree -T 20 ::: \
  $(ls reads/*_1.fastq* | sort) :::+ $(ls reads/*_2.fastq* | sort) &>> "rsv_short_step2.log" &

```

### **Command Parameters**

| **Parameter**      | **Description** |
|--------------------|--------------------------------------------------------------------------------------------------------------------------------|
| `-r, --reads`     | **Required.** Input reads file(s) in `fastq` or `fastq.gz` format. If multiple files are provided and `--read_type` is not `pe_short`, they will be concatenated, assuming they belong to the same sample. |
| `-t, --read_type` | **Required.** Read type: `se_short` (short single-end), `pe_short` (short paired-end), `pacbio`, `ont`. If `pe_short`, two input files are required in `--reads`. |
|`--root_dir`       | Root directory that contains step 1 results; all outputs are written under it. **Default:** current directory.|
| `--out_dir`      | Path to step-1 read2tree output (relative to --root_dir or absolute). **Default:** `read2tree_output`. |
| `--temp_dir`      | Temporary directory. If relative, it’s resolved under `--root_dir`. **Default:** `/tmp`. |
| `--stats_file`   | Name of the summary read statistics file. **Default:** `reads_statistics.tsv` | 
| `--dedup`        | Enables `czid-dedup` to remove duplicate reads. |
| `--dedup_l`      | Prefix length used for deduplication (requires `--dedup`). |
| `--downsample`   | Enables `rasusa` for read subsampling. It is required for all subsampling parameters |
| `--coverage`     | Minimum coverage for subsampling (integer or float: e.g., `250`, `0.1`). Requires `--genome_size`. |
| `--genome_size`  | Genome size for subsampling (integer or with a metric suffix: e.g., `15kb`, `4.1MB`). See [rasusa manual](https://github.com/mbhall88/rasusa?tab=readme-ov-file#genome-size) for more details. Requires `--coverage`. |
| `--num_bases`    | Target number of bases for subsampling (integer). Cannot to be used together with `--genome_size` and `--coverage`. |
| `--num_reads`    | Target number of reads for subsampling (integer). Cannot to be used together with `--genome_size` and `--coverage`. |
| `-T, --threads`   | Threads to use during step 2 of read2tree. **Default:** 4. |
| `--debug`        | Keeps temp directory with intermediate files. |
| `-h, --help`         | Show help. |

### **Output Files**

| **File**                       | **Description** |
|---------------------------------|------------------------------------------------------------------------------------------|
| `read2tree_output`       | Named according to the `--out_dir` parameter. Contains the output of step 2 of read2tree. |
| `reads_statistics.tsv`          | Summary of statistics for processed read samples, including initial state, deduplication, and downsampling. Reports the number of reads, average length, and total bases. |
| `temp/{sample}.fastq`           | Original reads. Uncompressed if initially compressed and concatenated if multiple input files were provided without `pe_short` option. |
| `temp/{sample}_dedup.fastq`     | Deduplicated reads. |
| `temp/{sample}_ds.fastq`        | Downsampled reads. |
| `temp/{sample}_dedup_ds.fastq`  | Deduplicated and downsampled reads. |

## Running step 3: Getting the tree

Finally, we run step 3 of read2tree to generate the tree in .nwk format.

```bash
read2tree --step 3combine --standalone_path marker_genes --dna_reference dna_ref.fa --output_path read2tree --tree --debug
```
## Downloading read samples

To easily get read samples from SRA database, we developed the script `v2t-sra`. The purpose of this script is to facilitate the download and conversion to FASTQ of SRA IDs, whether they correspond to RUN (e.g., SRR, ERR, DRR) or EXPERIMENT (e.g., SRX, ERX, DRX).

Depending on the IDs type, the script proceeds differently:

* RUN mode: The script downloads each RUN individually and converts it to FASTQ.
* EXPERIMENT mode: The script first identifies which RUNs are associated with each EXPERIMENT, and then downloads each resulting RUN.

By default, the script retrieves metadata from the SRA database using `esearch` and `efetch`, then checks column 16 of runinfo to automatically determine whether each RUN is SINGLE or PAIRED. However, if you specify `-l` or `--layout` to force one layout (SINGLE or PAIRED), the script will skip the metadata-based check and apply the specified layout to all downloaded RUNs. However, if your inputs are EXPERIMENT IDs, the script will still need to retrieve runinfo in order to map each experiment to its associated RUNs, even if a forced layout is specified.

The output consists of FASTQ files corresponding to each RUN, renamed to include the species name and, in the case of EXPERIMENTs, also the experiment accession.

### **Command Parameters**

| **Parameter**      | **Description** |
|--------------------|--------------------------------------------------------------------------------------------------------------------------------|
| `-i, --input`     | **Required.** Input file containing SRA IDs with one taxon per line. The format should be `<species_name>,SRA_ID1,SRA_ID2,...`. |
| `-o, --outdir` | Directory where downloaded read files will be saved.	Default: current directory (`pwd`). |
| `-c, --chunk-size`   | Number of SRA IDs per chunk when fetching metadata using esearch and efetch. Default: `350`. |
| `-w, --sleep-secs`      | Number of seconds to sleep between chunked metadata requests to avoid overloading the NCBI server. Default: `1`. |
| `-l, --layout`   | Force sequencing layout (`SINGLE` or `PAIRED`) for all runs. When specified, skips metadata fetching. | 
|`-d, --debug`	| Keep per-species temporary directories and intermediate files. |   
| `-h, --help`        | Displays help information and exits. |

#### **Example Command**

```bash
./v2t-sra -i sra_runs_rsv.csv --chunk-size 3 --sleep-secs 2 --outdir rsv_reads
```

It means to analyze the input file `sra_runs_rsv.csv`, extract the SRA IDs, fetch metadata in chunks of 3 SRA IDs at a time with a 2-second pause between chunks to avoid server overload, and download the resulting reads into the specified output directory `rsv_reads`.

> ⚠️ **Important**: You can adjust `--chunk-size` and `--sleep-secs` to avoid speed issues or overloading of the NCBI server

### **Input File Format**

The input file is a comma-separated values (CSV) text file. No header is required. Each line represents a taxon name (or any identifier of your choice) with one or more SRA IDs separated by commas. 

#### **Columns:**
1. **First column:** Species name (or any identifier). Spaces are allowed but should be avoided for simplicity. This is used for file naming and is sanitized to contain only alphanumeric characters
2. **Second column and onward:** One or more SRA IDs. All IDs in a line must be of the same type (RUN or EXPERIMENT). RUN: IDs starting with SRR, ERR, or DRR. EXPERIMENT: IDs starting with SRX, ERX, or DRX.

Commented lines starting with # are ignored.

#### **Example Input File**

```plaintext
SpeciesA,SRR123456,SRR123457
SpeciesB,SRR999999
# Example of a comment (this line is ignored)
SpeciesC,SRX000111
SpeciesD,ERX222333,ERX222334
```

### **Output Files**

At the end of processing, all generated fastq files are saved in the output directory:

For RUN mode (e.g., SRR123456 from Species A) it generates:
* `SpeciesA_SRR123456_1.fastq` and `SpeciesA_SRR123456_2.fastq` (if layout is PAIRED), or
* `SpeciesA_SRR123456.fastq` (if layout is SINGLE).

> **Note**: If the species name contains spaces (e.g., My Species), they will be deleted.

For EXPERIMENT mode (e.g., SRX000111 with RUN SRR000999 from Species A) it generates:
* `SpeciesA_SRX000111_SRR000999_1.fastq` and `SpeciesA_SRX000111_SRR000999_2.fastq` (if layout is PAIRED), or
* `SpeciesA_SRX000111_SRR000999.fastq` (if layout is SINGLE).

Also, a summary file with details about the SRA IDs downloaded is available:`{outdir}/summary_download.txt` (appends a per-species report + global totals).

At the end of execution, the script removes the directories containing the .sra and metadata files, leaving only the final FASTQ files and the summary file in the specified output directory.

