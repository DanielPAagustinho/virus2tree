# r2t_2025
This is about the development of a Read2tree version that works with viral and metagenomic samples.


# Running OMA

I am running OMA to get a good base to read2tree. Instructions on how to run OMA for viruses here:
```https://github.com/DessimozLab/read2tree/wiki/obtaining-marker-genes-for-viral-dataset```

## Instalation

```bash
cd ~/home/tools
# Downloading sif file from Docker
srun singularity pull docker://dessimozlab/oma_standalone:latest

```

## Data download
First, I created a little script for downloading the cDNA and translated cDNA of all the isolates from the dataset used to create the "New up to 2014" reference genomes. The script can be found at ```/stornext/snfs4/next-gen/scratch/daniel/projects/RSV/p1442/r2t/OMAdb/OMAdataDownloader.sh```
The lists of strains used for RSV/A and /B are in the files ``RSVAlistOfRefs2014.txt`` and ``RSVBlistOfRefs2014.txt`` respectively. I also created a list of 2 bovine RSV as outer groups ``outgroups.txt``. Finally, I created a python script ``protein_converter.py`` that grabs a multi fasta file and convert nucleotide sequence to aminoacid sequences, creating a new fasta. Lastly, I created a ``OMAdataDownloader.sh`` script that uses ``efetch`` to download nucleotide sequences from CDS and then uses ``protein_converter.py`` to create new fasta with AA info.

```bash
cd ~/home/projects/RSV/p1442/r2t/OMAdb
srun OMAdataDownloader.sh RSVAlistOfRefs2014.txt
srun OMAdataDownloader.sh RSVBlistOfRefs2014.txt
srun OMAdataDownloader.sh outgroups.txt
```

## Data cleaning 
First, we need to downoad the script ``clean_fasta_cdna_cds.py`` from ``https://raw.githubusercontent.com/DessimozLab/read2tree/main/archive/scripts/clean_fasta_cdna_cds.py``.

```bash
# Be sure to be in the correct conda environment
conda activate py39
# Downloading and running the script
cd /stornext/snfs4/next-gen/scratch/daniel/projects/RSV/p1442/r2t/OMAdb
wget https://github.com/DessimozLab/read2tree/blob/main/archive/scripts/clean_fasta_cdna_cds.py

srun python clean_fasta_cdna_cds.py db DB db/all_cdna_out.fa
```

## Inferring orthologous groups (gene markers) using OMA standalone.

```bash
# Creating a parameter file
cd ~/home/projects/RSV/p1442/r2t/OMAdb
~/home/tools/OMA.2.6.0/bin/oma -p
```
This will create the ``parameters.drw`` file in the folder. Now you have to access it and add the outgroups five letter codes to line 164 of that file. It should now be:
``OutgroupSpecies := ['s0020','s0160']`` in our case. You can see the codes that needed to be added in the file ```OMAdb/outgroups.txt``` made by me.
This was giving me a bunch of errors, so I changed it to ``OutgroupSpecies := ['s0020']`` per Sina's suggestion, and that solved the problem.

Now, be sure that the you are in the folder containing both the ``DB folder`` and the ``parameters.drw`` file. Being sure of that, run the OMA command as so:

```bash
sbatch OMAcaller.sh
# Moving the result folder to the main r2t working directory
cd ~home/projects/RSV/p1442/r2t
cp -r OMAdb/Output/OrthologousGroupsFasta marker_genes_RSV2014isolates

# Calling read2tree to create the reference.
cat  marker_genes_RSV2014isolates/*.fa > dna_ref.fa
srun -p medium -A proj-gm0001 --job-name="create_references" read2tree --standalone_path ./marker_genes_RSV2014isolates --output_path ./output --reference --dna_reference OMAdb/db/all_cdna_out.fa
```
That scrit contains basically one line, calling OMA. I need to run on sbatch because it talkes too long to do so using srun.


# Progress Adrián

This new version of **read2tree** enables the creation of a reference database via **OMA Standalone** using coding sequences from **NCBI**. It also supports read deduplication with **czid-dedup** and downsampling with **rasusa**.

## Installation of dependencies

This software relies on four external tools: [OMA Standalone](https://omabrowser.org/standalone/), [Rasusa](https://github.com/mbhall88/rasusa?tab=readme-ov-file#install), [czid-dedup](https://github.com/chanzuckerberg/czid-dedup?tab=readme-ov-file#installation), and [Read2Tree](https://github.com/DessimozLab/read2tree/tree/minimap2?tab=readme-ov-file#installation). It assumes all programs are in your Conda environment or `PATH`. 
Below are two general ways to install all the required dependencies. For more details, please visit the respective web pages.

### 1. Installation with Conda

[Conda](https://docs.anaconda.com/miniconda/) is a package manager that allows you to install all dependencies quickly and easily.

```bash
conda create -n my_env python=3.10.8 -y 
conda activate my_env && conda install -c bioconda oma rasusa read2tree -y
```
**Notes:** 
* `czid-dedup` is not available on Conda. See the "From Source" section below for installation.
* Installing `read2tree` via Conda does not include the minimap2 branch. If you need this branch, follow the "From Source" instructions.

### 2. Installation from source

If you prefer to install the tools manually from their source code, use the following commands:

**OMA Standalone**

```bash
wget -O oma.tgz https://omabrowser.org/standalone/OMA.2.6.0.tgz && tar xvzf oma.tgz && cd OMA.2.6.0 && ./install.sh /your/install/path
```

After installation, make sure the bin folder of OMA is in your PATH variable. Add the following line to your shell configuration file (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
export PATH=$PATH:/your/install/path/OMA/bin
```

Then, reload your shell configuration to apply the changes:

```bash
source ~/.bashrc  # or source ~/.zshrc
```

**Rasusa**
```bash
curl -sSL rasusa.mbh.sh | sh
```

**czid-dedup**

`czid-dedup` requires [rust/cargo](https://www.rust-lang.org/tools/install) for compilation

```bash
git clone https://github.com/chanzuckerberg/czid-dedup.git && cd czid-dedup && cargo build --release 
```
As with OMA, make sure that the executable is in your PATH variable. Add the following line to your shell configuration file (`~/.bashrc`, `~/.zshrc`, etc.):
```bash
export PATH=$PATH:your/install/path/czid-dedup/target/release/czid-dedup
```
Then, reload your shell configuration to apply the changes:
```bash
source ~/.bashrc  # or source ~/.zshrc
```
**Read2Tree**
```bash
## Create conda env
conda create -n r2t python=3.10.8 -y && conda activate r2t

## Get required python packages
conda install -c conda-forge biopython numpy Cython ete3 lxml tqdm scipy pyparsing requests natsort pyyaml filelock -y
conda install -c bioconda dendropy pysam -y

## Install required softwares
conda install -c bioconda mafft iqtree minimap2 samtools -y

## Clone minimap2 branch of read2tree
git clone --branch minimap2 https://github.com/DessimozLab/read2tree.git && cd read2tree && python setup.py install
```

### 3. Verify Installation

To verify that all tools are correctly installed and available in your Conda environment or `PATH`, run the following command:

```bash
oma -h && rasusa --help && czid-dedup --help && read2tree --help
```

## Running step 1: Creating the reference database

```bash
virus2tree_step1.sh -i rsv_accessions.txt -o rsv_outgroups.txt -T 25 --out_dir read2tree --temp_dir temp --debug &> def_rsv_long.log &
```
To create the reference database, two key input files are required:

`-i`, `--input` (mandatory): A file containing the NCBI accessions to be used for reference (`rsv_accessions.txt`).

`-o`, `--outgroup` (optional but recommended): A file containing taxon(s) to be used as outgroups by OMA Standalone during orthologous group detection (`rsv_outgroups.txt`).  
If this file is not specified, OMA Standalone will use midpoint rooting, which is likely incorrect and will significantly affect hierarchical orthologous groups (HOGs) inferred by OMA.

### **Command Parameters**

| **Parameter**       | **Description** |
|--------------------|-----------------------------|
| `-i`, `--input`    | **Mandatory.** Input file with NCBI accessions (CSV format). |
| `-o`, `--outgroup` | **Optional.** File with taxon names used as outgroups. |
| `-T`, `--threads`  | Number of threads to use for OmaStandalone and the first step of read2tree. |
| `--out_dir`        | Path to place the output of read2tree. **Default:** `read2tree_output`. |
| `--temp_dir`       | Directory for intermediate temporary files. |
| `--debug`         | Prevents the removal of the temporary directory upon script termination. |

### **Accession File Format**

The accession file must be a **comma-separated values (CSV) text file**, con la primera línea como encabezado. Each line represents a taxon/species/strain with associated accessions. The format varies depending on whether a five-letter code is included.

#### **Columns:**
1. **First column (mandatory):** Taxon/species/strain name. Header: taxon(s),species or strain(s).
2. **Second column (optional):** Five-letter code (must be exactly 5 alphanumeric characters). Header: code(s). If not provided, a random five-letter code for each taxa will be generated and saved in the file five_letter_taxon.tsv
3. **Third and onward (mandatory):** One or more accession numbers (comma-separated) to obtain coding sequences. Accepts NCBI Nucleotide database accessions and assembly identifiers (GCF_/GCA_). Header: accession(s).

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
| `db/{taxon}_cds_from_genomic.fna` | Nucleotide FASTA files for each taxon with coding sequences retrieved from NCBI. |
| `DB/{taxon}.fa`               | Amino acid FASTA files for each taxon, prepared for running with OMA Standalone and read2tree. |
| `dna_ref.fa`                  | Reference FASTA file with all nucleotide coding sequences for all taxa, prepared to be used as input for read2tree. |
| `five_letter_taxon.tsv`        | Table linking taxa with five-letter codes. |
| `parameters.drw`              | Parameter file for the OMA run, modified according to the input outgroup file and with the last 4 steps of OMA Standalone deactivated. |
| `Output/`                     | Folder containing the direct output from OMA Standalone. |
| `marker_genes/`               | Folder containing orthologous groups generated by OMA Standalone (contents of `Output/OrthologousGroupsFasta`). |
| `OG_genes.tsv`                | Summary file with all the features of each gene from each orthologous group identified by OMA. |
| `OG_genes-unique.tsv`         | Summary file listing each OG alongside its associated gene, protein, and the taxa in which it is found. |
| `read2tree output folder`     | Named according to the `--out_dir` parameter, this folder contains the output of step 1 of read2tree. |

## Running step 2: Processing sample reads and adding them to the read2tree folder

After generating the reference database of orthologous groups, we proceed to add the sample reads.

```bash
#For long nanopore reads
parallel -j 4 virus2tree_step2.sh \
  -r {1} -t ont --dedup --downsample --coverage 250 --genome_size 15kb --out_dir read2tree -T 20 ::: \
  $(ls reads/*fastq* | sort) &>> "rsv_long_step2.log" &

#For paired end illumina reads
parallel -j 4 virus2tree_step2.sh \
  -r {1} {2} -t pe_short --dedup --downsample --coverage 250 --genome_size 15kb --out_dir read2tree -T 20 ::: \
  $(ls reads/*_1.fastq* | sort) :::+ $(ls reads/*_2.fastq* | sort) &>> "rsv_short_step2.log" &

```

### **Command Parameters**

| **Parameter**      | **Description** |
|--------------------|--------------------------------------------------------------------------------------------------------------------------------|
| `-r, --reads`     | **Mandatory.** Input reads file(s) in `fastq` or `fastq.gz` format. If multiple files are provided and `--read_type` is not `pe_short`, they will be concatenated, assuming they belong to the same sample. |
| `-t, --read_type` | **Mandatory.** Read type: `se_short` (short single-end), `pe_short` (short paired-end), `pacbio`, `ont`. If `pe_short`, two input files are required in `--reads`. |
| `-T, --threads`   | Number of cores to use during step 2 of read2tree. |
| `--temp_dir`      | Directory for temporary files. |
| `--stats_file`   | Name of the summary read statistics file |    
| `--debug`        | Prevents the removal of the temporary directory upon script termination. |
| `--dedup`        | Enables `czid-dedup` to remove duplicate reads. |
| `--dedup_l`      | Prefix length used for deduplication (requires `--dedup`). |
| `--downsample`   | Enables `rasusa` for read subsampling. It is required for all subsampling parameters |
| `--coverage`     | Minimum coverage for subsampling (integer or float: e.g., `250`, `0.1`). Requires `--genome_size`. |
| `--genome_size`  | Genome size for subsampling (integer or with a metric suffix: e.g., `15kb`, `4.1MB`). See [rasusa manual](https://github.com/mbhall88/rasusa?tab=readme-ov-file#genome-size) for more details. Requires `--coverage`. |
| `--num_bases`    | Target number of bases for subsampling (integer). Does not require `--genome_size` and `--coverage`. |
| `--num_reads`    | Target number of reads for subsampling (integer). Does not require `--genome_size` and `--coverage`. |
| `--out_dir`      | Directory containing read2tree Step 1 output. |

### **Output Files**

| **File**                       | **Description** |
|---------------------------------|------------------------------------------------------------------------------------------|
| `read2tree output folder`       | Named according to the `--out_dir` parameter. Contains the output of step 2 of read2tree. |
| `reads_statistics.txt`          | Summary of statistics for processed read samples, including initial state, deduplication, and downsampling. Reports the number of reads, average length, and total bases. |
| `temp/{sample}.fastq`           | Original reads. Uncompressed if initially compressed and concatenated if multiple input files were provided without `pe_short` option. |
| `temp/{sample}_dedup.fastq`     | Deduplicated reads. |
| `temp/{sample}_ds.fastq`        | Downsampled reads. |
| `temp/{sample}_dedup_ds.fastq`  | Deduplicated and downsampled reads. |

## Running step 3: Getting the tree

Finally, we run step 3 of read2tree to generate the tree in .nwk format.

```bash
read2tree --step 3combine --standalone_path marker_genes --dna_reference dna_ref.fa --output_path read2tree --tree --debug
```

