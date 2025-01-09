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
