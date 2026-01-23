# Omni2treeView

Omni2treeView is a sub-tool from Omni2tree software to visualize phylogenetic trees with associated metadata in an interactive HTML format. It supports various customization options and is designed to handle large datasets efficiently.


## Features
- Interactive visualization of phylogenetic trees.
- Integration of metadata for enhanced analysis.
- Direct use of Newick tree files and CSV metadata files from Omni2tree.

## Installation

To install Omni2treeView, clone the repository and navigate to the `view` directory:

```bash
git clone git@github.com:DanielPAagustinho/omni2tree.git
cd omni2tree/view
```

## Usage

To generate an interactive tree view, use the following command:

```bash

python3  omni2treeview.py  -n demo/hcmv_tree.nwk -m demo/hcmv_meta.csv -c demo/hcmv_five_letter_taxon.tsv  -t template_v5.html -l hcmv  -o  demo/omni2treeview_hcmv
```

This command will create an HTML file in the specified output directory that visualizes the phylogenetic tree along with the provided metadata.

Options:
- `-n`: Path to the Newick tree file.
- `-m`: Path to the metadata CSV file.
- `-c`: (Optional) Path to the taxon mapping file.
- `-t`: Path to the HTML template file.
- `-l`: Label for the dataset.
- `-o`: Output prefix for the generated HTML file and other related files. 


## Input files:

**Newick tree file (.nwk):** output from Omni2tree.

**Metadata CSV file (.csv):** Metadata in CSV format for each sample. the first column should match the leaf names in the Newick tree file.

**Taxon mapping file (.tsv) (optional):** A TSV file mapping sample IDs to taxon names, if needed. the first column should match the metadata CSV file's first column, the second column should be the name in the tree.

**HTML template file (.html):** An HTML template file for rendering the tree view.


## Output files:

The output will include:
- An HTML file for interactive visualization (html).
- A JSON file for tree structure (tree.json).
- A JSON file for global tree settings (tree_meta.json).
- A integrated CSV file combining tree and metadata information (meta.csv).
