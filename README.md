# extract_genes_ABRicate
Small script to use ABRicate output to extract genes from genome assemblies, reverse complement if necessary, and print to a file

## Installation

This script needs Python 3 with the Pandas and BioPython libraries, as well as seqtk to run. ABRicate is not necessarily needed, although the ABRicate output should include a `STRAND` column with relevent information.

If you have Miniconda installed (https://docs.conda.io/en/latest/miniconda.html), these dependencies can be easily installed. First clone the directory to your machine:

```
# Clone and enter the directory
git clone https://github.com/boasvdp/extract_genes_ABRicate.git
cd extract_genes_ABRicate

# Create a conda environment with the necessary packages
conda env create -f env.yaml

# Activate the conda environment
conda activate env_extract_genes_ABRicate
```

Alternatively, these commands can be used to install the tools separately through conda (not in a separate environment!):

```
conda install -c conda-forge -c bioconda biopython pandas seqtk
```

## Usage

```
usage: extract_genes_abricate.py [-h] [-a ABRICATE FILE] [-g GENOMES DIR]
                                 [-o OUTPUT DIR] [-s SUFFIX]

Extract genes from genes based on ABRicate output.

optional arguments:
  -h, --help            show this help message and exit
  -a ABRICATE FILE, --abricatefile ABRICATE FILE
                        ABRicate file to parse genes
  -g GENOMES DIR, --genomedir GENOMES DIR
                        directory containing genomes
  -o OUTPUT DIR, --output OUTPUT DIR
                        directory for output
  -s SUFFIX, --suffix SUFFIX
                        Genome assembly file suffix (default: .fasta)
```

**IMPORTANT ASSUMPTIONS**

The script assumes the genome assemblies are named almost exactly as they are provided in the ABRicate output (`#FILE` column). The only thing that may differ is the suffix (default `.fasta`, unless otherwise provided using `--suffix`). The script is also at this time only able to handle a single suffix for genome assemblies at a time. 

If you have identified genes for all genomes in your `genomes/` directory (in which all genome assembly files end with `.fasta`) and your ABRicate output is present in `ABRicate_out/strainA.tsv`, run:

```
python extract_genes_ABRicate.py --abricatefile ABRicate_out/strainA.tsv --genomedir genomes/ --output extracted_genes/
```

## Extended usage

ABRicate files can also be combined to speed up things. To combine all files in ABRicate_out/, e.g. run:

```
cat <(head -n 1 ABRicate_out/strainA.tsv) <(for i in ABRicate_out/*.tsv; do tail -n +2 $i; done) > ABRicate_all.tsv
```

After which the extract_genes_ABRicate.py script has to be run only once:

```
python extract_genes_ABRicate.py --abricatefile ABRicate_all.tsv --genomedir genomes/ --output extracted_genes/
```
