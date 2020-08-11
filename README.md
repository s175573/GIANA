# GIANA: Geometry Isometry based ANtigen-specific TCR Alignment

GIANA is for fast alignment of up to 10<sup>7</sup> TCR hypervariable CDR3 sequences. 

GIANA is developed and maintained by [Li lab at UT Southwestern Medical Center](lilab-utsw.org). Please direct your questions regarding GIANA to [Bo Li](bo.li@utsouthwestern.edu).

GIANA is written in Python3, with the following dependencies:

- [Biopython](https://biopython.org)
- [faiss](https://github.com/facebookresearch/faiss)
- [Scikit-learn](https://scikit-learn.org/stable/)

After installing these dependencies, please download the source code of GIANA and the associated TRBV allele data (Imgt_Human_TRBV.fasta). 

## Usage

Type `python GIANA.py -h` to display all the commandline options. 

### 1. Input data format

Input of GIANA is flexible. The first column is kept for CDR3 amino acid sequence. If TRBV allele information is enabled (by default), the second column is required to be TRBV genes. As the TCR-seq data provided by the Adaptive Biotechnologies does not comply with the IMGT format, we provide the R code (ProcessAdaptiveGenes.R) to convert the Adaptive data input to standard format. In the output, GIANA inserts a column between the first and the second column as the cluster IDs. Other columns in the input data may contain any information, and will be kept in the final output. 

### 2. Clustering with TRBV variable gene

The following code performs standard TCR clustering on an input data, using the TRBV allele information:

`python GIANA.py -f tutorial.txt`

GIANA can also be applied to a folder containing a list of input files:

`python GIANA.py -d input_dir/`

### 3. Clustering without TRBV variable gene

Even with TRBV gene as the second column, GIANA can perform clustering without the TRBV allele information:

`python GIANA.py -f tutorial.txt -v`

This option is useful when TRBV gene information is not available (in that case, the input data can contain only one column of CDR3s). The user can choose a more stringent cut-off of Smith-Waterman alignment score (higher score is more stringent):

`python GIANA.py -f tutorial.txt -v -S 4`

### 4. Clustering in non-exact mode

By default, GIANA will run Smith-Waterman clustering for all the pre-clusters identified in the isometric nearest neighbor search, which is referred to as the 'exact' mode. The user can choose to disable the exact mode to gain over 10X computational efficiency, at the cost of less specific TCR clustering:

`python GIANA.py -f tutorial.txt -e`

This mode is not necessary when processing less than 1 million sequences. In non-exact mode, the users are recommended to apply a more stringent isometric distance cut-off to increase specificity:

`python GIANA.py -f tutorial.txt -e -t 5`

Here smaller `-t` value is more stringenet.

