# GIANA: Geometry Isometry based TCR AligNment Algorithm

GIANA is for fast alignment of up to 10<sup>7</sup> TCR hypervariable CDR3 sequences. GIANA applies a mathematical framework to perform isometric encoding of amino acid sequences. This encoding process is describted in the ipynb file.  

GIANA is developed and maintained by [Li lab at University of Pennsylvania](https://lilab-utsw.org). Please direct your questions regarding GIANA to Bo Li: lib3@chop.edu.

GIANA is written in Python3, with the following dependencies:

- [Biopython](https://biopython.org)
- [faiss](https://github.com/facebookresearch/faiss)
- [Scikit-learn](https://scikit-learn.org/stable/)

After installing these dependencies, please download the latest version of GIANA source code (currently v4), query.py and the associated TRBV allele data (Imgt_Human_TRBV.fasta). 

## Usage

Type `python GIANA.py -h` to display all the commandline options:

|Commands|Description|
|--|--|
|`-h, --help`|show this help message and exit|
|`-d DIRECTORY, --directory=DIRECTORY`|Input repertoire sequencing file directory. Please make sure that all the files in the directory are input files.| 
|`-f FILE, --file=FILE`|Input single file of CDR3 sequences for grouping|
|`-F FILES, --fileList=FILES`|Alternative input: a file containing the full path to all the files. If given, overwrite `-d` and `-f` option|
|`-t THR, --threshold=THR`|Isometric distance threshold for calling similar CDR3 groups. Without `-e`, smaller value will increase speed. With `-e`, smaller value will increase specificity. Must be smaller than 12.|
|`-S THR_S, --threshold_score=THR_S`|Threshold for Smith-Waterman alignment score (normalized by CDR3 length). Default 3.6|
|`-G THR_V, --threshold_vgene=THR_V`|Threshold for variable gene comparison. Default 3.7.|                       
|`-o OUTDIR, --output=OUTDIR`|Output directory for intermediate and final outputs.|  
|`-O OUTFILE, --outfile=OUTFILE`|Output file name. If not given, a file with --RotationEncoding will be added to the input file as the output file name.|                     
|`-T ST, --startPosition=ST`|Starting position of CDR3 sequence. The first ST letters are omitted. CDR3 sequence length L must be >= ST+7|
|`-V VFA, --VariableGeneFa=VFA`|IMGT Human beta variable gene sequences|                       
|`-v, --VariableGene`|="If False, GIANA will omit variable gene information and use CDR3 sequences only. This will yield reduced specificity. The cut-off will automatically become the current value-4.0|                     
|`-e, --Exact`|If False, GIANA will implement non-exact mode, and will not perform Smith-Waterman alignment after isometric encoding. Default: True|             
|`-N NN, --NumberOfThreads=NN`|Number of threads for multiple processing. Only applies to faiss search step.|  
|`-M --EncodingMatrix`|If true, GIANA will export the isometric encoding matrix for each TCR. Default: False|
|`-q QUERY, --queryFile=QUERY`|Input query file, if given, GIANA will run in query mode, also need to provide -r option.|
|`-r REF, --refFile=REF`|Input reference file. Only required in the query model.|  
|`-b, --Verbose`|Verbose option: if given, GIANA will print intermediate messages.|                      

Note: GIANA is successfully tested in Linux and Mac OS. In some operation systems, by default `python` is python2. To run GIANA correctly, user needs to use `python3` instead.

### Choosing desirable parameters for TCR clustering
`-S` is a critical parameter that balances clustering precision and recall. When running a new dataset, we recommand users to choose the appropriate cutoff depending on the research purposes. For example, when performing repertoire classification tasks for complex diseases, such as cancer, where multiple disease-associated antigens are anticipated, it is usually more desirable to choose a lenient cutoff to include as many antigen-specific TCRs as possible. When dealing with diseases bearing immunodominant epitopes (such as influenza), the users can usually afford to choose a more stringent cutoff to ensure higher specificity. 

The cutoff for isometric distance (`-t` option), is less important to the outcome of GIANA. By default, we chose a distance cutoff (10) that is large enough to include more dissimilar CDR3s, which will be curated with SW alignment. With a more stringent cutoff, GIANA will become faster, as it processes fewer TCRs, and will achieve higher precision, at the cost of lower sensitivity even with a more lenient cutoff (-S option) for the SW alignment score. In sum, the -t option marks a hard boundary on the TCRs to be passed forward for clustering, where the -S option further refines these TCRs. It is usually not necessary to reduce the default value of isometric distance cutoff unless there is an obvious need for higher (2X) speed. 

### 1. Input data format

Input of GIANA is flexible. The first column is kept for CDR3 amino acid sequence. If TRBV allele information is enabled (by default), the second column is required to be TRBV genes. As the TCR-seq data provided by the Adaptive Biotechnologies does not comply with the IMGT format, we provide the R code (ProcessAdaptiveFile.R) to convert the Adaptive data input to standard format. In the output, GIANA inserts a column between the first and the second column as the cluster IDs. Other columns in the input data may contain any information, and will be kept in the final output. 

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

### 5. Running GIANA in query mode

Assume the reference TCR data is ref.txt. After running clustering (for example, mode 2), GIANA produces a cluster file ref--RotationEncodingBL62.txt. Put this file in the same directory as ref.txt. GIANA will automatically search for this file when running in query mode, for example:

`python GIANA.py -q TestReal-ADIRP0000023_TCRB.tsv -r hc10s10.txt -S 3.3 -o tmp/`

The input query file is designated by `-q` option, which also accepts a file directory, if ending with '/'. Reference file is followed by the `-r` option. 

### 6. Multi-thread processing

GIANA supports multiple thread processing with the `-N` option. To date, this option only applies to the nearest neighbor search step. 

### 7. Visualizing output isometric cooridnates with tree plot

We provide an [example](https://htmlpreview.github.io/?https://github.com/s175573/GIANA/blob/master/GIANAtree.html) to use the neighbor joining tree to visualize the output clusters.
