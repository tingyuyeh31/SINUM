## SIngle-cell Network Using Mutual information (SINUM)
This repository contains the source code for single-cell networks (SCNs) construction by utilizing mutual information inference of gene-gene association.
The example datasets are stored inside `example` folder, as well the example outputs.
### Basic usage
#### sinum.py
```
$ python code/sinum.py -f example/dataset/ChuType_exampleGEM_log2.txt -b 0.2 -z 0.0 -o ./example -n ChuType_example
```
Argument | Variable | Description | Default value
------------ | ------------- | ------------- | -------------
-f | GEM_file | path to gene expression matrix (GEM) file | example/dataset/ChuType_exampleGEM_log2.txt
-b | bs | box size | 0.2
-z | cf | z-score threshold | 0.0
-o | outdir | output directory | ./example
-n | savename | save name | ChuType_example

There are two outputs: adjacency matrices & degree matrix.

The adjacency matrices (gene x gene) for each cell are in sparse matrix format(`.npz`).

The degree matrix (gene x cell) is in `.csv` format.

The outputs are stored in the `outdir` directory, within which two folders `output_adj_network` & `output_degree_matrix` are automatically created.
### Advanced usage example
#### Direct calling of `sinum_funcs` in your code.
```python
import pandas as pd
import sinum_funcs
# read the GEM as SINUM input
# GEM_file is the path to the gene expression matrix txt file
GEM = pd.read_csv(GEM_file,sep='\t',index_col=0)
# define boxsize parameter, default = 0.2
bs = 0.2
# define z-score cutoff as the determination of edges, default = 0.0
cf = 0.0
# give the path to save the result
outdir = "./output"
# give the name of the output file ex: dataset name
savename = "ChuType"

# execute the main SINUM algorithm
sinum_funcs.sinumnet(GEM, outdir, savename, boxsize = bs, cutoff=cf)
```
### Dependencies
The code is written in python3, and not guaranteed to work in python2, espeacially for the package versions. Additionally, the following packages (and the suggested version) must also be installed:
* [NumPy (1.20.1)](https://numpy.org/)
* [Pandas (1.4.2)](https://pandas.pydata.org/)
* [SciPy (1.8.1)](https://scipy.org/)
* [scikit-learn (1.1.1)](https://scikit-learn.org/stable/)
