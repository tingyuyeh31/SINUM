## SIngle-cell Network Using Mutual information (SINUM)
### Basic usage
#### sinum.py
```
$ python code/sinum.py
```
Argument | Variable | Description | Default value
------------ | ------------- | ------------- | -------------
-f | GEM_file | gene expression matrix (GEM) file | example/dataset/ChuType_exampleGEM_log2.txt
-b | bs | box size | 0.2
-z | cf | z-score threshold | 0.0
-o | outdir | output directory | ./example
-n | savename | save name | ChuType_example

There are two outputs: adjacency matrix & degree matrix.

The adjacency matrices (gene x gene) for each cell are in sparse matrix format(`.npz`).

The degree matrix (gene x cell) is in `.csv` format.

The outputs are stored in the `outdir` directory, within which two folders `output_adj_network` & `output_degree_matrix` are automatically created.
