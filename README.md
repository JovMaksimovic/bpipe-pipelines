# bpipe-pipelines
A collection of bpipe pipeline stages and pipeline examples for common genomics workflows. Groovy files with "stages" in the filename contain definitions of bpipe pipeline stages. All other groovy files provide examples of pipelines for common analyses such as rna-seq, atac-seq and chip-seq.

The "bpipe.config" is an example of a bpipe configuration file for running the pipelines on a cluster; specifically on meerkat at MCRI.

For example, the rnaseq pipeline can be run on the meerkat cluster using the following command:

bpipe run -n 100 /path/to/groovy/file/rnaseq.groovy /path/to/input/files/\*.sra
