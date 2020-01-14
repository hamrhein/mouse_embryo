# Generating Blossom and Transcription Factor STRING diagrams

This repository is to provide reproducible figures in full resolution for the
manuscript:

## The changing mouse embryo transcriptome at whole tissue and single-cell resolution

Peng He, Brian A. Williams, Georgi K. Marinov, Diane Trout, Henry Amrhein, Libera Berghella, Say-Tar Goh, Ingrid Plajzer-Frick, Veena Afzal, Len A. Pennacchio, Diane E. Dickel, Axel Visel, Bing Ren, Ross C. Hardison, Yu Zhang and Barbara J. Wold

### Instructions for running locally

1. Clone or download the repository into a working directory

```
# Clone using git
https://github.com/hamrhein/mouse_embryo.git

# Download directly
wget https://github.com/hamrhein/mouse_embryo/archive/master.zip
unzip master.zip
rm master.zip
```

To run the utilities, you'll either need to be in this directory or set your
PYTHONPATH environment variable to reference this directory.

2. Download the Mus musculus annotation from STRING

```
/usr/bin/wget -q https://stringdb-static.org/download/protein.aliases.v11.0/10090.protein.aliases.v11.0.txt.gz
/usr/bin/wget -q https://stringdb-static.org/download/protein.links.detailed.v11.0/10090.protein.links.detailed.v11.0.txt.gz
/usr/bin/wget -q https://stringdb-static.org/download/protein.actions.v11.0/10090.protein.actions.v11.0.txt.gz
```

3. Build the database

```
python3 util/build_sql_stringdb_database.py 10090.protein.aliases.v11.0.txt.gz 10090.protein.links.detailed.v11.0.txt.gz 10090.protein.actions.v11.0.txt.gz mus_musculus_stringdb_v11.0.db
```

4. Remove the downloaded flat files

```
rm 10090.protein.*.v11.0.txt.gz
```

5. Make directories for output files

There are going to be many files produced, so it's best to give them their own
folders.

```
mkdir figure2 figure10
```

6. Build the Blossom graph for Figure 2

```
util/build_blossom_graph.py /data/peng_bloom_adjacency_matrix.tsv /data/peng_bloom_cluster_size_table.txt --output_base figure2_blossom_graph --output_dir figure2
```

7. Build the TF networks for extended Figure 10

```
util/build_10x_tf_graphs.py /data/MouseLimbData.h5 mus_musculus_stringdb_v11.0.db --output_dir figure10
```
