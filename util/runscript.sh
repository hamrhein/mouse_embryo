#!/bin/sh

/bin/echo -n "Downloading STRING database files..."

/usr/bin/wget -q https://stringdb-static.org/download/protein.aliases.v11.0/10090.protein.aliases.v11.0.txt.gz
/usr/bin/wget -q https://stringdb-static.org/download/protein.links.detailed.v11.0/10090.protein.links.detailed.v11.0.txt.gz
/usr/bin/wget -q https://stringdb-static.org/download/protein.actions.v11.0/10090.protein.actions.v11.0.txt.gz

/bin/echo "done"
/bin/echo "Building database"
/bin/echo "-----------------"

/software/hamrhein/build_sql_stringdb_database.py 10090.protein.aliases.v11.0.txt.gz 10090.protein.links.detailed.v11.0.txt.gz 10090.protein.actions.v11.0.txt.gz mus_musculus_stringdb_v11.0.db

/bin/echo "-----------------"

/bin/rm 10090.protein.aliases.v11.0.txt.gz
/bin/rm 10090.protein.links.detailed.v11.0.txt.gz
/bin/rm 10090.protein.actions.v11.0.txt.gz

/bin/mkdir figure2
/bin/mkdir figure10

/bin/echo -n "Building Fig 2 Blossom Plot..."

/software/hamrhein/build_blossom_graph.py /data/peng_bloom_adjacency_matrix.tsv /data/peng_bloom_cluster_size_table.txt --output_base figure2_blossom_graph --output_dir figure2

/bin/echo "done"
/bin/echo -n "Building Fig 4 and 10 10X TF graphs..."

/software/hamrhein/build_10x_tf_graphs.py /data/MouseLimbData.h5 mus_musculus_stringdb_v11.0.db --output_dir figure10

/bin/echo "done"

exit 0
