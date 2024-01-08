# Pangenome graph building

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/1000G-ONT-F100-PGGB

PGGB=/home/guarracino/tools/pggb/pggb-13482bd06359a7ad8e3d3e0dd6eb6d9399f26046
ODGI=/home/guarracino/tools/odgi/bin/odgi-861b1c04f5622c5bb916a161c1abe812c213f1a5
WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-0b191bb84ffdfd257354c1aa82a7f1e13dc536d0
```

Run `pggb`:

```shell
mkdir -p $DIR_BASE/graphs
cd $DIR_BASE/graphs

ls $DIR_BASE/partitions/1000G-ONT.100x2+4.chr*.fa.gz | grep chrM -v | grep chr2 | tail -n 1 | while read FASTA; do
    NAME=$(basename $FASTA .fa.gz)
    sbatch -c 96 -p tux --job-name pggb-$NAME --wrap "$PGGB -i $FASTA -o $DIR_BASE/graphs/$NAME -s 10k -p 98 -D /scratch -t 96 -T 48"
done
```

Compressed visualizations:

```shell
ls $DIR_BASE/graphs/*/*.og | grep chrM -v | while read GRAPH; do
    NAME=$(basename $GRAPH | cut -f 1,2,3,4 -d '.')
    echo $NAME
    sbatch -c 48 -p workers --job-name odgi-viz-$NAME --wrap "$ODGI paths -i $GRAPH -L | cut -f 1,2 -d '#' | sort | uniq > $GRAPH.prefixes.txt; $ODGI viz -i $GRAPH -o $GRAPH.viz_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500; $ODGI viz -i $GRAPH -o $GRAPH.viz_pos_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500 -u -d; $ODGI viz -i $GRAPH -o $GRAPH.viz_depth_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500 -m; $ODGI viz -i $GRAPH -o $GRAPH.viz_inv_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500 -z; $ODGI viz -i $GRAPH -o $GRAPH.viz_O_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500 -O; $ODGI viz -i $GRAPH -o $GRAPH.viz_uncalled_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500 -N"
done
```
