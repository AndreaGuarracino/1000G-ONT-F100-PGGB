# Pangenome graph building

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/1000G-ONT-F100-PGGB

PGGB=/home/guarracino/tools/pggb/pggb-0317e7f8e20bed2c2c9766fc0227310353889e92
ODGI=/home/guarracino/tools/odgi/bin/odgi-861b1c04f5622c5bb916a161c1abe812c213f1a5
WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-0b191bb84ffdfd257354c1aa82a7f1e13dc536d0
```

Run `pggb`:

```shell
mkdir -p $DIR_BASE/graphs
cd $DIR_BASE/graphs

ls $DIR_BASE/partitions/1000G-ONT.100x2+4.chr*.fa.gz | grep 30kbp | grep chr9 | while read FASTA; do
    NAME=$(basename $FASTA .fa.gz)

    if echo "$FASTA" | grep -q "chrM"; then
        s=500
        p=90
        DIR_OUTPUT=$DIR_BASE/graphs/$NAME.p$p.s$s
    else
        s=10k
        p=98
        DIR_OUTPUT=$DIR_BASE/graphs/$NAME
    fi

    if echo "$FASTA" | grep -q 'chr1\.'; then
        THREADS=96
        POA_THREADS=48
        GROUP=tux
        SEQWISH_K=""
    elif echo "$FASTA" | grep -q 'chr9\.'; then
        s=20k
        THREADS=96
        POA_THREADS=64
        GROUP=tux
        SEQWISH_K="-k 47"
    else
        THREADS=96
        POA_THREADS=80
        GROUP=tux
        SEQWISH_K=""
    fi

    sbatch -c $THREADS -p $GROUP --job-name pggb-$NAME --wrap "hostname; $PGGB -i $FASTA -o $DIR_OUTPUT -s $s -p $p -D /scratch -t $THREADS -T $POA_THREADS $SEQWISH_K --resume"
done
```

Call variants:

```shell
cd $DIR_BASE/graphs

ls $DIR_BASE/partitions/1000G-ONT.100x2+4.chr*.fa.gz | grep 30kbp | grep 'chr9\.' | while read FASTA; do
    NAME=$(basename $FASTA .fa.gz)

    if echo "$FASTA" | grep -q "chrM"; then
        s=500
        p=90
        DIR_OUTPUT=$DIR_BASE/graphs/$NAME.p$p.s$s
        REFS=hg002:10000,chm13:10000
    else
        s=10k
        p=98
        DIR_OUTPUT=$DIR_BASE/graphs/$NAME
        REFS=grch38:10000,chm13:10000
    fi

    if echo "$FASTA" | grep -q 'chr1\.'; then
        SEQWISH_K=""
    elif echo "$FASTA" | grep -q 'chr9\.'; then
        s=20k
        SEQWISH_K="-k 47"
    else
        SEQWISH_K=""
    fi

    sbatch -c 96 -p tux -x octopus03,octopus05 --job-name vg-deconstruct-$NAME --wrap "hostname; alias vcfwave=/home/guarracino/tools/vcflib/build/vcfwave-7c1a31a430d339adcb9a0c2fd3fd02d3b30e3549; $PGGB -i $FASTA -o $DIR_OUTPUT -s $s -p $p -D /scratch -t 96 -T 48 $SEQWISH_K -V $REFS --resume"
done
```

Compressed visualizations:

```shell
cd $DIR_BASE/graphs

ls $DIR_BASE/graphs/*/*.og | sort -V | grep 30kbp | while read GRAPH; do
    NAME=$(basename $GRAPH | cut -f 1,2,3,4 -d '.')
    echo $NAME
    sbatch -c 48 -p workers --job-name odgi-viz-$NAME --wrap "
    $ODGI paths -i $GRAPH -L | cut -f 1,2 -d '#' | sort | uniq > $GRAPH.prefixes.txt;
    $ODGI viz -i $GRAPH -o $GRAPH.viz_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500;
    $ODGI viz -i $GRAPH -o $GRAPH.viz_pos_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500 -u -d;
    $ODGI viz -i $GRAPH -o $GRAPH.viz_depth_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500 -m;
    $ODGI viz -i $GRAPH -o $GRAPH.viz_inv_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500 -z;
    $ODGI viz -i $GRAPH -o $GRAPH.viz_O_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500 -O;
    $ODGI viz -i $GRAPH -o $GRAPH.viz_uncalled_multiqc.merged.png -M $GRAPH.prefixes.txt -x 1500 -y 500 -N"
done
```

```shell
ls 1000G-ONT.100x2+4.*.30kbp/*gfa.zst | while read f; do NAME=$(basename $f); aws s3 cp $f s3://garrisonlab/1000G-ONT-F100-PGGB/$NAME; done
```