# Pangenome graph building

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/1000G-ONT-F100-PGGB

PGGB=/home/guarracino/tools/pggb/pggb-8131f28c3085a241b88a4a4de2ae3b64c2ddaaf0
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

    # Check if FASTA contains "chr1." or "chr2."
    if echo "$FASTA" | grep -q 'chr1\.'; then
        THREADS=96
        POA_THREADS=48
        GROUP=tux
        SEQWISH_K=""
    elif echo "$FASTA" | grep -q 'chr9\.'; then
        s=20k
        THREADS=48
        POA_THREADS=24
        GROUP=1tbmem
        SEQWISH_K="-k 47"
    else
        THREADS=96
        POA_THREADS=80
        GROUP=tux
        SEQWISH_K=""
    fi

    sbatch -c $THREADS -p    --job-name pggb-$NAME --wrap "hostname; $PGGB -i $FASTA -o $DIR_OUTPUT -s $s -p $p -D /scratch -t $THREADS -T $POA_THREADS $SEQWISH_K --resume"
done
```

Call variants:

```shell
cd $DIR_BASE/graphs

ls $DIR_BASE/partitions/1000G-ONT.100x2+4.chr*.fa.gz | grep 30kbp | grep 'chr9\.' -v | while read FASTA; do
    NAME=$(basename $FASTA .fa.gz)

    if echo "$FASTA" | grep -q "chrM"; then
        s=5k
        p=98
        REFS=hg002:10000,chm13:10000
        DIR_OUTPUT=$DIR_BASE/graphs/$NAME.p$p.s$s
    else
        s=10k
        p=98
        REFS=grch38:10000,chm13:10000
        DIR_OUTPUT=$DIR_BASE/graphs/$NAME
    fi

    echo sbatch -c 48 -p allnodes -x octopus03,octopus05 --job-name vg-deconstruct-$NAME --wrap "hostname; alias vcfwave=/home/guarracino/tools/vcflib/build/vcfwave-7c1a31a430d339adcb9a0c2fd3fd02d3b30e3549; $PGGB -i $FASTA -o $DIR_OUTPUT -s $s -p $p -D /scratch -t 48 -T 48 -V $REFS --resume"
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

## Comparisons

Prepare HPRC graphs:

```shell
mkdir -p $DIR_BASE/graphs_hprc
cd $DIR_BASE/graphs_hprc

wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2021_11_16_pggb_wgg.88/chroms/chr1.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa.gz
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2021_11_16_pggb_wgg.88/chroms/chr6.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa.gz
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2021_11_16_pggb_wgg.88/chroms/chr22.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa.gz

gunzip *.gfa.gz

ls *.gfa | while read GFA; do
    NAME=$(basename $GFA .gfa)
    $ODGI build -g $GFA -o $DIR_BASE/graphs_hprc/$NAME.og -t 48 -P
done

rm *gfa
```

Prepare ranges:

```shell
mkdir -p $DIR_BASE/loci
cd $DIR_BASE/loci

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
zgrep 'gene_id "C4A"\|gene_id "C4B"' hg38.ncbiRefSeq.gtf.gz |
  awk '$1 == "chr6"' | cut -f 1,4,5 |
  bedtools sort | bedtools merge -d 15000 | bedtools slop -l 10000 -r 20000 -g hg38.chrom.sizes |
  sed 's/chr6/grch38#1#chr6/g' > hg38.ncbiRefSeq.C4.coordinates.bed
rm hg38.chrom.sizes hg38.ncbiRefSeq.gtf.gz

wget -c https://raw.githubusercontent.com/pangenome/odgi/master/test/chr6.HLA_genes.bed
sed 's/grch38#chr6/grch38#1#chr6/g' chr6.HLA_genes.bed -i
bedtools merge -i chr6.HLA_genes.bed -d 10000000 > chr6.MHC.bed
```

Graphs:

```shell
# Function to process graph
process_graph() {
    local locus=$1
    local input_graph_ont=$2
    local input_graph_hprc=$3
    local region=$4

    mkdir -p "$DIR_BASE/loci/$locus"

    NAME="${input_graph_ont%%.fa.gz*}"
    NAME=$(basename $NAME)
    $ODGI extract -i "$input_graph_ont" -r "$region" -o - -t 48 -P | \
        $ODGI sort -i - -o "$DIR_BASE/loci/$locus/$NAME.$locus.og" -O -p gYs -t 48 -P --temp-dir /scratch

    CHR=$(basename $input_graph_hprc | cut -f 1 -d '.')
    $ODGI extract -i "$input_graph_hprc" -r $(echo "$region" | sed 's/#1#/#/g') -o - -t 48 -P | \
        $ODGI sort -i - -o "$DIR_BASE/loci/$locus/HPRC.$CHR.$locus.og" -O -p gYs -t 48 -P --temp-dir /scratch
}

# Process graphs
process_graph "RHD" \
              "$DIR_BASE/graphs/1000G-ONT.100x2+4.chr1.30kbp/1000G-ONT.100x2+4.chr1.30kbp.fa.gz.a8a102b.eb0f3d3.b691e61.smooth.final.og" \
              "$DIR_BASE/graphs_hprc/chr1.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og" \
              "grch38#1#chr1:25272393-25330445"
process_graph "RHCE" \
              "$DIR_BASE/graphs/1000G-ONT.100x2+4.chr1.30kbp/1000G-ONT.100x2+4.chr1.30kbp.fa.gz.a8a102b.eb0f3d3.b691e61.smooth.final.og" \
              "$DIR_BASE/graphs_hprc/chr1.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og" \
              "grch38#1#chr1:25362249-25430192"
process_graph "RHD+RHCE" \
              "$DIR_BASE/graphs/1000G-ONT.100x2+4.chr1.30kbp/1000G-ONT.100x2+4.chr1.30kbp.fa.gz.a8a102b.eb0f3d3.b691e61.smooth.final.og" \
              "$DIR_BASE/graphs_hprc/chr1.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og" \
              "grch38#1#chr1:25272393-25430192"

process_graph "MHC" \
              "$DIR_BASE/graphs/1000G-ONT.100x2+4.chr6.30kbp/1000G-ONT.100x2+4.chr6.30kbp.fa.gz.a8a102b.eb0f3d3.b691e61.smooth.final.og" \
              "$DIR_BASE/graphs_hprc/chr6.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og" \
              "grch38#1#chr6:29722775-33089696"

process_graph "C4" \
              "$DIR_BASE/graphs/1000G-ONT.100x2+4.chr6.30kbp/1000G-ONT.100x2+4.chr6.30kbp.fa.gz.a8a102b.eb0f3d3.b691e61.smooth.final.og" \
              "$DIR_BASE/graphs_hprc/chr6.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og" \
              "grch38#1#chr6:31972057-32055418"

process_graph "CYP2D6" \
              "$DIR_BASE/graphs/1000G-ONT.100x2+4.chr22.30kbp/1000G-ONT.100x2+4.chr22.30kbp.fa.gz.a8a102b.eb0f3d3.b691e61.smooth.final.og" \
              "$DIR_BASE/graphs_hprc/chr22.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og" \
              "grch38#1#chr22:42126499-42130865"
process_graph "CYP2D7" \
              "$DIR_BASE/graphs/1000G-ONT.100x2+4.chr22.30kbp/1000G-ONT.100x2+4.chr22.30kbp.fa.gz.a8a102b.eb0f3d3.b691e61.smooth.final.og" \
              "$DIR_BASE/graphs_hprc/chr22.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og" \
              "grch38#1#chr22:42140203-42149500"
process_graph "CYP2D6+CYP2D7" \
              "$DIR_BASE/graphs/1000G-ONT.100x2+4.chr22.30kbp/1000G-ONT.100x2+4.chr22.30kbp.fa.gz.a8a102b.eb0f3d3.b691e61.smooth.final.og" \
              "$DIR_BASE/graphs_hprc/chr22.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og" \
              "grch38#1#chr22:42126499-42149500"

# Visualizations
ls $DIR_BASE/loci/*/*.og |grep CYP2D6+CYP2D7 |  while read GRAPH; do
    DIR_PARENT="$(dirname $GRAPH)"
    echo $DIR_PARENT $GRAPH
    NAME=$(basename $GRAPH .og)
    $ODGI viz -i $GRAPH -o $DIR_PARENT/$NAME.png -m -B Spectral:4
    $ODGI layout -i $GRAPH -o $DIR_PARENT/$NAME.lay -t 48 --temp-dir /scratch -P
    $ODGI draw -i $GRAPH -c $DIR_PARENT/$NAME.lay -p $DIR_PARENT/$NAME.2D.png
done

# New subgraphs
cd $DIR_BASE/loci
ls $DIR_BASE/loci/*/*.og | grep RHD+RHCE | while read GRAPH; do
    DIR_PARENT="$(dirname $GRAPH)"
    echo $DIR_PARENT $GRAPH
    NAME=$(basename $GRAPH .og)
    $ODGI paths -i $GRAPH -f | bgzip -@ 48 -l 9 > $DIR_PARENT/$NAME.fa.gz && samtools faidx $DIR_PARENT/$NAME.fa.gz
    sbatch -c 48 -p workers -x octopus03,octopus05 --job-name pggb-$NAME --wrap "$PGGB -i $DIR_PARENT/$NAME.fa.gz -o $DIR_PARENT/pggb.$NAME.p90.s1k.k23 -p 90 -s 1k -k 23 -D /scratch/pggb.$NAME.tmp"
    sbatch -c 48 -p workers -x octopus03,octopus05 --job-name pggb-$NAME --wrap "$PGGB -i $DIR_PARENT/$NAME.fa.gz -o $DIR_PARENT/pggb.$NAME.p90.s1k.k47 -p 90 -s 1k -k 47 -D /scratch/pggb.$NAME.tmp"
    # sbatch -c 48 -p workers -x octopus03,octopus05 --job-name pggb-$NAME --wrap "$PGGB -i $DIR_PARENT/$NAME.fa.gz -o $DIR_PARENT/pggb.$NAME.p90.s5k.k23 -p 90 -s 5k -k 23 -D /scratch/pggb.$NAME.tmp"
    # sbatch -c 48 -p workers -x octopus03,octopus05 --job-name pggb-$NAME --wrap "$PGGB -i $DIR_PARENT/$NAME.fa.gz -o $DIR_PARENT/pggb.$NAME.p90.s5k.k47 -p 90 -s 5k -k 47 -D /scratch/pggb.$NAME.tmp"
    # sbatch -c 48 -p workers -x octopus03,octopus05 --job-name pggb-$NAME --wrap "$PGGB -i $DIR_PARENT/$NAME.fa.gz -o $DIR_PARENT/pggb.$NAME.p98.s5k.k23 -p 98 -s 5k -k 23 -D /scratch/pggb.$NAME.tmp"
    # sbatch -c 48 -p workers -x octopus03,octopus05 --job-name pggb-$NAME --wrap "$PGGB -i $DIR_PARENT/$NAME.fa.gz -o $DIR_PARENT/pggb.$NAME.p98.s5k.k47 -p 98 -s 5k -k 47 -D /scratch/pggb.$NAME.tmp"
done
```

Mixed graphs:

```shell
ls $DIR_BASE/loci/*/*.og | rev | cut -f 2 -d '/' | rev | sort | uniq | while read LOCUS; do
    echo $LOCUS
    
    PATH_1000G_FASTA=$DIR_BASE/loci/$LOCUS/1000G-ONT.100x2+4.*.30kbp.$LOCUS.fa.gz
    PATH_1000G_FASTA=$(eval echo $PATH_1000G_FASTA) # force '*' expansion
    
    PATH_HPRC_FASTA=$DIR_BASE/loci/$LOCUS/HPRC.*.$LOCUS.fa.gz
    PATH_HPRC_FASTA=$(eval echo $PATH_HPRC_FASTA) # force '*' expansion

    CHR=$(basename $PATH_HPRC_FASTA .fa.gz | cut -f 2 -d '.')
    
    # Do not take chm13 and grch38 twice
    cat \
        <(zcat $PATH_1000G_FASTA) \
        <(samtools faidx $PATH_HPRC_FASTA $(cut -f 1 $PATH_HPRC_FASTA.fai | grep 'chm13\|grch38' -v)) | \
        bgzip -@ 48 -l 9 > $DIR_BASE/loci/$LOCUS/1000G-30kbp+HPRC.$CHR.$LOCUS.fa.gz && samtools faidx $DIR_BASE/loci/$LOCUS/1000G-30kbp+HPRC.$CHR.$LOCUS.fa.gz
done

# pggb
NAME=1000G-30kbp+HPRC.chr6.C4
$PGGB -i $DIR_BASE/loci/C4/$NAME.fa.gz -o $DIR_BASE/loci/C4/pggb.$NAME.p90.s5k.k23             -D /scratch/pggb.$NAME.tmp

NAME=1000G-30kbp+HPRC.chr6.MHC
sbatch -c 48 -p workers -x octopus03,octopus05 --job-name pggb-$NAME --wrap "$PGGB -i $DIR_BASE/loci/MHC/$NAME.fa.gz -o $DIR_BASE/loci/MHC/pggb.$NAME.p90.s5k.k23             -D /scratch/pggb.$NAME.tmp"

# pca/tree
GRAPH=$DIR_BASE/loci/MHC/pggb.$NAME.p90.s5k.k23/1000G-30kbp+HPRC.*.smooth.final.og
GRAPH=$(eval echo $GRAPH) # force '*' expansion
$ODGI similarity -i $GRAPH -d -D '#' -p 1 -t 48 > $DIR_BASE/loci/MHC/pggb.$NAME.p90.s5k.k23/$NAME.sample.dist.tsv
$ODGI similarity -i $GRAPH -d -D '#' -p 2 -t 48 > $DIR_BASE/loci/MHC/pggb.$NAME.p90.s5k.k23/$NAME.haplotype.dist.tsv
```
