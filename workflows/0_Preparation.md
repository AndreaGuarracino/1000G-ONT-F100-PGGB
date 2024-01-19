# Preparation

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/1000G-ONT-F100-PGGB

WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-0b191bb84ffdfd257354c1aa82a7f1e13dc536d0
```

Download assemblies:

```shell
mkdir -p $DIR_BASE/assemblies/F100
cd $DIR_BASE/assemblies/F100

# From https://s3.amazonaws.com/1000g-ont/index.html?prefix=FIRST_100_FREEZE/CARD_assembly_pipeline_data/
cut -f 3 $DIR_BASE/data/1000G-ONT-F100.urls.tsv | parallel -j 4 'wget -q -O - {} | bgzip -@ 11 -l 9 > $(basename {}).gz && echo Compressed {}'
```

PanSN:

```shell
ls $DIR_BASE/assemblies/F100/*fasta.gz | while read FASTA; do
    echo $FASTA
    SAMPLE=$(basename $FASTA .fasta.gz | cut -f 1 -d '.')
    HAP=$(basename $FASTA .fasta.gz | cut -f 3 -d '.' | cut -f 4 -d '_')
    echo $SAMPLE $HAP
    zcat $FASTA | sed "s/^>/>$SAMPLE#$HAP#/g" | bgzip -@ 48 -l 9 > $DIR_BASE/assemblies/F100/$SAMPLE.$HAP.fa.gz
    samtools faidx $DIR_BASE/assemblies/F100/$SAMPLE.$HAP.fa.gz
done

rm $DIR_BASE/assemblies/F100/*fasta.gz
```

Map against 4 references:

```shell
cd $DIR_BASE/assemblies

# Prepare references
zcat \
    /lizardfs/guarracino/pggb-paper/assemblies/primates/chm13v2.0.fa.gz \
    /lizardfs/guarracino/pggb-paper/assemblies/primates/grch38.fa.gz \
    /lizardfs/guarracino/pggb-paper/assemblies/primates/hg002v101.fa.gz | \
    sed -e 's/_MATERNAL//g' -e 's/_PATERNAL//g' | bgzip -@ 48 -l 9 > $DIR_BASE/assemblies/chm13v2+grch38+hg002v101.fa.gz
samtools faidx $DIR_BASE/assemblies/chm13v2+grch38+hg002v101.fa.gz

mkdir -p $DIR_BASE/partitions
cd $DIR_BASE/partitions

ls $DIR_BASE/assemblies/F100/*fa.gz | while read FASTA; do
    NAME=$(basename $FASTA .fa.gz)
    echo $NAME
    
    # -c/-t 24 leads to a few OOMs
    sbatch -c 48 -p allnodes --job-name $NAME --wrap "hostname; $WFMASH $DIR_BASE/assemblies/chm13v2+grch38+hg002v101.fa.gz $FASTA -m -t 48 > $DIR_BASE/partitions/$NAME.mapping.paf"
done
```

Check the length threshold to apply without losing too much sequence:

```shell
# We just want the chromosome information, not the specific reference
sed -e 's/chm13#1#//g' -e 's/grch38#1#//g' -e 's/hg002#P#//g' -e 's/hg002#M#//g' $DIR_BASE/partitions/*.mapping.paf | python3 /lizardfs/guarracino/robertsonian_translocation/scripts/majority_rule_partitioning.py | sort -k 1,1 -k 2,2 -V > $DIR_BASE/partitions/F100.partitions.tsv

# Check how much sequence was partitioned
ls $DIR_BASE/assemblies/F100/*.fai | while read f; do
    DIM1=$(cat $f | awk '{sum+=$2}END{print(sum)}')
    DIM2=$(cat $f | grep -Ff <(cut -f 1 $DIR_BASE/partitions/F100.partitions.tsv) | awk '{sum+=$2}END{print(sum)}')
    DIM3=$(cat $f | awk '$2 >= 30000' | grep -Ff <(cut -f 1 $DIR_BASE/partitions/F100.partitions.tsv) | awk '{sum+=$2}END{print(sum)}') # without short contigs
    N1=$(cat $f | wc -l)
    N2=$(cat $f | grep -Ff <(cut -f 1 $DIR_BASE/partitions/F100.partitions.tsv) | wc -l)
    N3=$(cat $f | awk '$2 >= 30000' | grep -Ff <(cut -f 1 $DIR_BASE/partitions/F100.partitions.tsv) | wc -l)
    RATIO21=$(echo "scale=4; $DIM2 / $DIM1" | bc)
    RATIO31=$(echo "scale=4; $DIM3 / $DIM1" | bc)
    echo $(basename $f .fa.gz.fai) $DIM1 $DIM2 $RATIO21 $RATIO31 $N1 $N2 $N3
done | sort -k 5,5nr

# Worst sample:
# - min. 50 kbps: HG01615.1 2924607582 2917573978 .9975 .9856 3353 3123 985  (max loss: 1.44%)
# - min. 30 kbps: HG03279.1 2899727015 2892122280 .9973 .9904 3183 2949 1250 (max loss: 0.96%)
# - min. 25 kbps: HG03279.1 2899727015 2892122280 .9973 .9918 3183 2949 1396 (max loss: 0.82%)
```

Partition by chromosome, filtering by length:

```shell
# Create partitions
(seq 1 22; echo X; echo Y; echo M) | while read i; do
    CHR=chr$i
    echo $CHR

    # Set output file name conditionally
    if [[ "$CHR" == "chrM" ]]; then
        NAME_CHR_FASTA="1000G-ONT.100x2+4.$CHR.fa"
    else
        NAME_CHR_FASTA="1000G-ONT.100x2+4.$CHR.30kbp.fa"
    fi

    # Add references
    samtools faidx $DIR_BASE/assemblies/chm13v2+grch38+hg002v101.fa.gz $(cut -f 1 $DIR_BASE/assemblies/chm13v2+grch38+hg002v101.fa.gz.fai | grep chr$i -w) > /scratch/$NAME_CHR_FASTA

    # Add F100 assemblies
    sed -e 's/chm13#1#//g' -e 's/grch38#1#//g' -e 's/hg002#P#//g' -e 's/hg002#M#//g' $DIR_BASE/partitions/F100.partitions.tsv | grep chr$i -w | cut -f 1 > $DIR_BASE/partitions/$CHR.contigs.txt

    ls $DIR_BASE/assemblies/F100/*fa.gz | while read FASTA; do       
        if [[ "$CHR" == "chrM" ]]; then
            samtools faidx $FASTA $(cut -f 1 $FASTA.fai | grep -Ff $DIR_BASE/partitions/$CHR.contigs.txt -w) >> /scratch/$NAME_CHR_FASTA
        else
            samtools faidx $FASTA $(awk '$2 >= 30000' $FASTA.fai | cut -f 1 | grep -Ff $DIR_BASE/partitions/$CHR.contigs.txt -w) >> /scratch/$NAME_CHR_FASTA
        fi
    done
    rm $DIR_BASE/partitions/$CHR.contigs.txt

    bgzip -@ 48 -l 9 /scratch/$NAME_CHR_FASTA
    samtools faidx /scratch/$NAME_CHR_FASTA.gz

    mv /scratch/$NAME_CHR_FASTA.gz* $DIR_BASE/partitions/
done
```
