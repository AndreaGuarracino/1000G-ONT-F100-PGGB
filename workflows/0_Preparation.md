# Preparation

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/1000G-ONT-F100-PGGB

PGGB=/home/guarracino/tools/pggb/pggb-13482bd06359a7ad8e3d3e0dd6eb6d9399f26046
ODGI=/home/guarracino/tools/odgi/bin/odgi-861b1c04f5622c5bb916a161c1abe812c213f1a5
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

Partition by chromosome, filtering by length:

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

# Map assemblies against references (-t 24 leads to a few OOMs)
ls $DIR_BASE/assemblies/F100/*fa.gz | while read FASTA; do
    NAME=$(basename $FASTA .fa.gz)
    echo $NAME
    
    sbatch -c 48 -p allnodes --job-name $NAME --wrap "hostname; $WFMASH $DIR_BASE/assemblies/chm13v2+grch38+hg002v101.fa.gz $FASTA -m -t 48 > $DIR_BASE/partitions/$NAME.mapping.paf"
done

# We just want the chromosome information, not the specific reference
sed -e 's/chm13#1#//g' -e 's/grch38#1#//g' -e 's/hg002#P#//g' -e 's/hg002#M#//g' $DIR_BASE/partitions/*.mapping.paf | python3 /lizardfs/guarracino/robertsonian_translocation/scripts/majority_rule_partitioning.py | sort -k 1,1 -k 2,2 -V > $DIR_BASE/partitions/F100.partitions.tsv

# Check how much sequence was partitioned
ls $DIR_BASE/assemblies/F100/*.fai | while read f; do
    DIM1=$(cat $f | awk '{sum+=$2}END{print(sum)}')
    DIM2=$(cat $f | grep -Ff <(cut -f 1 $DIR_BASE/partitions/F100.partitions.tsv) | awk '{sum+=$2}END{print(sum)}')
    DIM3=$(cat $f | awk '$2 >= 50000' | grep -Ff <(cut -f 1 $DIR_BASE/partitions/F100.partitions.tsv) | awk '{sum+=$2}END{print(sum)}') # without short contigs
    N1=$(cat $f | wc -l)
    N2=$(cat $f | grep -Ff <(cut -f 1 $DIR_BASE/partitions/F100.partitions.tsv) | wc -l)
    N3=$(cat $f | awk '$2 >= 50000' | grep -Ff <(cut -f 1 $DIR_BASE/partitions/F100.partitions.tsv) | wc -l)
    RATIO21=$(echo "scale=4; $DIM2 / $DIM1" | bc)
    RATIO31=$(echo "scale=4; $DIM3 / $DIM1" | bc)
    echo $(basename $f .fa.gz.fai) $DIM1 $DIM2 $RATIO21 $RATIO31 $N1 $N2 $N3
done | sort -k 5,5nr

# Worst sample:
# - min. 50 kbps: HG01615.1 2924607582 2917573978 .9975 .9856 3353 3123 985  (max loss: 1.44%)
# - min. 30 kbps: HG03279.1 2899727015 2892122280 .9973 .9904 3183 2949 1250 (max loss: 0.96%)
# - min. 25 kbps: HG03279.1 2899727015 2892122280 .9973 .9918 3183 2949 1396 (max loss: 0.82%)

# Create partitions
(seq 1 22; echo X; echo Y; echo MT) | while read i; do
    CHR=chr$i
    echo $CHR

    sed -e 's/chm13#1#//g' -e 's/grch38#1#//g' -e 's/hg002#P#//g' -e 's/hg002#M#//g' $DIR_BASE/partitions/F100.partitions.tsv | grep chr$i -w | cut -f 1 > $DIR_BASE/partitions/$CHR.contigs.txt

    # Add references
    samtools faidx $DIR_BASE/assemblies/chm13v2+grch38+hg002v101.fa.gz $(cut -f 1 $DIR_BASE/assemblies/chm13v2+grch38+hg002v101.fa.gz.fai | grep chr1 -w) > $DIR_BASE/partitions/1000G-ONT.100x2+4.$CHR.50bkp.fa

    # Add F100 assemblies
    ls $DIR_BASE/assemblies/F100/*fa.gz | while read FASTA; do       
        if [[ "$CHR" == "chrM" ]]; then
            samtools faidx $FASTA $(cut -f 1 $FASTA.fai | grep -Ff $DIR_BASE/partitions/$CHR.contigs.txt -w) >> $DIR_BASE/partitions/1000G-ONT.100x2+4.$CHR.50bkp.fa
        else
            samtools faidx $FASTA $(awk '$2 >= 50000' $FASTA.fai | cut -f 1 | grep -Ff $DIR_BASE/partitions/$CHR.contigs.txt -w) >> $DIR_BASE/partitions/1000G-ONT.100x2+4.$CHR.50bkp.fa
        fi
    done

    bgzip -@ 48 -l 9 $DIR_BASE/partitions/1000G-ONT.100x2+4.$CHR.50bkp.fa
    samtools faidx $DIR_BASE/partitions/1000G-ONT.100x2+4.$CHR.50bkp.fa.gz
    
    rm $DIR_BASE/partitions/$CHR.contigs.txt
done
```
