#!/usr/bin/bash

args=("$@")
bam=${args[0]}
id=${args[1]}
path=${args[2]}
samtools=${args[3]}

genome=$($samtools depth $bam |  awk '{sum+=$3} END { print sum/NR}')

echo "Ratio Genome = "$genome > $path/genome_coverage_$id.txt


regions=$(for i in `seq 1 1000 15073001`; do j=$(($i+1000));echo -e "I:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_I_ratio_depth_1kb.txt"

regions=$(for i in `seq 1 1000 15280001`; do j=$(($i+1000));echo -e "II:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_II_ratio_depth_1kb.txt"

regions=$(for i in `seq 1 1000 13790001`; do j=$(($i+1000));echo -e "III:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_III_ratio_depth_1kb.txt"

regions=$(for i in `seq 1 1000 17494001`; do j=$(($i+1000));echo -e "IV:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_IV_ratio_depth_1kb.txt"

regions=$(for i in `seq 1 1000 20925001`; do j=$(($i+1000));echo -e "V:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_V_ratio_depth_1kb.txt"

regions=$(for i in `seq 1 1000 17719001`; do j=$(($i+1000));echo -e "X:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_X_ratio_depth_1kb.txt"

chrs="I II III IV V X"
for chr in $chrs; do file=$id"/"$id"_chr"$chr"_ratio_depth_1kb.txt"; output=$id"_ratio_1kb.tsv"; sed 's/:/\t/g' $file | sed 's/-/\t/g' | sed -z 's/\nRatio Genome =  /\t/g' | awk '{if ($4=="") print "chr"$1"\t"$2"\t"$3"\t0"; else print "chr"$0}' >> $output; done; done

