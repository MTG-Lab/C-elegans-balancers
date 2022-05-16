#!/usr/bin/bash

args=("$@")
bam=${args[0]}
id=${args[1]}
path=${args[2]}
cluster=${args[3]}


if [ $cluster == "arc" ]; then
	samtools="";
fi

if [ $cluster == "marc" ]; then
	samtools="/project/M-mtgraovac182840/tools/samtools-1.3.1/samtools";
fi

if [ $cluster == "gpcc-node01" ]; then
	source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh; 
	source /cvmfs/soft.computecanada.ca/config/profile/bash.sh; 
	module load StdEnv/2020; 
	module load samtools/1.11;
	samtools="samtools";
fi

if [ $cluster == "gpcc-node02" ]; then
	source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh; 
	source /cvmfs/soft.computecanada.ca/config/profile/bash.sh; 
	module load StdEnv/2020; 
	module load samtools/1.11;
	samtools="samtools";
fi


genome=$($samtools depth $bam |  awk '{sum+=$3} END { print sum/NR}')

echo "Ratio Genome = "$genome > $path/genome_coverage_$id.txt


regions=$(for i in `seq 1 10000 245300000`; do j=$(($i+10000));echo -e "1:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr1_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 243400000`; do j=$(($i+10000));echo -e "2:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr2_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 199500000`; do j=$(($i+10000));echo -e "3:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr3_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 191700000`; do j=$(($i+10000));echo -e "4:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr4_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 181000000`; do j=$(($i+10000));echo -e "5:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr5_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 170800000`; do j=$(($i+10000));echo -e "6:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr6_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 158500000`; do j=$(($i+10000));echo -e "7:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr7_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 146000000`; do j=$(($i+10000));echo -e "8:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr8_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 134600000`; do j=$(($i+10000));echo -e "9:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr9_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 135500000`; do j=$(($i+10000));echo -e "10:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr10_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 135000000`; do j=$(($i+10000));echo -e "11:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr11_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 133500000`; do j=$(($i+10000));echo -e "12:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr12_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 114200000`; do j=$(($i+10000));echo -e "13:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr13_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 105400000`; do j=$(($i+10000));echo -e "14:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr14_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 100200000`; do j=$(($i+10000));echo -e "15:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr15_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 90000000`; do j=$(($i+10000));echo -e "16:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr16_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 81700000`; do j=$(($i+10000));echo -e "17:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr17_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 77800000`; do j=$(($i+10000));echo -e "18:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr18_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 63800000`; do j=$(($i+10000));echo -e "19:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr19_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 63700000`; do j=$(($i+10000));echo -e "20:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr20_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 47000000`; do j=$(($i+10000));echo -e "21:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr21_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 49500000`; do j=$(($i+10000));echo -e "22:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chr22_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 156000000`; do j=$(($i+10000));echo -e "X:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chrX_ratio_depth_10kb.txt"

regions=$(for i in `seq 1 10000 51000000`; do j=$(($i+10000));echo -e "Y:$i-$j"; done)
for r in $regions; do echo $r; $samtools depth -r $r $bam | awk -v gen="$genome" '{sum+=$3} END { print "Ratio Genome = ",sum/NR/gen}' ; done > $path/$id"_chrY_ratio_depth_10kb.txt"



