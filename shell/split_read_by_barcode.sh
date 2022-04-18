#! /bin/bash

DIR="$( cd "$( dirname "$0"  )" && pwd  )"
readfile1=$1
readfile2=$2
ref=$3
max_workers=$4

mkdir barcode_floder
mkdir merge_barcode_1
mkdir merge_barcode_2

echo "==== Start split read by barcode $(date +%F%t%T) ===="
gzip -dc ${readfile1} |  awk -F '@|#|/|\t' '{if(NR%4==1&&NF>1)printf("%s#%s/%s\n",$2,$3,$4) >> "barcode_floder/"$3".txt"}'

i=0
for file in barcode_floder/*;
do
    less $file >> merge_barcode_1/$i"_barcode_list.txt"
	i=$((i=i+1))
	if ((i%max_workers==0))
	then
		i=0
	fi
done

for((j=0;j<max_workers;j++));
do
	awk -F '/' '{printf("%s/2\n",$1)}' merge_barcode_1/$j"_barcode_list.txt" >> merge_barcode_2/$j"_barcode_list.txt"
done

mkdir read_1_floder
mkdir read_2_floder

for((j=0;j<max_workers;j++));
do
	/ldfssz1/ST_OCEAN/USER/xiaogaohong/xiaogaohong/software/anaconda3/bin/seqtk subseq ${readfile1} merge_barcode_1/$j"_barcode_list.txt" > read_1_floder/$j"_1.fq" &
done
wait

for((j=0;j<max_workers;j++));
do
	/ldfssz1/ST_OCEAN/USER/xiaogaohong/xiaogaohong/software/anaconda3/bin/seqtk subseq ${readfile2} merge_barcode_2/$j"_barcode_list.txt" > read_2_floder/$j"_2.fq" &
done
wait
echo "==== END split read $(date +%F%t%T) ===="

mkdir r1_bam
mkdir r2_bam

echo "==== Start to align $(date +%F%t%T) ===="
for((j=0;j<max_workers;j++));
do
/share/app/bwa-0.7.12/bwa mem -a -t 1 $ref read_1_floder/$j"_1.fq"  1>r1_bam/$j".r1.sam" 2> r1_bam/aln.err &
done
wait

for((j=0;j<max_workers;j++));
do
/share/app/samtools/1.11/bin/samtools sort -@ 1 -o r1_bam/$j".r1.sorted.bam" -O bam -T $j"_tmp" r1_bam/$j".r1.sam" &
done
wait

for((j=0;j<max_workers;j++));
do
/share/app/bwa-0.7.12/bwa mem -a -t 1 $ref read_2_floder/$j"_2.fq"  1>r2_bam/$j".r2.sam" 2> r2_bam/aln.err &
done
wait

for((j=0;j<max_workers;j++));
do
/share/app/samtools/1.11/bin/samtools sort -@ 1 -o r2_bam/$j".r2.sorted.bam" -O bam -T $j"_tmp" r2_bam/$j".r2.sam" &
done
wait

for((j=0;j<max_workers;j++));
do
rm -rf r1_bam/$j".r1.sam"
rm -rf r2_bam/$j".r2.sam"
done
wait
echo "==== End aligning $(date +%F%t%T) ===="

mkdir merge_bam
cd merge_bam
echo "==== Start merge bam $(date +%F%t%T) ===="
for((j=0;j<max_workers;j++));
do
	bash $DIR/merge_bam.sh ../r1_bam/$j".r1.sorted.bam" ../r2_bam/$j".r2.sorted.bam" $j"_ref" 1 &
done
wait
echo "==== End merge $(date +%F%t%T) ===="
