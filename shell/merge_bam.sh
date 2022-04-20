#!/bin/bash
helpdoc()
{ cat <<EOF
Usage:

    The shell is used to merge bam by read name in the designated floders

Options:

    -1   The read_1 alignment bam file path
    -2   The read_2 alignment bam file path
    -o   The output floder designated
    -m   The worker number
    -s   The samtools path
EOF
}

while getopts ":t:1:2:m:o:s:h" opt
do
    case $opt in
		m) worker=$OPTARG;;
        1) bam_floder_1=$OPTARG;;
        2) bam_floder_2=$OPTARG;;
        o) output=$OPTARG;;
        s) samtools=$OPTARG;;
        h|help) helpdoc
        exit 1;;
        ?) echo "$OPTARG Unknown parameter"
        exit 1;;
    esac
done

if [ $# = 0 ];
then
    helpdoc
    exit 1
fi

if [ -z $bam_floder_1 ];
then
	echo "The read_1 aligment bam file path are required"
	exit 1
fi

if [ -z $bam_floder_2 ];
then
	echo "The read_2 aligment bam file path are required"
	exit 1
fi

$samtools merge --no-PG -r --threads 1 ${output}/${worker}.bam ${bam_floder_1}/${worker}.r1.sorted.bam ${bam_floder_2}/${worker}.r2.sorted.bam
$samtools sort -@ 1 -n -o ${output}/${worker}.sorted.bam -O bam -T ${output}/${worker}_tmp ${output}/${worker}.bam
rm -f ${output}/${worker}.bam
