#!/bin/bash
helpdoc()
{ cat <<EOF
Usage:
    The shell is use bwa mem align single_read in the floer are gived and output to designated floder

Option:
    -r   The reference file path
    -f   The read file path
    -o   The output floder designated
    -m   The worker number
    -b   The BWA path
    -s   The samtools path
    -i   The index of read [1 or 2]

EOF
}

while getopts ":r:t:m:o:f:b:s:i:h" opt
do
    case $opt in
        m) worker=$OPTARG;;
        r) reference=$OPTARG;;
        f) floder=$OPTARG;;
        o) output=$OPTARG;;       
        i) index=$OPTARG;;
        b) bwa=$OPTARG;;
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

if [ -z $reference ];
then
    echo "The reference are required"
    exit 1
fi

$bwa mem -R '@RG\tID:'${worker}.r${index}.sorted'' -a -t 1 $reference $floder 1>${output}/${worker}.r${index}.sam 2>${output}/${worker}_aln.err
$samtools sort -@ 1 -o ${output}/${worker}.r${index}.sorted.bam -O bam -T ${output}/${worker}_tmp ${output}/${worker}.r${index}.sam
rm -rf ${output}/${worker}.r${index}.sam
