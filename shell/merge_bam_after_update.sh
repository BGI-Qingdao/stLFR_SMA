#!/bin/bash

DIR="$( cd "$( dirname "$0"  )" && pwd  )"
helpdoc()
{ cat <<EOF
Usage:

    The shell is used to merge bam after solved multip align in the designated floders

Options:

    -f   The update alignment bam path
    -o   The output file
    -t   The thread number
    -s   The samtools path
EOF
}

while getopts ":t:f:m:o:s:h" opt
do
    case $opt in
        t) threads=$OPTARG;;
        f) floder=$OPTARG;;
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

if [ -z $floder ];
then
	echo "The aligment bam floder path are required"
	exit 1
fi

if [ -z $threads ];
then
	echo "The threads are required"
	exit 1
fi

`find ${floder}/ -type f > name_list.txt`

$samtools merge -f --no-PG --threads $threads ${output}.bam -b name_list.txt
