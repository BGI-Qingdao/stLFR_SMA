#! /bin/bash

DIR="$( cd "$( dirname "$0"  )" && pwd  )"
helpdoc()
{ cat <<EOF

    The shell is used to extract read in name list to the file designated

    -r   The readfile path
    -n	 The name list floder
    -o	 The output floder designated
    -t   The thread number
    -m	 The workers number
    -i	 The index of read [1 or 2]

EOF
}

while getopts ":r:t:o:m:n:i:h" opt
do
	case $opt in
		t) threads=$OPTARG;;
		m) workers=$OPTARG;;
		r) readfile=$OPTARG;;
		n) namelist=$OPTARG;;
		o) output=$OPTARG;;
		i) index=$OPTARG;;
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

if [ -z $readfile ];
then
	echo "The readfile are required"
	exit 1
fi

if [ -z $namelist ];
then
	echo "The read name list are required"
	exit 1
fi

if [ -z $output ];
then 
	echo "The output floder are required"
	exit 1
fi

if [ -z $index ];
then 
	echo "The index of read are required 1 or 2 "
	exit 1
fi

#echo $DIR
#echo "==== Start extract read$index by barcode $(date +%F%t%T) ===="

j=0
while (($j<workers))
do
	for ((i=0;i<threads;i++))
	do
		$DIR/tools/bin/seqtk subseq ${readfile} $namelist/$j"_read_list.txt" >$output/${j}_${index}.fq &
		j=$(($j+1))
		if (($j==workers))
		then 
			break
		fi
	done
	wait
done
wait

#echo "==== END extract read$index $(date +%F%t%T) ===="
