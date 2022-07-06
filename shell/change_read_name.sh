#! /bin/bash

helpdoc()
{ cat <<EOF

Usage:

    The shell is used to change read name in name list floder

Options:

    -n   The name list floder
    -t   The thread number
    -m   The worker number
EOF
}

while getopts ":t:n:m:h" opt
do
    case $opt in
        t) threads=$OPTARG;;
		m) workers=$OPTARG;;
        n) namelist=$OPTARG;;
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

if [ -z $namelist ];
then
    echo "The name list floder are required"
    exit 1
fi

j=0
while (($j<workers))
do
	for((i=0;i<threads;i++));
	do
		sed -i "s/\/1/\/2/g" ${namelist}/${j}_read_list.txt &
		j=$(($j+1))
		if (($j==workers))
		then
			break
		fi
	done
	wait
done
wait
