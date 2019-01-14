#!/bin/bash
pbsnodes > tmp
# sed -i '1,16d' tmp
echo "###########################################" > tmpp
echo "Total node number:   8" >> tmpp
echo "CPU number per node: 2" >> tmpp
echo "Core number per CPU: 14" >> tmpp
echo "###########################################" >> tmpp
echo "Node   CPU  free_core        state" >> tmpp
i=1
for ((node=3;node<=9;node++))
do
    for sub_node in 0 1
    do
	state=`cat tmp|grep 'state ='|cut -d = -f 2|sed -n "${i}p"`
	# str_=`awk -v i=$node 'BEGIN {printf "00%d", i}'`
	line=`grep -nr "node00${node}-${sub_node}" tmp|cut -d : -f 1`
	line=$[ $line + 4]
	free_core=`sed -n "${line}p" tmp|grep "jobs ="|grep -o '.node'|wc -l`
	free_core=`echo "14-${free_core}"|bc`
	echo $node "    " $sub_node "    " $free_core "            " $state >> tmpp
	# i=$[$i+1]
	((i++))
    done
done
i=1
for ((node=10;node<=10;node++))
do
    for sub_node in 0 1
    do
	state=`cat tmp|grep 'state ='|cut -d = -f 2|sed -n "${i}p"`
	# str_=`awk -v i=$node 'BEGIN {printf "00%d", i}'`
	line=`grep -nr "node0${node}-${sub_node}" tmp|cut -d : -f 1`
	echo $line
	line=$[ $line + 4]
	free_core=`sed -n "${line}p" tmp|grep "jobs ="|grep -o '.node'|wc -l`
	free_core=`echo "14-${free_core}"|bc`
	echo $node "    " $sub_node "    " $free_core "            " $state >> tmpp
	# i=$[$i+1]
	((i++))
    done
done
echo "###########################################" >> tmpp
#
cat tmpp
cat tmpp |grep -v 'exclusive'|grep -v down
rm -f tmp tmpp
