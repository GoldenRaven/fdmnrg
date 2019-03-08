#!/bin/bash
pbsnodes > tmp
# sed -i '1,16d' tmp
curr_time=`date`
echo "###########################################" > tmpp
echo 'current time:' $curr_time >> tmpp
echo "Total node number:   8" >> tmpp
echo "CPU number per node: 2" >> tmpp
echo "Core number per CPU: 14" >> tmpp
echo "###########################################" >> tmpp
echo "Node   CPU  free_core        state" >> tmpp
i=1
for ((node=3;node<=10;node++))
do
    for sub_node in 0 1
    do
	state=`cat tmp|grep 'state ='|cut -d = -f 2|sed -n "${i}p"`
	# str_=`awk -v i=$node 'BEGIN {printf "00%d", i}'`
	if [ $node == 10 ]
	then
	    line=`grep -nr "node0${node}-${sub_node}" tmp|cut -d : -f 1`
	else
	    line=`grep -nr "node00${node}-${sub_node}" tmp|cut -d : -f 1`
	fi
	line=$[ $line + 4]
	free_core=`sed -n "${line}p" tmp|grep "jobs ="|grep -o '.node'|wc -l`
	free_core=`echo "14-${free_core}"|bc`
	echo $node $sub_node $free_core $state |awk '{printf "%-8d%-8d%-13d%-8s\n", $1,$2,$3,$4}' >> tmpp
	# i=$[$i+1]
	((i++))
    done
done
cat tmpp
# cat tmpp |grep -v 'exclusive'|grep -v down
rm -f tmp tmpp
