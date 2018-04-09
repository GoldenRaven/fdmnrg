#!/bin/bash
pbsnodes > tmp
echo "###########################################" > tmpp
echo "Total node number:   8" >> tmpp
echo "CPU number per node: 2" >> tmpp
echo "Core number per CPU: 14" >> tmpp
echo "###########################################" >> tmpp
echo "CPU    free_core        state" >> tmpp
i=1
for ((node=3;node<=10;node++))
do
    for sub_node in 0 1
    do
	state=`cat tmp|grep 'state ='|cut -d = -f 2|sed -n "${i}p"`
	free_core=`cat tmp|grep 'jobs ='|sed -n "${i}p"|grep -o '.node'|wc -l`
	free_core=`echo "14-${free_core}"|bc`
	# exit 1
	echo $i "    " $free_core "            " $state >> tmpp
	# i=$[$i+1]
	((i++))
    done
done
echo "###########################################" >> tmpp
#
cat tmpp |grep -v 'exclusive'|grep -v down
rm -f tmp tmpp