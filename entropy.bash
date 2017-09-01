cd /home/ligy/my-tool-box/sort/
make clean                             
make                                   
cd -                                   
dir=`pwd`
#for file in Vg*
#do
#	cd $file
#	pwd
#	Vg=`echo $file |cut -d "_" -f 1 |cut -d g -f 2`
#	h=`echo $file |cut -d "_" -f 2 |cut -d h -f 2`
#	for z in 0*
#	do
#		cd $z
#		cat `find .|grep out` |grep ent |grep -vi job> ent
#		awk '{print $2,"      ",$5/0.69314718}' ent > s
#		cat s|wc -l > raw_num
#		cp /home/ligy/my-tool-box/sort/sort.x .
#		./sort.x
#		cp s ../s$z
#		cd ..
#	done
#	paste s0* > s
#	awk '{print $1,"      ",($2+$4+$6+$8)/4}' s > S
#	cp S ../S_Vg${Vg}_h${h}
#	cd ..
#done
#
for Vg in 0
do
	for h in 0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28
	do
		#echo $file $Vg  $h $T1 $T2 $dT > tmpp
		cd Vg${Vg}_h${h}
		for z in 0*
		do
			cd $z
			cat `find .|grep out` |grep ent |grep -vi job> ent
			awk '{print $2,"      ",$5/0.69314718}' ent > s
			cat s|wc -l > raw_num
			cp /home/ligy/my-tool-box/sort/sort.x .
			./sort.x
			cp s ../s$z
			cd ..
		done
		paste s0* > s
		awk '{print $1,"      ",($2+$4+$6+$8)/4}' s > S
		cp S ../S_Vg${Vg}_h${h}
		cd $dir
	done
done
