cd /home/ligy/my-tool-box/sort/
make clean                             
make                                   
cd -                                   
dir=`pwd`
for file in Vg*
do
	cd $file
	pwd
	Vg=`echo $file |cut -d "_" -f 1 |cut -d g -f 2`
	h=`echo $file |cut -d "_" -f 2 |cut -d h -f 2`
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
	cd ..
done
