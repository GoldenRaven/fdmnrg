#rm -fr `find .|grep -v "\<job.bash\>"`

Lambda=10
Gamma=0.001
U=0.01
Ed_up=-0.005
Ed_down=-0.005
num_kept=1024
alpha=0.69
occupation=1
Nz=4
dir=`pwd`
cp /home/ligy/NRG/entropy.bash .

for z in 0.125 0.375 0.5 0.75
do
	rm -fr $z
	mkdir $z
	cd $z
	#for T in 5.15e-2 3e-2 1e-2 7e-3 5.15e-3 2e-3 9e-4 5.15e-4 2e-4 9e-5 5.15e-5 1e-5 8e-6 5.15e-6 1e-6 8e-7 5.15e-7 1e-7 8e-8 5.15e-8
	for ((i=0;i<110;i++))
	do
		#T=echo "(0.8^$i)*1000"|bc -l
		dir2=`pwd`
		T=`echo "0.8 $i 1000"| awk '{ printf "%0.15f\n" ,$1^$2*$3}'`
		rm -fr $T
		mkdir $T
		cd $T
		echo $U          > U
		echo $Ed_up      > Ed_up
		echo $Ed_down    > Ed_down
		echo $T          > temperature
		echo $Lambda     > Lambda
		echo $num_kept   > num_kept
		echo $alpha      > alpha
		echo $occupation > occupation

		cp /home/ligy/NRG/fdmnrg.x .
		cp /home/ligy/NRG/makeinput.cpp .
		cp /home/ligy/NRG/makefreq.cpp .
		cp /home/ligy/NRG/test.submit .
		
		icpc makeinput.cpp -o makeinput.x
		./makeinput.x
		
		icpc makefreq.cpp -o makefreq.x
		./makefreq.x

		cp ~/data_chain/chain_Lambda${Lambda}_Gamma${Gamma}_z${z}.dat chain_band.dat
		cp ~/data_chain/chain_Lambda${Lambda}_Gamma${Gamma}_z${z}.dat chain_total.dat
		job_name=$z-$T
		qsub -N ${job_name} test.submit
		cd $dir2
	done
	cd $dir
done
exit 0
