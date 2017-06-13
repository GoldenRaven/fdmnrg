#rm -fr `find .|grep -v "\<job.bash\>"`

U=0.01
Ed_up=-0.005
Ed_down=-0.005
Lambda=10
alpha=0.69
num_kept=512
smear=1
unsmear=1
dim_imp=4
dim_dot=4
Beta_bar=0.6
Q=0
Q_Sz=0
N_up_N_down=1
occupation=1
dos=1

Gamma=0.001
Nz=4
dir=`pwd`

cp /home/ligy/NRG/entropy.bash .

for z in 0.125 0.375 0.5 0.75
do
	rm -fr $z
	mkdir $z
	cd $z
	#for temperature in 5.15e-2 3e-2 1e-2 7e-3 5.15e-3 2e-3 9e-4 5.15e-4 2e-4 9e-5 5.15e-5 1e-5 8e-6 5.15e-6 1e-6 8e-7 5.15e-7 1e-7 8e-8 5.15e-8
	for ((i=0;i<110;i++))
	do
		#temperature=echo "(0.8^$i)*1000"|bc -l
		dir2=`pwd`
		temperature=`echo "0.8 $i 1000"| awk '{ printf "%0.15f\n" ,$1^$2*$3}'`
		omega0=$temperature
		rm -fr $temperature
		mkdir $temperature
		cd $temperature

		cp /home/ligy/NRG/fdmnrg.x .

		echo $U           > input_total
		echo $Ed_up       >> input_total
		echo $Ed_down     >> input_total
		echo $temperature >> input_total
		echo $Lambda      >> input_total
		echo $alpha       >> input_total
		echo $num_kept    >> input_total 
		echo $smear       >> input_total 
		echo $unsmear     >> input_total 
		echo $omega0      >> input_total 
		echo $dim_imp     >> input_total 
		echo $dim_dot     >> input_total 
		echo $Beta_bar    >> input_total 
		echo $Q           >> input_total 
		echo $Q_Sz        >> input_total 
		echo $N_up_N_down >> input_total 
		echo $occupation  >> input_total 
		echo $dos         >> input_total 

		echo $temperature >> input_band
		echo $Lambda      >> input_band
		echo $num_kept    >> input_band 
		echo $dim_imp     >> input_band 
		echo $dim_dot     >> input_band 
		echo $Beta_bar    >> input_band 
		echo $Q           >> input_band 
		echo $Q_Sz        >> input_band 
		echo $N_up_N_down >> input_band 

		cp /home/ligy/NRG/makefreq.cpp .
		cp /home/ligy/NRG/test.submit .
		
		icpc makefreq.cpp -o makefreq.x
		./makefreq.x

		cp ~/data_chain/chain_Lambda${Lambda}_Gamma${Gamma}_z${z}.dat chain_band.dat
		cp ~/data_chain/chain_Lambda${Lambda}_Gamma${Gamma}_z${z}.dat chain_total.dat
		job_name=$z-$temperature
		qsub -N ${job_name} test.submit
		cd $dir2
	done
	cd $dir
done
exit 0
