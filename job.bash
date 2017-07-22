#改正了之前实际上z全为0的错误.
#rm -fr `find .|grep -v "\<job.bash\>"`

job_id=1
U=0.2
Ed_up=-0.1
Ed_down=-0.1
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
N_max=0

Gamma0=0.01
Nz=4
dir=`pwd`

cp /home/ligy/NRG/entropy.bash .

for Vg in 0.3
do
	for h in 0.17 0.22
	do
		rm -fr Vg${Vg}_h$h
		mkdir Vg${Vg}_h$h
		cd Vg${Vg}_h$h
		dir1=`pwd`
		for z in 0 0.25 0.5 0.75
		do
			rm -fr $z
			mkdir $z
			cd $z
			#for temperature in 5.15e-2 3e-2 1e-2 7e-3 5.15e-3 2e-3 9e-4 5.15e-4 2e-4 9e-5 5.15e-5 1e-5 8e-6 5.15e-6 1e-6 8e-7 5.15e-7 1e-7 8e-8 5.15e-8
			for ((i=30;i<50;i++))
			do
				#temperature=echo "(0.8^$i)*1000"|bc -l
				dir2=`pwd`
				temperature=`echo "0.8 $i 10"| awk '{ printf "%0.15f\n" ,$1^$2*$3}'`
				omega0=$temperature
				rm -fr $temperature
				mkdir $temperature
				cd $temperature
		
				cp /home/ligy/NRG/fdmnrg.x .
		
				echo $job_id      >  job_id
		
				echo $U           >> input_total
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
				echo $N_max       >> input_total
		
				echo $temperature >> input_band
				echo $Lambda      >> input_band
				echo $num_kept    >> input_band 
				echo $dim_imp     >> input_band 
				echo $dim_dot     >> input_band 
				echo $Beta_bar    >> input_band 
				echo $Q           >> input_band 
				echo $Q_Sz        >> input_band 
				echo $N_up_N_down >> input_band 
				echo $N_max       >> input_band
		
				cp /home/ligy/NRG/makefreq.cpp .
				cp /home/ligy/NRG/test.submit .
				
				icpc makefreq.cpp -o makefreq.x
				./makefreq.x
		
				cp ~/data_chain/chain_Lambda${Lambda}_r1_k${Gamma0}_Vg${Vg}_h${h}_z${z}.dat chain_band.dat
				cp ~/data_chain/chain_Lambda${Lambda}_r1_k${Gamma0}_Vg${Vg}_h${h}_z${z}.dat chain_total.dat
				job_name=$h-$z-$i
				qsub -N ${job_name} test.submit
				cd $dir2
			done
			cd $dir1
		done
		cd $dir
	done
done
exit 0
