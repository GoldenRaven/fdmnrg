#rm -fr `find .|grep -v "\<job.bash\>"`

job_id=2
U=0.2
Ed=-0.1
Lambda=1.8
#alpha=0.5
alpha=0.4
num_kept=512
smear=0
unsmear=1
dim_imp=4
dim_dot=4
Beta_bar=0.6
Q=0
Q_Sz=0
N_up_N_down=1
occupation=1
imp_dos=1
stm_dos=0
tc=0.5
td=0.5
N_max=0
eig=1
smooth='newsc' #'wvd'

Nz=4
dir=`pwd`
Gamma0=0.1
Vg=0
t=7
for h in 0.14
do
    for zeeman in -0.0176015
    do
	echo $h  $zeeman
	mkdir h${h}_zeeman${zeeman}_T=1e-${t} #_Ns=1024
	cd h${h}_zeeman${zeeman}_T=1e-${t} #_Ns=1024
	dir1=`pwd`
	for temperature in 1e-$t
	do
	    omega0=$temperature
	    Omega=$temperature
			#rm -fr $temperature
	    mkdir $temperature
	    cd $temperature
	    for z in 0 0.25 0.5 0.75
			# for z in 0.75
	    do
				#rm -fr $z
		mkdir $z
		cd $z
		cp /home/ligy/NRG/fdmnrg.x .

		echo $job_id      >  job_id

		echo $U           >> input_total
		Ed_up=`echo "${zeeman}+${Ed}"|bc -l`
		Ed_down=`echo "-1*${zeeman}+${Ed}"|bc -l`
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
		echo $imp_dos     >> input_total
		echo $stm_dos     >> input_total
		echo $tc          >> input_total
		echo $td          >> input_total
		echo $N_max       >> input_total
		echo $eig         >> input_total
		echo $Omega       >> input_total

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

				# cp /home/ligy/NRG/makefreq.cpp .
		cp  $dir/freqency .
		cp $dir/test.submit .

				# icpc makefreq.cpp -o makefreq.x
				# ./makefreq.x

		cp ~/data_chain_fitting/chain_Lambda${Lambda}_r1_k${Gamma0}_Vg${Vg}_h${h}_z${z}.dat chain_band.dat
		cp ~/data_chain_fitting/chain_Lambda${Lambda}_r1_k${Gamma0}_Vg${Vg}_h${h}_z${z}.dat chain_total.dat
		job_name=zero_${zeeman}_$z
		qsub -N ${job_name} test.submit
		cd ..
	    done
	    cd ..
	done
	cd ..
    done
    cd ..
done
exit 0
