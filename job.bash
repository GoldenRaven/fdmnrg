rm -fr `find .|grep -v "\<job.bash\>"`

Lambda=2
Gamma=0.02
dir=`pwd`

for p in 0 #0.2 0.4 0.6 0.8
do
	for T in 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11
	do
		mkdir $T
		cd $T
		echo 0.12  > U
		echo -0.04 > Ed_up
		echo -0.04 > Ed_down
		echo $T    > temperature
		echo 2     > Lambda
		echo 0.69  > alpha
		echo 1     > occupation

		cp /home/ligy/NRG/fdmnrg.x .
		cp /home/ligy/NRG/makeinput.cpp .
		cp /home/ligy/NRG/makefreq.cpp .
		cp /home/ligy/NRG/test.submit .
		
		icpc makeinput.cpp -o makeinput.x
		./makeinput.x
		
		icpc makefreq.cpp -o makefreq.x
		./makefreq.x

		cp ~/data_chain/chain_Lambda${Lambda}_Gamma${Gamma}_p${p}.dat chain_band.dat
		cp ~/data_chain/chain_Lambda${Lambda}_Gamma${Gamma}_p${p}.dat chain_total.dat
		job_name=occu-$p
		qsub -N ${job_name} test.submit
		cd $dir
		done
	
done
#cat `find .|grep out`|grep entropy > entr
#awk '{print $2/5.15e-5,"      ",$5/0.69314718}' entr > s
exit 0
