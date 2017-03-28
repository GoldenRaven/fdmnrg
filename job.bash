#rm -f `find .|grep -v "\<job.bash\>"`
cp /home/ligy/NRG/fdmnrg.x .
cp /home/ligy/NRG/makeinput.cpp .
#cp /home/ligy/NRG/makefreq.cpp .
cp /home/ligy/NRG/test.submit .

echo 0     > findpeak 
echo -0.15 > Ed
echo 1e-8  > temperature
echo 1e-6  > x_min   # to be modified!
echo 1e-1  > x_max   # to be modified!
echo 1e-9  > x_error
echo 0.3   > U
echo 1.8   > Lambda
echo 0.588 > alpha
echo 512 > num_kept 

Lambda=1.8
k=0.1
Vg=0.3
h=
Nz=1

icpc makeinput.cpp -o makeinput.x
./makeinput.x

icpc makefreq.cpp -o makefreq.x
./makefreq.x

cp ~/chain_Lambda${Lambda}_r1_k${k}_Vg${Vg}_h${h}.dat chain.dat
#cp ~/chain_Lambda${Lambda}_r1_k${k}_Vg${Vg}_h${h}_Nz${Nz}.dat chain.dat
job_name=10Tk
qsub -N ${job_name} test.submit
exit 0
