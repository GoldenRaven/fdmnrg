rm -f `find .|grep -v "\<job.bash\>"`
cp /home/ligy/NRG/fdmnrg.x .
cp /home/ligy/NRG/makeinput.cpp .
cp /home/ligy/NRG/makefreq.cpp .
cp /home/ligy/NRG/test.submit .

echo 0.5   > U
echo -0.25 > Ed
echo 1e-7  > temperature
echo 1.8   > Lambda
echo 0.588 > alpha

Lambda=1.8
k=0.1
Vg=0.1
h=0.02

icpc makeinput.cpp -o makeinput.x
./makeinput.x

icpc makefreq.cpp -o makefreq.x
./makefreq.x

cp ~/chain_Lambda${Lambda}_r1_k${k}_Vg${Vg}_h$h.dat chain.dat
job_name=0.5_test
qsub -N ${job_name} test.submit
exit 0
