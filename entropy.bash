dir=`pwd`
for file in 0*
do
	echo $file
	cd $file
	pwd
	cat `find .|grep out` |grep ent > ent$file
	awk '{print $2/5.15e-5,"      ",$5/0.69314718}' ent$file > s$file
	cp s$file $dir
	cd $dir
done
paste s* > s
awk '{print $1,"      ",($2+$4+$6+$8)/4}' s > S
#
rm -f gnu.txt
echo "set term pdf font 'Times New Roman,6'" >> gnu.txt                                                           
echo "set output 'S.pdf'" >> gnu.txt
echo "set grid" >> gnu.txt
echo "set title 'U=0.01,Ed=-0.005,Gamma=0.001,D=1,Lambda=5,Nz=4,Tk=5.15e-5'" >> gnu.txt
echo "set ylabel 'S/ln2'" >> gnu.txt
echo "set xlabel 'T/T_k'" >> gnu.txt
echo "#set yrange [0:2.5]" >> gnu.txt
echo "set logscale x" >> gnu.txt
echo "plot '~/Documents/fig30.png.dat' w lp t 'PRB-86-75153'  ps 0.2 pt 7,'S' w p t 'mine' ps 0.2" >> gnu.txt
echo "set output" >> gnu.txt
echo "q" >> gnu.txt
#                 
cat gnu.txt|gnuplot
