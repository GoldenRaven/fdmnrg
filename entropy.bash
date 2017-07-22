cd /home/ligy/my-tool-box/sort/
make clean                             
make                                   
cd -                                   
dir=`pwd`
for file in 0*
do
	echo $file
	cd $file
	pwd
	cat `find .|grep out` |grep ent |grep -vi job> ent
	awk '{print $2,"      ",$5/0.69314718}' ent > s
	cat s|wc -l > raw_num
	cp /home/ligy/my-tool-box/sort/sort.x .
	./sort.x
	cp s $dir/s$file
	cd $dir
done
paste s0* > s
awk '{print $1,"      ",($2+$4+$6+$8)/4}' s > S
#
rm -f gnu.txt
echo "set term pdf font 'Times New Roman,6'" >> gnu.txt                                                           
echo "set output 'S.pdf'" >> gnu.txt
echo "set grid" >> gnu.txt
echo "set title 'U=0.01,Ed=-0.005,Gamma=0.001,D=1,Lambda=10,Nz=4'" >> gnu.txt
echo "set ylabel 'S/ln2'" >> gnu.txt
echo "set xlabel 'T'" >> gnu.txt
echo "#set yrange [0:2.5]" >> gnu.txt
echo "set logscale x" >> gnu.txt
echo "#plot '~/Documents/fig30.png.dat' w lp t 'PRB-86-75153'  ps 0.2 pt 7,'S' w p t 'mine' ps 0.2" >> gnu.txt
echo "plot 'S' w p t 'mine' ps 0.2" >> gnu.txt
echo "set output" >> gnu.txt
echo "q" >> gnu.txt
#                 
cat gnu.txt|gnuplot
