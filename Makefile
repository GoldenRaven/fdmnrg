#.SILENT:
cc=icpc
name=fdmnrg.x 
#name=work_dir/fdmnrg.x 
#objects=main.o setup.o genoutput.o iterative_dia.o dos.o date_time.o deallocate.o func_wn.o
#objects=main.o setup.o genoutput.o iterative_dia.o dos2.o date_time.o deallocate.o func_wn.o
objects=main.o setup.o genoutput.o iterative_dia.o dos3.o date_time.o deallocate.o func_wn.o

#CPPFLAGS=-O0 -g
#CPPFLAGS=-O3
MKL_FLAGS=-mkl=sequential
OPENMP_FLAGS=-qopenmp

$(name): $(objects)
	$(cc) $(CPPFLAGS) $(MKL_FLAGS) $(OPENMP_FLAGS) -o $(name) $(objects) 

main.o: main.cpp 
	$(cc)  $(CPPFLAGS) -c main.cpp

setup.o: setup.cpp setup.h
	$(cc)  $(CPPFLAGS) -c setup.cpp

genoutput.o: genoutput.cpp setup.h
	$(cc)  $(CPPFLAGS) -c genoutput.cpp

date_time.o: date_time.cpp
	$(cc)  $(CPPFLAGS) -c date_time.cpp

iterative_dia.o: iterative_dia.cpp setup.h
	$(cc)  $(CPPFLAGS)  $(MKL_FLAGS) $(OPENMP_FLAGS) -c iterative_dia.cpp 

func_wn.o: func_wn.cpp setup.h
	$(cc)  $(CPPFLAGS) $(OPENMP_FLAGS) -c func_wn.cpp

dos.o: dos.cpp setup.h
	$(cc)  $(CPPFLAGS) $(OPENMP_FLAGS) -c dos.cpp

dos2.o: dos2.cpp setup.h
	$(cc)  $(CPPFLAGS) $(OPENMP_FLAGS) -c dos2.cpp

dos3.o: dos3.cpp setup.h
	$(cc)  $(CPPFLAGS) $(OPENMP_FLAGS) -c dos3.cpp

deallocate.o: deallocate.cpp setup.h
	$(cc)  $(CPPFLAGS)  -c deallocate.cpp

clean:
	rm $(name) $(objects)

