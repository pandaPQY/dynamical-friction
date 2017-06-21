FC = ifort
FFLAG = -O3 -qopenmp 
RFLAG = #-r8

all: main

main:  init.o params_cosmology.o work.o output.o main.o
	$(FC) $(FFLAG) $(RFLAG) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAG) $(RFLAG) -c $< 
clean: 
	rm *.o
	rm *.mod
#	rm main
