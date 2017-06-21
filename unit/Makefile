FC = ifort

FFLAGS = -r8

all : output

output: init.o unit.o
	$(FC) $(FFLAGS) -o $@ $^

init.o: init.f90
	$(FC) $(FFLAGS) -c $<

unit.o: unit.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm *.o
	rm *.mod
	rm output
