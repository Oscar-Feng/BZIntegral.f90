FCOM=gfortran
EXE=temp

SOURCES=temp.f90
MODS=lintetra.f90 quadtetra.f90 tetrasupply.f90 splitmesh.f90 bzint3d.f90

OBJS=$(SOURCES:.f90=.o)
MODOBJS=$(MODS:.f90=.o)

all: $(MODOBJS) $(OBJS)
	$(FCOM) $(OBJS) $(MODOBJS) -o $(EXE) -fopenmp -O3

$(MODOBJS): %.o: %.f90
	$(FCOM) -c -o $@ $*.f90 -fopenmp -O3

$(OBJS): %.o: %.f90
	$(FCOM) -c -o $@ $*.f90 -fopenmp -O3

.PHONY: clean
clean:
	rm *.o *.mod
