FFLAGS = -Wall -Wextra -g -fdefault-real-8 -fdefault-double-8 -fbacktrace -O3

OBJECTS = LinAl.o Driver_LinAl.o

MODULES = LinAl.mod

.PHONY: clean

output.txt: LinAl.exe
	./LinAl.exe > output.txt

LinAl.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o LinAl.exe -llapack

%.o: %.f90
	gfortran $(FFLAGS) -c $<

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(MODULES) LinAl.exe
