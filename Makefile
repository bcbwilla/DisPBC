
OBJECTS = util.o dispe.o dispf.o discon.o dispbc.o main.o
.PHONY: clean

out: main.exe
	./main.exe

main.exe: $(OBJECTS)
	gfortran $(OBJECTS) -o main.exe

%.o : %.F
	gfortran -c $<

clean:
	rm -f $(OBJECTS) main.exe
