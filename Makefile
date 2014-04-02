OBJECTS = util.o dispe.o dispf.o discon.o dispft.o dispbc.o tests.o main.o
.PHONY: clean

out: main.exe
	./main.exe

main.exe: $(OBJECTS)
	gfortran $(OBJECTS) -o main.exe

%.o : %.F
	gfortran -c $<

clean:
	rm -f $(OBJECTS) main.exe
