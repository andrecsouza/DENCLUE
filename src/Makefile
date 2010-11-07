CPP=g++-4.0
INCLUDE=-I../include/
FLAGS=-Wall -ggdb #-O2 -ffast-math
OBJECTS= dataset.o hypercube.o hyperspace.o denclue_functions.o denclue.o
DEFINE=
LIBS=#-lefence
EXE=denclue
.SUFFIXES : .cpp .o .h

all: $(EXE)

.cpp.o: %.cpp %.h Makefile
	$(CPP) $(FLAGS) $(INCLUDE) $(DEFINE) -c $< -o $@

$(EXE): $(OBJECTS) Makefile
	$(CPP) $(FLAGS) $(LIBS) $(OBJECTS) -o $(EXE)

clean:
	rm -f core $(OBJECTS) *~ $(EXE) out.txt

run: $(OBJECTS) $(EXE) Makefile
	./$(EXE) -d 2 -s 5 -x 2 -i in.txt -o out.txt 2>&1

#./$(EXE) -d 2 -s 0.5 -x 1 -i in.txt -o out.txt 2>&1


