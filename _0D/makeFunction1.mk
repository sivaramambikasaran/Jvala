CC	=g++-5 -std=c++11
CFLAGS	=-c -Wall -Ofast -ffast-math -ffinite-math-only -I header/
LDFLAGS	=
SOURCES	=./Anukanala-master/src/RK_NonAdaptive.cpp ./Anukanala-master/src/RK_Adaptive.cpp ./Anukanala-master/src/Adambashforth.cpp ./Anukanala-master/src/Euler_Implicit.cpp ./Anukanala-master/src/LeapFrog.cpp ./src/Isobaric1.cpp ./src/Isochoric1.cpp ./src/Properties1.cpp ./examples/test.cpp
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./exec/test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $(DIM) $< -o $@

clean:
	rm -rf *.out ./examples/*.o ./src/*.o ./exec/*

tar:
	tar -zcvf Anukalana.tar.gz ./makefile.mk ./exec ./src ./header ./examples ./README.md ./LICENSE.md
