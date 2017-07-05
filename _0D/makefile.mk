CC	=g++-5 -std=c++11
CFLAGS	=-c -g -Wall -Ofast -ffast-math -ffinite-math-only -I header/ 
LDFLAGS	=
SOURCES	= ./src/Isobaric1.cpp ./src/Isochoric1.cpp ./src/Properties1.cpp ./examples/test.cpp
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./exec/test

subsystem:
	$(MAKE) -C ./../Anukalana
all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $(DIM) $< -o $@

clean:
	rm -rf *.out ./examples/*.o ./src/*.o ./exec/* 

tar:
	tar -zcvf Anukalana.tar.gz ./makefile.mk ./exec ./src ./header ./examples ./README.md ./LICENSE.md
