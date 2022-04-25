CXX = g++ -std=c++17 -O3 -o
SRC = main.cc
DEP = stiff.h stiff.cc
TAR := $(SRC:.cc=.x)

all : $(TAR)

$(TAR) : $(SRC) $(DEP)
	$(CXX) $(TAR) $<

vi:
	vim $(SRC)

run:
	./$(TAR)

plot:
	python plot.py

clean:
	$(RM) *.dat *.x
