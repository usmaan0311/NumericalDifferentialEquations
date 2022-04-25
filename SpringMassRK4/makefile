CXX=g++ -Wall -O3 -o
SRC=main.cc
TAR:=$(SRC:.cc=.x)

all: $(TAR)

$(TAR) : $(SRC)
	$(CXX) $(TAR) $<

run:
	./$(TAR)

plot:
	python plot.py

clean:
	$(RM) *.dat *.x

