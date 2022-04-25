CXX = g++ -pthread -o 
SRC = main.cc
TAR = $(SRC:.cc=.x)
DEP = Laplace.cc Laplace.h

all: $(TAR)

$(TAR):$(SRC) $(DEP)
	$(CXX) $(TAR) $<

run:
	./$(TAR)

plot:
	python plot.py

clean:
	$(RM) *.x *.txt
