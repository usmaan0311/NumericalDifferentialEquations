CXX = g++ -Wall -O3 -o
SRC = main.cc
IN = Input.txt
PLT = plot.py
TAR := $(SRC:.cc=.x)

all : $(TAR)

$(TAR) : $(SRC)

	$(CXX) $(TAR) $<
vi:
	vim $(SRC)

run:
	./$(TAR)

in:
	vim $(IN)

plot:
	python $(PLT)

clean:
	$(RM) *.x *.dat
