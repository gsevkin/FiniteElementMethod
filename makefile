CFLAGS = -Wall -std=c++11 -pedantic

make:
	g++ $(CFLAGS) ./main.cpp ./waveguide.cpp ./CG.cpp ./IPI.cpp -o waveguide

clean:
	rm -rf *.o ksq.txt A.txt M.txt eigenmode.txt waveguide

.PHONY: clean
