CFLAGS = -std=c99

default: MWT_Calculator

main.o: main.c

read_input.o: read_input.c read_input.h

MWT_Calculator: main.o read_input.o
	$(CC) -lm -o MWT_Calculator main.o read_input.o

clean:
	rm -f main.o read_input.o MWT_Calculator

tar: MWT_Calculator.tgz

MWT_Calculator.tgz: Makefile *.h *.c
	mkdir .MWT_Calculator_tmp && cd .MWT_Calculator_tmp &&\
	mkdir MWT_Calculator && cp ../Makefile ../*.h ../*.c ../parameter_input ../readme.txt ./MWT_Calculator &&\
	tar -czvf MWT_Calculator.tgz MWT_Calculator && mv MWT_Calculator.tgz .. &&\
	cd ../ && rm -rf .MWT_Calculator_tmp
	
.PHONY: clean tar
