.c.o:
	cc -o $*.o -c $*.c

test: src/algorithms.o src/test.c
	cc -o test src/algorithms.o src/test.c
