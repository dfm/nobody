.c.o:
	cc -o $*.o -c $*.c

test: test.o nobody.o kepler.o variational.o
	cc $(CFLAGS) -o test nobody/nobody.o nobody/kepler.o nobody/variational.o\
		nobody/test.o

clean:
	rm -rf *.o test
