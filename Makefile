CFLAGS = -Inobody

.c.o:
	cc $(CFLAGS) -o $*.o -c $*.c

test: nobody/test.o nobody/nobody.o nobody/kepler.o
	cc $(CFLAGS) -o test nobody/nobody.o nobody/kepler.o nobody/test.o
