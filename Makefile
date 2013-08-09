SRC        = src
TEST_FILES = $(SRC)/test.o $(SRC)/nobody.o $(SRC)/kepler.o \
			 $(SRC)/variational.o
CFLAGS     = -Iinclude

.c.o:
	cc $(CFLAGS) -o $*.o -c $*.c

test: $(TEST_FILES)
	cc $(CFLAGS) -o test $(TEST_FILES)
	./test

clean:
	rm -rf $(TEST_FILES) test
