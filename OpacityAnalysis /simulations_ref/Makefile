
ifeq ($(USER),linaissa)
	CC = gcc-8
        INCLUDE=/Users/linaissa/gsl/include/
        LIB=/Users/linaissa/gsl/lib/

	CFLAGS = -O2 -Wall -I$(INCLUDE)
	LDFLAGS = -L$(LIB) -lgsl -lgslcblas -lm
endif

ifeq ($(USER),lissa)
	CC = gcc
	INCLUDE=/home/gsl/include/
	LIB=/home/lissa/gsl/lib/

	CFLAGS = -O2 -Wall -I$(INCLUDE)
	LDFLAGS = -L$(LIB) -lgsl -lgslcblas -lm
endif


DEPS = mylib.h

OBJ =  possis.o mylib.o

all: possis

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDFLAGS)

possis: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f *o