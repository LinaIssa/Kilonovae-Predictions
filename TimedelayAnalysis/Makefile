ifeq ($(USER),linaissa)
	CC = gcc-8
        INCLUDE=/Users/linaissa/gsl/include/
        LIB=/Users/linaissa/gsl/lib/

	CFLAGS = -O2 -Wall -I$(INCLUDE) 
	LDFLAGS = -L$(LIB) -lgsl -lgslcblas -lm
endif

ifeq ($(USER),lissa)
	CC = gcc
	INCLUDE=/gsl/include/
	LIB=/home/lissa/gsl/lib/
#/home/lissa/possis/
	CFLAGS = -O2 -Wall -I$(INCLUDE)
	LDFLAGS = -L$(LIB) -lgsl  -lgslcblas -lm 

endif

ifeq ($(USER),gwmilano)
	CC = gcc
	CFLAGS = -O2 -Wall
	LDFLAGS =  -lgsl -lgslcblas -lm
endif


DEPS = mylib.h

OBJ = kn_specpol.o mylib.o

all: kn_specpol

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDFLAGS)

kn_specpol: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f *o