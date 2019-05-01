IDIR = ../include
CC=gcc
CFLAGS=-I$(IDIR) -O3

ODIR=obj
LDIR=../lib

LIBS=-lm -lfftw3

_DEPS = delivery.h kitchen.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = manager.o delivery.o kitchen.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

../ipfb: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o ../ipfb
