CFLAGS=-g -Wall -O2
LDFLAGS = -shared -fPIC
LIBS=-lz -lm -lpthread
PROG=kaci

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

all: $(PROG) lib$(PROG)

kaci: kaci.c kaci.h khashl.h kseq.h logaddexp.h
	$(CC) $(CFLAGS) -o $@ kaci.c kaci.h khashl.h kseq.h logaddexp.h $(LIBS)

libkaci: kaci.c kaci.h khashl.h kseq.h logaddexp.h
	$(CC) $(CFLAGS) -o $@.so kaci.c kaci.h khashl.h kseq.h logaddexp.h $(LIBS) $(LDFLAGS)

clean:
	rm -fr lib$(PROG).so $(PROG)
