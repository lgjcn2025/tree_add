CCOMP = gcc
#CFLAGS = -O4 -DNDEBUG -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall -lm
CFLAGS = -g -DPRINT_FLOW -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall -lm

all: ladd 
ladd: main.c ladd.c parser_ladd.c timer.c
	$(CCOMP) $(CFLAGS) -o main_ladd main.c libm.so 
clean: 
	rm -f main_ladd *~
