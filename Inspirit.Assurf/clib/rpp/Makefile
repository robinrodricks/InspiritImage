CC      = gcc
LCC     = ar
LDFLAGS = -O3 -Wall

PROG  = librpp.a

# Files extensions .cpp, .c, .o
SUFFIXES = .cpp .c .o 
.SUFFIXES: $(SUFFIXES) .

# Object files .o necessary to build the main program
OBJS  = librpp.o rpp.o rpp_quintic.o rpp_svd.o rpp_vecmat.o
 
all: $(PROG)

# Compilation and link
$(PROG): $(OBJS)
	$(LCC) rcs $(PROG) $(OBJS)

.c.o:
	$(CC) $(LDFLAGS) -c $< -o $@

clean:
	-rm -f $(PROG)
	-rm -f *.l.bc
	-rm -f *.o
