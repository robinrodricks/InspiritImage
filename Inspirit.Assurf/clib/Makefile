CC      = gcc
CC2     = g++
WINFLAGS  = -O3 -Wall
MACFLAGS  = -DOSX -O3 -Wall
LDFLAGS = -O3 -Wall
LDLIBS  = -Llibfast -lfast -Lkdtree -lkdtree -Lrpp -lrpp

PROG  = assurf_clib.swc

# Files extensions .cpp, .c, .o
SUFFIXES = .cpp .c .o 
.SUFFIXES: $(SUFFIXES) .

OBJS  = main.o utils.o rpp.o lkpyramid.o
 
all: $(PROG)

# Compilation and link
$(PROG): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -swc -o $(PROG) $(LDLIBS)
	mv -f $(PROG) ../libs
	-rm -f *.o
	
.c.o:
	$(CC) $(WINFLAGS) -c $< -o $@
	
.cpp.o:
	$(CC2) $(WINFLAGS) -c $< -o $@

clean:
	-rm -f $(PROG)
	-rm -f *.o