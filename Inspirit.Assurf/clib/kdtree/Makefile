CC      = gcc
WINFLAGS  = -O3 -Wall
LDLIBS  = 

PROG  = libkdtree.a

SUFFIXES = .cpp .c .o 
.SUFFIXES: $(SUFFIXES) .

OBJS  = kdtree.o random.o
 
all: $(PROG)

$(PROG): $(OBJS)
	ar r $(PROG) $(OBJS)

.c.o:
	$(CC) $(WINFLAGS) -c $< -o $@

clean:
	-rm -f $(PROG)
	-rm -f *.o
