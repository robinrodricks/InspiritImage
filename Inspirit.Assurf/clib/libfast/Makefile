CC      = gcc
WINFLAGS  = -O3 -Wall
LDLIBS  = 

PROG  = libfast.a

SUFFIXES = .cpp .c .o 
.SUFFIXES: $(SUFFIXES) .

OBJS  = fast.o fast9.o fast_10.o fast_11.o nonmax.o
 
all: $(PROG)

$(PROG): $(OBJS)
	ar r $(PROG) $(OBJS)

.c.o:
	$(CC) $(WINFLAGS) -c $< -o $@

clean:
	-rm -f $(PROG)
	-rm -f *.o
