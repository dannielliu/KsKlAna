ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
#INCDIR = -I./include -I$(shell echo ${G4INCLUDE})
vpath %.h ./include
vpath %.o obj
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) #$(INCDIR)

all: Ks_info


%: %.C
	@echo "compiling $@"
	$(CC) -lRooFitCore -lRooFit -lMathMore $^ -o $@

.PHONY:clean
clean:
	-rm -f Ks_info

