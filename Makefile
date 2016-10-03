# J.B. Rommel, October 2016

SHELL   = sh

TARGET  = TauLeapwithGly.x

MODS	= TauLeapwithGly 
SRCS	= $(MODS:=.f90)
OBJS	= $(MODS:=.o)

#F90     = nagfor
#OPTS    = -C -mtrace -gline -g -O0 -ieee=nonstd -nan
F90     = gfortran
OPTS    = -g -O0 -Wall -fbacktrace -Wextra -ffpe-trap=invalid,zero,overflow,underflow -C 
#F90     = ifort
#OPTS	= -O0 -g -check all -fpe0 -warn -traceback -debug extended
 
#For profiling
#OPTS = -profile-functions -profile-loops=all -profile-loops-report=2
#LINK	= -lmkl_lapack -lmkl -lguide -lpthread -lmi -L../lib/ -lgnufor

#LINK   = -L../../lib/ -lgnufor -lasa047

all: $(TARGET)

$(TARGET): $(OBJS)
	$(F90) -o $@ $(OBJS) $(LINK)

%.o: %.f90
	$(F90) -o $@ $(OPTS) -c $<

.PHONY: clean
clean:
	@-$(RM) -f $(OBJS) $(TARGET) fort.* *~
