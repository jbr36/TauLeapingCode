#*********************************************************************************
# *
# * TauLeaping for GLY on Cu 
# *
# * Copyright (c) 2013-2016 J.B. Rommel
# * All rights reserved.
# *
# * Redistribution and use in source and binary forms, with or without
# * modification, are permitted provided that the following conditions are met:
# *
# * 1. Redistributions of source code must retain the above copyright notice, this
# * list of conditions and the following disclaimer.
# * 2. Redistributions in binary form must reproduce the above copyright notice,
# * this list of conditions and the following disclaimer in the documentation
# * and/or other materials provided with the distribution.
# *
# * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# * 
# *********************************************************************************/

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
