# pour djgpp
# COMP=gcc -c -O3
# LINK=gxx
# LFLAGS= -lm

# pour gnu-win32
# COMP=c++ -c -O3 -I\cygnus\H-i386-cygwin32\i386-cygwin32\include
# LINK=g++
# LFLAGS= -lm

# pour min-gnu-win32
# COMP=c++ -c -O 
# LINK=c++
# LFLAGS= -lm

HOMEEasyDyn=.
INCLUDEPATH=$(HOMEEasyDyn)/include
HPATH=$(HOMEEasyDyn)/include/EasyDyn

# pour gcc UNIX
#COMP=c++ -c -O -I$(INCLUDEPATH) -DEIGEN
COMP=c++ -c -O -I$(INCLUDEPATH)
LINK=c++ -L$(HOMEEasyDyn)
LFLAGS=-L. -lEasyDyn

EDOBJS = toeol.o vec.o spline.o visu.o sim.o mbs.o

EasyDyn.a : $(EDOBJS)
	ar r libEasyDyn.a $(EDOBJS)

mbs.o : mbs.cpp $(HPATH)/mbs.h $(HPATH)/sim.h $(HPATH)/vec.h
	$(COMP) mbs.cpp

sim.o : sim.cpp $(HPATH)/sim.h
	$(COMP) sim.cpp

visu.o : visu.cpp $(HPATH)/visu.h $(HPATH)/vec.h $(HPATH)/spline.h
	$(COMP) visu.cpp

spline.o : spline.cpp $(HPATH)/spline.h $(HPATH)/toeol.h
	$(COMP) spline.cpp

vec.o : vec.cpp $(HPATH)/vec.h
	$(COMP) vec.cpp

toeol.o : toeol.cpp $(HPATH)/toeol.h
	$(COMP) toeol.cpp
