# Makefile to build vanim under X/Windows
# Either the Athena (EasyAnimx) and Motif (EasyAnimm) versions are built
# Tested only under Linux
# You will need to install and build the V library

# Here is our preferred compiler

CC	=	c++

# Adapt HOMEV to your own configuration
HOMEV	=	../../v

# For Motif
LIBSM	=	-lV -lVgl  -lGLU -lGL -lXmu -lXext -lXm -lXt -lX11
# For Athena
LIBSX	=	-lVx -lVxgl -lGLU -lGL -lXmu -lXext -lXaw -lXt -lX11

X11LIB  =       /usr/X11R6/lib
VLibDir	=	$(HOMEV)/lib
oDir	=	.
Bin	=	.

# Flags for includes and libraries

#CFLAGS	=	-O -I$(HOMEV)/include
CFLAGS	=	-O -I$(HOMEV)/includex -I/usr/X11R6/include -fpermissive
LFLAGS	=	-O -L$(VLibDir) -L/usr/X11R6/lib
EXOBJS  =       $(oDir)/EasyAnim.o \
                $(oDir)/EasyAnimcnv.o \
                $(oDir)/EasyAnimcmdw.o \
                $(oDir)/EasyAnimmdlg.o \
                $(oDir)/EasyAnimGL.o

# Let's go to program building

all:    $(Bin)/EasyAnimx $(Bin)/EasyAnimm

objs:	$(EXOBJS)

$(Bin)/EasyAnimm:        $(EXOBJS)
	$(CC) -o $@ $(LFLAGS) $(EXOBJS) $(LIBSM)

$(Bin)/EasyAnimx:        $(EXOBJS)
	$(CC) -o $@ $(LFLAGS) $(EXOBJS) $(LIBSX)

$(oDir)/EasyAnim.o:      EasyAnim.cpp  EasyAnim.h EasyAnimcnv.h EasyAnimcmdw.h EasyAnimmdlg.h
	$(CC) -c $(CFLAGS) -o $@ $<

$(oDir)/EasyAnimcnv.o:   EasyAnimcnv.cpp  EasyAnimcnv.h EasyAnimGL.h 
	$(CC) -c $(CFLAGS) -o $@ $<

$(oDir)/EasyAnimcmdw.o:      EasyAnimcmdw.cpp EasyAnimcmdw.h EasyAnimcnv.h EasyAnimmdlg.h EasyAnimGL.h 
	$(CC) -c $(CFLAGS) -o $@ $<

$(oDir)/EasyAnimGL.o:      EasyAnimGL.cpp EasyAnimGL.h EasyAnimcnv.h
	$(CC) -c $(CFLAGS) -o $@ $<

$(oDir)/EasyAnimmdlg.o:      EasyAnimmdlg.cpp EasyAnimmdlg.h EasyAnimcmdw.h EasyAnimcnv.h
	$(CC) -c $(CFLAGS) -o $@ $<
