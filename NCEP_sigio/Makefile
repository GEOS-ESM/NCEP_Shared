#!/bin/make
##########################################################################
#!BOP
#
# !ROUTINE:        Makefile for 'NCEP's sigio module' 
# 
#
# !DESCRIPTION: 
#       create 'NCEP's sigio module' 
#
# !REVISION HISTORY: 
#
# 97.12.16  Meta Sienkiewicz   Initial Makefile for bufron29
# 26Apr2001 da Silva           Made fvDAS compatible
# 
#!EOP
##########################################################################

COREROOT = ../../..
COREBIN  = $(COREROOT)/bin
CORELIB  = $(COREROOT)/lib
COREINC  = $(COREROOT)/include
COREETC  = $(COREROOT)/etc

SHELL		= sh
MAKE		= make
MAKEFILE	= Makefile
MAKECONF	= Makefile.conf

CPP		= /lib/cpp
CPPFLAGS	= -P
FPP	        = $(CPP)

FC              = f90
FFLAGS          = -O3 -extend_source -64
F90		= $(FC)
F90FLAGS	= $(FFLAGS)
F90PP		= $(FPP)
F90PPFLAGS	= $(FPPFLAGS)

LN		= ln -s
RM		= rm -f
CP		= cp -p
MKDIR  = /bin/mkdir -p

ZIP		= compress
ZIPTAG		= .Z
ZCAT		= zcat
TAR		= tar
AR		= ar cq

#..................................................................
.SUFFIXES:
.SUFFIXES: .H .h .F90 .f90 .F .f .o

.f.o:
	$(FC) -c $(FFLAGS) $*.f

.F.o:
	@$(RM) $*.f
	$(FPP) $(FPPFLAGS) $(DEFS) $*.F $*.f
	$(FC) -c $(FFLAGS) $*.f
	@$(RM) $*.f

.f90.o:
	@if [ -f $*.[Ff] ]; then \
	   echo "make: WARNING! coexistent $*.[Ff] and $*.f90" ; fi
	$(F90) -c $(F90FLAGS) $*.f90

.F90.o:
	@ if [ -f $*.[Ff] ]; then \
	   echo "make: WARNING! coexistent $*.[Ff] and $*.F90" ; fi
	@ $(RM) $*.f90
	$(F90PP) $(F90PPFLAGS) $(DEFS) $*.F90 $*.f90
	$(F90) -c $(F90FLAGS) $*.f90
	@ $(RM) $*.f90

# Architecture/site specific costumization
# ----------------------------------------
include $(MAKECONF)

#...................................................................

LIB =	libsigio.a

SRCS  = $(OBSOBJ:.o=.f) 

OBJS  = sigio_module.o 

#...................................................................

lib: $(LIB)

$(LIB): $(OBJS)
	$(RM) $(LIB)
	$(AR) $@ $(OBJS)

export: $(LIB)
	$(MKDIR) $(CORELIB) $(COREINC)/sigio
	$(CP) $(LIB)                  $(CORELIB)
	$(CP) *.[Mm][Oo][Dd]          $(COREINC)/sigio

clean:
	$(RM) -r *~ *.[axol] *.mod *.[fF]*~ Makefile*~ rii_files/ \
	         *.aux *.dvi *.ps  *.toc *.log 

distclean:
	make clean
	$(RM) Makefile.conf


