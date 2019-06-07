#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 3mar2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva  Standardization
# 16Mar2007  Kokron    Disable default O3 optimization; add LOOP_VECT optimization
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../..
endif
this_make := $(lastword $(MAKEFILE_LIST))

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

test:
	@ echo "$(THIS)_$@: ARCH = $(ARCH)"
	@ echo "$(THIS)_$@: MACH = $(MACH)"
	@ echo "$(THIS)_$@: SITE = $(SITE)"
	@ echo "$(THIS)_$@: NODE = $(NODE)"
	@ echo "$(THIS)_$@:"
	@ echo "$(THIS)_$@: ESMADIR = $(ESMADIR)"
	@ echo "$(THIS)_$@: BASEDIR = $(BASEDIR)"
	@ echo "$(THIS)_$@:"
	@ echo "$(THIS)_$@: FC             = $(FC)"
	@ echo "$(THIS)_$@: FCFLAGS_local  = $(FCFLAGS_local)"
	@ echo "$(THIS)_$@: FPPFLAGS_local = $(FPPFLAGS_local)"

%.echo:
	@ echo " $* = $($*)"

THIS = $(shell basename `pwd`)
LIB  = lib$(THIS).a
.PHONY: lib

esma_install install: $(LIB) $(BIN)
esma_install install:
	@ echo ">>> $(THIS)_$@, with \"$(LIB)\" etc. <<<"
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a         $(ESMALIB)
	$(CP) -p *.mod       $(ESMAINC)/$(THIS)

esma_clean clean:
	$(RM) *~ *.[aox] *.[Mm][Oo][Dd]
	$(RM) CRTM_Module.F90

esma_distclean distclean:
	$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]
	$(RM) CRTM_Module.F90

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"


esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"
	@echo "        FREAL = $(FREAL)"

#                  --------------------
#                  User Defined Targets
#                  --------------------

include make.filelist
SRCS := $(SRC_FILES)

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

list:
	@ for d in $(DEPS); do echo $$d; done

#ifeq ("$(FOPT)","-O3")
#   FOPT += $(LOOP_VECT)
#endif
#FREAL = $(FREAL4)
#
#FFLAGS =
#USER_FFLAGS = -IPF_fp_relaxed -ip $(BIG_ENDIAN) $(BYTERECLEN)
#USER_FFLAGS = -IPF_fp_relaxed -ip $(BYTERECLEN)

#MOD_DIRS = . $(INC_MPEU)

###=========================================================================
### Try to override default rules with a local FCFLAGS_local, and local
### "%.o: %.[Ff]90" pattern rules, to ensure the same compilation setting
### as in its previous tag.

FCFLAGS_local = -free -I. -O2 -convert little_endian -assume byterecl

%.F90: %.fpp
	$(CP) $*.fpp $@

%.o: %.F90
	@ echo -n "$(THIS)_$(MAKE): "
	$(FC) $(FCFLAGS_local) $(FPPFLAGS_local) -c $<

%.o: %.f90
	@ echo -n "$(THIS)_$(MAKE): "
	$(FC) $(FCFLAGS_local) -c $<

# For comparison, this is conceptually equivalent to,
#
#	FCFLAGS_local	= $(EXTENDED_SOURCE) $(FINCS) $(FMODS) $(FOPT) $(LITTLE_ENDIAN) $(BYTERECLEN)
#
# while F90FLAGS itself is defined in ESMA_base.mk, more or less
#	F90FLAGS = $(FDEFS) $(FINCS) $(FMODS) $(FOPT) $(FREAL) $(FINT) $(XFLAGS) $(USER_FLAGS)
#

MOD_DIRS = .
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

vpath % $(MOD_DIRS)

%.x : $(LIB) %.o 
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB) $(LIB_SYS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

lib: $(LIB)
$(LIB): $(OBJS)
	@ echo ">>> $(THIS)_$@ <<<"
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

__ignored__ :=  iso_fortran_env.mod
-include ignore.mk
-include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
