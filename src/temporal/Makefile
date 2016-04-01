#.........................................................................
# Version "@(#)$Header$"
# Copyright (C) 2009 Baron Advanced Meteorological Systems, LLC.
# Distributed under the GNU GENERAL PUBLIC LICENSE version 2
# See file "GPL.txt" for conditions of use.
#.........................................................................
#  Environment Variables:
#       BIN     machine/OS/compiler/mode type. Shows up as suffix
#               for "Makeinclude.${BIN}" to determine compilation
#               flags, and in ${OBJDIR} and $(INSTALL) to determine
#               binary directories
#       INSTALL installation-directory root, used for "make install":
#               SMOKE executables will be installed in $(INSTALL)/${BIN}
#.........................................................................
#  Directories:
#       ${BASEDIR}  is the root directory for the SMOKE source and
#                   the (machine/compiler/flag-specific) binary
#                   object/library/executable directories.
#       ${SRCDIR}   is the source directory for this Makefile
#       ${INCDIR}   is the source directory for SMOKE INCLUDE-files
#       ${IODIR}    is the source directory for the I/O API library
#       ${OBJDIR}   is the current machine/compiler/flag-specific
#                   build-directory
#.........................................................................
#
#       ---------------     Definitions:   -------------------------

.SUFFIXES: .m4 .c .F .f .mod

include ../Makeinclude

SRCDIR = ${BASEDIR}/temporal

FSRC = \
asgntpro.f  genhemis.f  mktmat.f    opentmp.f   opentmpin.f \
proctpro.f  setdaylt.f  temporal.f  tmnamunt.f  updtmat.f   \
wrtsup.f

OBJ  =  $(FSRC:.f=.o)

EXE  = temporal


#      ----------------------   TOP-LEVEL TARGETS:   ------------------

all: ${OBJDIR} ${EXE}

clean:
	cd ${OBJDIR}; /bin/rm ${OBJ} ${EXE}

rmexe:
	cd ${OBJDIR}; rm ${EXE}

relink:
	make -i rmexe ; make all

install:  all
	cd ${OBJDIR}; cp ${EXE} ${INSTDIR}


#      -----------------------   RULES:   -----------------------------

%.o : %.mod        #  Disable "gmake"s obnoxious implicit Modula-2 rule !!

.f.o:
	cd ${OBJDIR}; $(FC) $(FFLAGS) -c $(SRCDIR)/$<


#      -----------------------   DEPENDENCIES:------------------------

asgntpro.o : $(MODSOURC)
asgntpro.o : $(MODTMPRL)
asgntpro.o : $(MODINFO)
genhemis.o : $(MODSOURC)
genhemis.o : $(MODXREF)
genhemis.o : $(MODTMPRL)
genhemis.o : $(MODDAYHR)
genhemis.o : $(MODINFO)
mktmat.o : $(MODSOURC)
mktmat.o : $(MODTMPRL)
mktmat.o : $(MODDAYHR)
opentmp.o : $(MODINFO)
opentmp.f: $(MODTMPRL)
opentmp.o : $(MODFILESET)
opentmpin.o : $(MODINFO)
proctpro.o : $(MODINFO)
proctpro.o : $(MODSOURC)
proctpro.o : $(MODTMPRL)
setdaylt.o : $(MODSOURC)
setdaylt.o : $(MODSTCY)
setdaylt.o : $(MODINFO)
temporal.o : $(MODSOURC)
temporal.o : $(MODTMPRL)
temporal.o : $(MODEMFAC)
temporal.o : $(MODDAYHR)
temporal.o : $(MODLISTS)
temporal.o : $(MODINFO)
temporal.o : $(MODMBSET)
tmnamunt.o : $(MODINFO)
updtmat.o : $(MODSOURC)
updtmat.o : $(MODTMPRL)
updtmat.o : $(MODDAYHR)
wrtsup.o : $(MODINFO)
wrtsup.o : $(MODSOURC)
wrtsup.o : $(MODTMPRL)

#  ---------------------------  ${EXE} Program builds:  -----------------

${OBJDIR}:
	mkdir -p ${OBJDIR}

temporal: asgntpro.o genhemis.o mktmat.o opentmp.o opentmpin.o  \
          proctpro.o setdaylt.o temporal.o tmnamunt.o updtmat.o \
          wrtsup.o 
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

