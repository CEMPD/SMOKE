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

SRCDIR = ${BASEDIR}/smkmerge

FSRC =  \
allocmrg.f   bldmrgidx.f   getmrgev.f  initstcy.f    mrgasciielev.f mrgbio.f   \
mrgelev.f    mrggrid.f     mrgmult.f   mrgonams.f    mrgunits.f     mrgvnams.f \
openmrgin.f  openmrgout.f  rdmrginv.f  setoutdate.f  smkmerge.f     wmrgelev.f \
wmrgemis.f   wrmrgrep.f

OBJ  =  $(FSRC:.f=.o)

EXE  = smkmerge mrggrid mrgelev


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

allocmrg.o : $(MODMERGE)
allocmrg.o : $(MODELEV)
allocmrg.o : $(MODCNTRL)
allocmrg.o : $(MODSTCY)
allocmrg.o : $(MODGRID)
bldmrgidx.o : $(MODMERGE)
getmrgev.o : $(MODMERGE)
getmrgev.o : $(MODINFO)
initstcy.o : $(MODMERGE)
initstcy.o : $(MODSTCY)
mrgelev.o : $(MODMERGE)
mrgelev.o : $(MODCNTRL)
mrgelev.o : $(MODELEV)
mrggrid.o : $(MODGRID)
mrggrid.o : $(MODFILESET)
mrgmult.o : $(MODMERGE)
mrgmult.o : $(MODELEV)
mrgmult.o : $(MODSTCY)
mrgonams.o : $(MODMERGE)
mrgunits.o : $(MODMERGE)
mrgvnams.o : $(MODMERGE)
mrgvnams.o : $(MODLISTS)
openmrgin.o : $(MODMERGE)
openmrgin.o : $(MODINFO)
openmrgin.o : $(MODELEV)
openmrgin.o : $(MODGRID)
openmrgin.o : $(MODFILESET)
openmrgout.o : $(MODMERGE)
openmrgout.o : $(MODELEV)
openmrgout.o : $(MODGRID)
openmrgout.o : $(MODFILESET)
openmrgout.o : $(MODMERGE)
rdmrginv.o : $(MODMERGE)
rdmrginv.o : $(MODLISTS)
rdmrginv.o : $(MODSOURC)
smkmerge.o : $(MODMERGE)
smkmerge.o : $(MODCNTRL)
smkmerge.o : $(MODELEV)
smkmerge.o : $(MODLISTS)
smkmerge.o : $(MODSTCY)
smkmerge.o : $(MODSURG)
smkmerge.o : $(MODGRID)
wmrgelev.o : $(MODMERGE)
wmrgelev.o : $(MODSOURC)
wmrgelev.o : $(MODELEV)
wmrgelev.o : $(MODGRID)
wmrgemis.o : $(MODMERGE)
wmrgemis.o : $(MODELEV)
wrmrgrep.o : $(MODMERGE)
wrmrgrep.o : $(MODSTCY)
wrmrgrep.o : $(MODGRID)

#  ---------------------------  ${EXE} Progrgm builds:  -----------------

${OBJDIR}:
	mkdir -p ${OBJDIR}

smkmerge: smkmerge.o allocmrg.o bldmrgidx.o getmrgev.o mrgmult.o mrgonams.o \
          mrgvnams.o openmrgin.o openmrgout.o wmrgemis.o wrmrgrep.o initstcy.o \
          rdmrginv.o mrgelev.o mrgunits.o mrgbio.o wmrgelev.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

mrggrid: mrggrid.o setoutdate.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

mrgelev: mrgasciielev.o setoutdate.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

