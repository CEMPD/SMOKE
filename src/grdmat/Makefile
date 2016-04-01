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

SRCDIR = ${BASEDIR}/grdmat

FSRC =  \
asgnsurg.f genagmat.f genggmat.f  genlgmat.f genmgmat.f genpgmat.f \
grdmat.f   opengmat.f rdsrg4grd.f sizgmat.f  setfrac.f

OBJ  =  $(FSRC:.f=.o)

MODO = modsourc.o modsurg.o modlists.o modinfo.o modgrid.o modmobil.o modtag.o modspro.o modxref.o modfileset.o

EXE  = grdmat


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

asgnsurg.o : $(MODSOURC)
asgnsurg.o : $(MODXREF)
asgnsurg.o : $(MODSURG)
asgnsurg.o : $(MODINFO)
genagmat.o : $(MODSOURC)
genagmat.o : $(MODSURG)
genagmat.o : $(MODGRID)
genagmat.o : $(MODLISTS)
genagmat.o : $(MODINFO)
genagmat.o : $(MODXREF)
genggmat.o : $(MODSOURC)
genggmat.o : $(MODGRID)
genggmat.o : $(MODINFO)
genlgmat.o : $(MODSOURC)
genlgmat.o : $(MODINFO)
genlgmat.o : $(MODGRID)
genlgmat.o : $(MODGCTP)
genlgmat.o : $(M3UTILIO)
genmgmat.o : $(MODSOURC)
genmgmat.o : $(MODGRID)
genmgmat.o : $(MODSURG)
genmgmat.o : $(MODLISTS)
genmgmat.o : $(MODINFO)
genmgmat.o : $(MODXREF)
grdmat.o : $(MODSOURC)
grdmat.o : $(MODXREF)
grdmat.o : $(MODGRID)
grdmat.o : $(MODINFO)
grdmat.o : $(MODFILESET)
grdmat.o : $(MODSURG)
opengmat.o : $(MODINFO)
rdsrg4grd.o : $(MODSURG)
rdsrg4grd.o : $(MODGRID)
setfrac.o : $(MODSURG)
setfrac.o : $(MODXREF)
sizgmat.o : $(MODSOURC)
sizgmat.o : $(MODXREF)
sizgmat.o : $(MODGRID)
sizgmat.o : $(MODSURG)

#  ---------------------------  ${EXE} Program builds:  -----------------

grdmat: $(OBJ) $(MODO)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

