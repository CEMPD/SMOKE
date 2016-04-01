#.........................................................................
# Version "@(#)$Header$"
# WRF-Chem/RTOFS SMOKE
# Copyright (C) 2005 Baron Advanced Meteorological Systems, LLC.
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

SRCDIR = ${BASEDIR}/emutil

FSRC =  \
aggwndw.f     beld3to2.f    bluesky2inv.f   cemscan.f    \
extractida.f  gentpro.f     geofac.f        invsplit.f   \
layalloc.f    metcombine.f  metscan.f       pktreduc.f   \
smk2emis.f    surgtool.f    uam2ncf.f       gcntl4carb.f

OBJ  =  $(FSRC:.f=.o)

EXE  = aggwndw  beld3to2  bluesky2inv  cemscan  extractida \
       geofac   invsplit  metcombine   metscan  pktreduc   \
       smk2emis surgtool  uam2ncf      gentpro  gcntl4carb \
       layalloc

#      ----------------------   TOP-LEVEL TARGETS:   ------------------

all: ${EXE}

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

aggwndw.o : $(MODGRID)
bluesky2inv.o : $(MODSTCY)
gentpro.o : $(MODINFO)
gentpro.o : $(MODSOURC)
gentpro.o : $(MODGRID)
gentpro.o : $(MODSURG)
gentpro.o : $(MODSTCY)
geofac.o : $(MODFILESET)
invsplit.o : $(MODINFO)
smk2emis.o : $(MODINFO)
surgtool.o : $(MODSURG)
surgtool.o : $(MODGRID)

#  ---------------------------  ${EXE} Program builds:  -----------------

aggwndw: aggwndw.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

beld3to2: beld3to2.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

bluesky2inv: bluesky2inv.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

cemscan: cemscan.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

extractida: extractida.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

gcntl4carb: gcntl4carb.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

gentpro: gentpro.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

geofac: geofac.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

invsplit: invsplit.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

layalloc: layalloc.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

metcombine: metcombine.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

metscan: metscan.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

pktreduc: pktreduc.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

smk2emis: smk2emis.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

surgtool: surgtool.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

uam2ncf: uam2ncf.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)
