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

SRCDIR = ${BASEDIR}/point

FSRC =  \
asgngrps.f  elevpoint.f   fire_plmris.f  fire_postplm.f  fire_preplm.f  \
laypoint.f  mxgrpemis.f   openeout.f     openlayout.f    plmris.f       \
plsprd.f    plumesti.f    postplm.f      preplm.f        rpelvcfg.f     \
wpingstk.f

OBJ  =  $(FSRC:.f=.o)

EXE  = laypoint elevpoint

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

asgngrps.o : $(MODSOURC)
asgngrps.o : $(MODELEV)
asgngrps.o : $(MODINFO)
elevpoint.o : $(MODSOURC)
elevpoint.o : $(MODELEV)
elevpoint.o : $(MODINFO)
elevpoint.o : $(MODGRID)
laypoint.o : $(MODSOURC)
laypoint.o : $(MODELEV)
laypoint.o : $(MODINFO)
laypoint.o : $(MODGRID)
laypoint.o : $(MODFILESET)
mxgrpemis.o : $(MODSOURC)
mxgrpemis.o : $(MODELEV)
mxgrpemis.o : $(MODINFO)
mxgrpemis.o : $(MODGRID)
mxgrpemis.o : $(MODFILESET)
openeout.o : $(MODINFO)
openeout.o : $(MODELEV)
openlayout.o : $(MODELEV)
openlayout.o : $(MODGRID)
openlayout.o : $(MODINFO)
rpelvcfg.o : $(MODREPRT)
rpelvcfg.o : $(MODELEV)
rpelvcfg.o : $(MODINFO)
wpingstk.o : $(MODELEV)

#  ---------------------------  ${EXE} Program builds:  -----------------

${OBJDIR}:
	mkdir -p ${OBJDIR}

laypoint: laypoint.o openlayout.o preplm.o plmris.o postplm.o plsprd.o \
          fire_preplm.o fire_plmris.o fire_postplm.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

elevpoint: elevpoint.o asgngrps.o openeout.o plumesti.o wpingstk.o  \
           mxgrpemis.o rpelvcfg.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)
