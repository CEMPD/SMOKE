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

SRCDIR = ${BASEDIR}/spcmat

FSRC = spcmat.f  asgnspro.f  asgntag.f  chknonhap.f  opensmat.f  rdcombo.f rdtag.f

OBJ  =  $(FSRC:.f=.o)

EXE  = spcmat


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

s01/depts/ie/cempd/bbaek/SMOKEv40_git/subsys/smoke/src/temporal> grep ' USE MOD' *.f|asgnspro.o : $(MODSOURC)
asgnspro.o : $(MODLISTS)
asgnspro.o : $(MODXREF)
asgnspro.o : $(MODSPRO)
asgnspro.o : $(MODINFO)
asgntag.o : $(MODSOURC)
asgntag.o : $(MODTAG)
asgntag.o : $(MODINFO)
chknonhap.o : $(MODLISTS)
chknonhap.o : $(MODSPRO)
opensmat.o : $(MODINFO)
opensmat.o : $(MODTAG)
opensmat.o : $(MODFILESET)
rdcombo.o : $(MODSPRO)
rdcombo.o : $(MODLISTS)
rdtag.o : $(MODXREF)
rdtag.o : $(MODSPRO)
rdtag.o : $(MODINFO)
rdtag.o : $(MODLISTS)
rdtag.o : $(MODTAG)
spcmat.o : $(MODSPRO)
spcmat.o : $(MODEMFAC)
spcmat.o : $(MODINFO)
spcmat.o : $(MODLISTS)
spcmat.o : $(MODTAG)
spcmat.o : $(MODFILESET)

#  ---------------------------  ${EXE} Program builds:  -----------------

${OBJDIR}:
	mkdir -p ${OBJDIR}

spcmat: $(OBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)
