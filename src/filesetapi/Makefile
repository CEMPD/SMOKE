#
#.........................................................................
# Version "@(#)$Header$"
# WRF-Chem/RTOFS SMOKE
# Copyright (C) 2005 Baron Advanced Meteorological Systems, LLC.
# Distributed under the GNU GENERAL PUBLIC LICENSE version 2
# See file "GPL.txt" for conditions of use.
#.........................................................................
# Makefile for FileSetAPI library
#
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

.SUFFIXES: .f .F .m4 .c .mod .MOD

include ../Makeinclude

SRCDIR = ${BASEDIR}/filesetapi

MSRC = modfileset.f

fSRC = \
  appendname.f  chkfileset.f  chksetdesc.f  cleanup.f     closeset.f  \
  createset.f   descset.f     openset.f     promptset.f

FSRC = readset.F  writeset.F

OBJ =  $(MSRC:.f=.o) $(fSRC:.f=.o) $(FSRC:.F=.o)

MOD = $(MSRC:.f=.mod)

LIB  = libfileset.a


#      ----------------------   TOP-LEVEL TARGETS:   ------------------

all: ${LIB}

clean:
	cd ${OBJDIR} ; rm ${OBJ} ${MOD} ${LIB}

install: $(INSTDIR)
	cd ${OBJDIR}; cp ${MOD} ${LIB} $(INSTDIR)

${LIB}: ${OBJ}
	cd ${OBJDIR}; $(AR) ${ARFLAGS} ${LIB} ${OBJ} ; ranlib ${LIB}


#      -----------------------   RULES:   -------------------------

%.o : %.mod        #  Disable "gmake"s obnoxious implicit Modula-2 rule !!

.F.o:
	cd ${OBJDIR}; $(FC) $(FPPFLAGS) $(FFLAGS) -c $(SRCDIR)/$<

.f.o:
	cd ${OBJDIR}; $(FC) $(FFLAGS) -c $(SRCDIR)/$<

.f.mod:
	cd ${OBJDIR}; $(FC) $(FFLAGS) -c $(SRCDIR)/$<


#      -----------------------   DEPENDENCIES:------------------------

chkfileset.o: $(MODFILESET)
chksetdesc.o: $(MODFILESET)
cleanup.o:    $(MODFILESET)
closeset.o:   $(MODFILESET)
createset.o:  $(MODFILESET)
descset.o:    $(MODFILESET)
openset.o:    $(MODFILESET)
readset.o:    $(MODFILESET)
writeset.o:   $(MODFILESET)
