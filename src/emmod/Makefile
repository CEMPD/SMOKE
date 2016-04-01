#
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

.SUFFIXES: .f .F .m4 .c .mod .MOD

include ../Makeinclude

SRCDIR = ${BASEDIR}/emmod

FSRC =  \
modar2pt.f  modbeis3.f  modbiog.f   modcntrl.f \
moddayhr.f  modelev.f   modemfac.f  modgrid.f  \
modinfo.f   modlists.f  modmbset.f  modmerge.f \
modmet.f    modmobil.f  modmvsmrg.f modrepbn.f \
modreprt.f  modsourc.f  modspro.f   modstcy.f  \
modsurg.f   modtag.f    modtmprl.f  modxref.f

OBJ = $(FSRC:.f=.o)

MOD = $(FSRC:.f=.mod)

LIB  = libemmod.a


#      ----------------------   TOP-LEVEL TARGETS:   ------------------

all: ${LIB}

${LIB}:  ${OBJ}
	cd ${OBJDIR}; $(AR) ${ARFLAGS} ${LIB} ${OBJ} ; ranlib ${LIB}

${OBJDIR}:
	cd ${BASEDIR} ; mkdir -p ${BIN}

clean:
	cd ${OBJDIR} ; rm ${OBJ} ${MOD} ${LIB}

install: $(INSTDIR)
	cd ${OBJDIR}; cp ${MOD} ${LIB} $(INSTDIR)


#      -----------------------   RULES:   -----------------------------

%.o : %.mod        #  Disable "gmake"s obnoxious implicit Modula-2 rule !!

.f.o .f.mod:
	cd ${OBJDIR}; $(FC) $(FFLAGS) -c $(SRCDIR)/$<

