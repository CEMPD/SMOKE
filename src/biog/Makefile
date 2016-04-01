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

SRCDIR = ${BASEDIR}/biog


FSRC = \
  czangle.f    getpar.f   getparb.f    hrbeis_361.f    hrbeis3.f    prebmet.f \
   rdbpro.f  tmpbeis3.f tmpbeis314.f   tmpbeis361.f       hrno.f  normbeis3.f \
  rdb3fac.f   rdb4fac.f rdb4fac_csv.f normbeis314.f  normbeis361.f

normbeis314.o : $(MODBEIS3)

OBJ  =  $(FSRC:.f=.o)

EXE  = normbeis3  tmpbeis3


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

normbeis314.o : $(MODBEIS3)
normbeis361.o : $(MODBEIS3)
tmpbeis314.o : $(MODSPRO)
tmpbeis314.o : $(MODBEIS3)
tmpbeis361.o : $(MODSPRO)
tmpbeis361.o : $(MODBEIS3)

#  ---------------------------  ${EXE} Program builds:  -----------------

tmpbeis3:  czangle.o getpar.o getparb.o hrbeis_361.o hrbeis3.o prebmet.o rdbpro.o tmpbeis3.o \
           tmpbeis314.o tmpbeis361.o hrno.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

normbeis3: normbeis3.o rdb3fac.o rdb4fac.o rdb4fac_csv.o normbeis314.o normbeis361.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)
