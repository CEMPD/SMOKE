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

SRCDIR = ${BASEDIR}/cntlmat

FSRC =  \
cntlmat.f   alocpkts.f  asgncntl.f  errpkts.f   \
fillcdat.f  fillcntl.f  genmultc.f  genproj.f   genreact.f  \
opencmat.f  openctmp.f  openpmat.f  openrmat.f  pktloop.f   \
procpkts.f  rdpacket.f  wcntlrep.f  wrctmp.f    wrrmat.f

OBJ  =  $(FSRC:.f=.o)

EXE  = cntlmat


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

alocpkts.o : $(MODCNTRL)
asgncntl.o : $(MODSOURC)
asgncntl.o : $(MODXREF)
asgncntl.o : $(MODINFO)
asgncntl.o : $(MODCNTRL)
cntlmat.o : $(MODINFO)
cntlmat.o : $(MODFILESET)
fillcdat.o : $(MODCNTRL)
fillcntl.o : $(MODXREF)
fillcntl.o : $(MODINFO)
genmultc.o : $(MODSOURC)
genmultc.o : $(MODCNTRL)
genmultc.o : $(MODINFO)
genproj.o : $(MODSOURC)
genproj.o : $(MODCNTRL)
genproj.o : $(MODINFO)
genreact.o : $(MODSOURC)
genreact.o : $(MODCNTRL)
genreact.o : $(MODSPRO)
genreact.o : $(MODINFO)
opencmat.o : $(MODCNTRL)
opencmat.o : $(MODINFO)
opencmat.o : $(MODFILESET)
openpmat.o : $(MODINFO)
openpmat.o : $(MODCNTRL)
openpmat.o : $(MODFILESET)
openrmat.o : $(MODINFO)
openrmat.o : $(MODFILESET)
pktloop.o : $(MODXREF)
pktloop.o : $(MODINFO)
procpkts.o : $(MODSOURC)
procpkts.o : $(MODXREF)
procpkts.o : $(MODCNTRL)
procpkts.o : $(MODINFO)
rdpacket.o : $(MODINFO)
wcntlrep.o : $(MODSOURC)
wcntlrep.o : $(MODCNTRL)
wcntlrep.o : $(MODINFO)
wrctmp.o : $(MODINFO)

#  ---------------------------  ${EXE} Program builds:  -----------------

${OBJDIR}:
	mkdir -p ${OBJDIR}

cntlmat: $(OBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

