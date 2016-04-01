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

SRCDIR = ${BASEDIR}/emqa

FSRC =  \
asgnbins.f   bldrepidx.f  genrprt.f   openrepin.f  openrepout.f  qarepin.f      \
rdgrps.f     rdrepin.f    rdrprts.f   repmrggrd.f  repunits.f    scanrepc.f     \
selectsrc.f  smkreport.f wrrephdr.f   wrrepout.f

OBJ  =  $(FSRC:.f=.o)

EXE  = smkreport


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

asgnbins.o : $(MODSOURC)
asgnbins.o : $(MODLISTS)
asgnbins.o : $(MODREPRT)
asgnbins.o : $(MODREPBN)
asgnbins.o : $(MODGRID)
asgnbins.o : $(MODELEV)
asgnbins.o : $(MODSTCY)
asgnbins.o : $(MODINFO)
bldrepidx.o : $(MODREPRT)
bldrepidx.o : $(MODREPBN)
bldrepidx.o : $(MODTMPRL)
bldrepidx.o : $(MODCNTRL)
bldrepidx.o : $(MODINFO)
bldrepidx.o : $(MODFILESET)
genrprt.o : $(MODSOURC)
genrprt.o : $(MODREPRT)
genrprt.o : $(MODREPBN)
genrprt.o : $(MODTMPRL)
genrprt.o : $(MODCNTRL)
genrprt.o : $(MODINFO)
openrepin.o : $(MODREPRT)
openrepin.o : $(MODTMPRL)
openrepin.o : $(MODREPBN)
openrepin.o : $(MODCNTRL)
openrepin.o : $(MODINFO)
openrepin.o : $(MODFILESET)
openrepout.o : $(MODREPRT)
qarepin.o : $(MODREPRT)
qarepin.o : $(MODSTCY)
qarepin.o : $(MODINFO)
rdgrps.o : $(MODREPRT)
rdgrps.o : $(MODGRID)
rdgrps.o : $(MODSTCY)
rdrepin.o : $(MODSOURC)
rdrepin.o : $(MODREPRT)
rdrepin.o : $(MODREPBN)
rdrepin.o : $(MODCNTRL)
rdrepin.o : $(MODELEV)
rdrepin.o : $(MODLISTS)
rdrepin.o : $(MODGRID)
rdrepin.o : $(MODINFO)
rdrepin.o : $(MODFILESET)
rdrprts.o : $(MODREPRT)
rdrprts.o : $(MODINFO)
repmrggrd.o : $(MODSOURC)
repmrggrd.o : $(MODREPRT)
repmrggrd.o : $(MODREPBN)
repmrggrd.o : $(MODGRID)
repunits.o : $(MODREPRT)
repunits.o : $(MODREPBN)
repunits.o : $(MODTMPRL)
repunits.o : $(MODINFO)
scanrepc.o : $(MODREPRT)
scanrepc.o : $(MODINFO)
selectsrc.o : $(MODSOURC)
selectsrc.o : $(MODREPRT)
selectsrc.o : $(MODREPBN)
selectsrc.o : $(MODELEV)
selectsrc.o : $(MODINFO)
smkreport.o : $(MODREPRT)
smkreport.o : $(MODREPBN)
smkreport.o : $(MODGRID)
smkreport.o : $(MODINFO)
wrrephdr.o : $(MODSOURC)
wrrephdr.o : $(MODLISTS)
wrrephdr.o : $(MODREPRT)
wrrephdr.o : $(MODREPBN)
wrrephdr.o : $(MODSTCY)
wrrephdr.o : $(MODGRID)
wrrephdr.o : $(MODINFO)
wrrepout.o : $(MODSOURC)
wrrepout.o : $(MODLISTS)
wrrepout.o : $(MODREPRT)
wrrepout.o : $(MODREPBN)
wrrepout.o : $(MODSTCY)
wrrepout.o : $(MODINFO)

#  ---------------------------  ${EXE} Program builds:  -----------------

${OBJDIR}:
	mkdir -p ${OBJDIR}

smkreport: $(OBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

