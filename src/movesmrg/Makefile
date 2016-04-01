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

SRCDIR = ${BASEDIR}/movesmrg

FSRC =  \
mbldmrgidx.f    bldprocidx.f    bldsrccell.f    mgetmrgev.f     \
hourmet.f       minitstcy.f     met4moves.f     movesmrg.f      \
mmrgonams.f     mmrgunits.f     mmrgvnams.f     mopenmrgin.f    \
mopenmrgout.f   rdcfpro.f       rdfmref.f       \
rdmrclist.f     rdmxref.f       rdrpdemfacs.f   rdrppemfacs.f   \
rdrpvemfacs.f   rdrphemfacs.f   rdspdpro.f      setrefcnty.f    \
mwrmrgrep.f     wrtemprof.f

OBJ  =  $(FSRC:.f=.o)

EXE  =  movesmrg    met4moves


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

mbldmrgidx.o : $(MODMERGE)
mbldmrgidx.o : $(MODMVSMRG)
bldprocidx.o : $(MODMVSMRG)
bldprocidx.o : $(MODMERGE)
bldsrccell.o : $(MODMVSMRG)
mgetmrgev.o : $(MODINFO)
mgetmrgev.o : $(MODMERGE)
mgetmrgev.o : $(MODMVSMRG)
hourmet.o : $(MODMBSET)
hourmet.o : $(MODMET)
hourmet.o : $(MODSURG)
minitstcy.o : $(MODMERGE)
minitstcy.o : $(MODSTCY)
minitstcy.o : $(MODSOURC)
minitstcy.o : $(MODMVSMRG)
minitstcy.o : $(MODLISTS)
met4moves.o : $(MODMBSET)
met4moves.o : $(MODGRID)
met4moves.o : $(MODSURG)
met4moves.o : $(MODSTCY)
met4moves.o : $(MODMET)
movesmrg.o : $(MODMERGE)
movesmrg.o : $(MODMVSMRG)
movesmrg.o : $(MODLISTS)
movesmrg.o : $(MODSTCY)
movesmrg.o : $(MODGRID)
movesmrg.o : $(MODMBSET)
movesmrg.o : $(MODSOURC)
mmrgonams.o : $(MODMERGE)
mmrgunits.o : $(MODMERGE)
mmrgunits.o : $(MODMVSMRG)
mmrgvnams.o : $(MODMERGE)
mmrgvnams.o : $(MODLISTS)
mmrgvnams.o : $(MODMVSMRG)
mopenmrgin.o : $(MODMERGE)
mopenmrgin.o : $(MODMVSMRG)
mopenmrgin.o : $(MODINFO)
mopenmrgin.o : $(MODSOURC)
mopenmrgin.o : $(MODEMFAC)
mopenmrgin.o : $(MODGRID)
mopenmrgin.o : $(MODFILESET)
mopenmrgout.o : $(MODMERGE)
mopenmrgout.o : $(MODGRID)
mopenmrgout.o : $(MODMVSMRG)
mopenmrgout.o : $(MODFILESET)
rdcfpro.o : $(MODMERGE)
rdcfpro.o : $(MODMVSMRG)
rdcfpro.o : $(MODLISTS)
rdcfpro.o : $(MODMBSET)
rdfmref.o : $(MODMBSET)
rdmrclist.o : $(MODMBSET)
rdmrclist.o : $(MODMVSMRG)
rdmxref.o : $(MODMBSET)
rdrpdemfacs.o : $(MODMBSET)
rdrpdemfacs.o : $(MODMVSMRG)
rdrpdemfacs.o : $(MODMERGE)
rdrpdemfacs.o : $(MODLISTS)
rdrphemfacs.o : $(MODMBSET)
rdrphemfacs.o : $(MODMVSMRG)
rdrphemfacs.o : $(MODMERGE)
rdrphemfacs.o : $(MODLISTS)
rdrppemfacs.o : $(MODMBSET)
rdrppemfacs.o : $(MODMVSMRG)
rdrppemfacs.o : $(MODMERGE)
rdrppemfacs.o : $(MODLISTS)
rdrpvemfacs.o : $(MODMBSET)
rdrpvemfacs.o : $(MODMVSMRG)
rdrpvemfacs.o : $(MODMERGE)
rdrpvemfacs.o : $(MODLISTS)
rdspdpro.o : $(MODMVSMRG)
rdspdpro.o : $(MODLISTS)
rdspdpro.o : $(MODMOBIL)
setrefcnty.o : $(MODMVSMRG)
setrefcnty.o : $(MODINFO)
setrefcnty.o : $(MODSOURC)
setrefcnty.o : $(MODLISTS)
setrefcnty.o : $(MODMBSET)
mwrmrgrep.o : $(MODMERGE)
mwrmrgrep.o : $(MODSTCY)
mwrmrgrep.o : $(MODMVSMRG)
mwrmrgrep.o : $(MODLISTS)
wrmrgrep.o : $(MODGRID)
wrtemprof.o : $(MODMET)
wrtemprof.o : $(MODMBSET)

#  ---------------------------  $(EXE) Program builds:  -----------------

movesmrg: movesmrg.o mbldmrgidx.o bldprocidx.o bldsrccell.o \
          mgetmrgev.o minitstcy.o mmrgonams.o mmrgunits.o  mmrgvnams.o   \
          mopenmrgin.o mopenmrgout.o rdfmref.o  rdmrclist.o rdmxref.o \
          rdrpdemfacs.o rdrppemfacs.o rdrpvemfacs.o rdrphemfacs.o rdspdpro.o     \
          rdcfpro.o setrefcnty.o mwrmrgrep.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) $^ ${LIBS} -o $@

met4moves: met4moves.o rdmxref.o rdfmref.o hourmet.o wrtemprof.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) $^ ${LIBS} -o $@


