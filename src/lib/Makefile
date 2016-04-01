#.........................................................................
# Version "$Id$"
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
# Makefile for SMOKE library
#
#       ---------------     Definitions:   -------------------------
#
.SUFFIXES: .m4 .c .F .f

include ../Makeinclude

SRCDIR = ${BASEDIR}/lib

#       HACK:
#       Partition the source and object lists to get around
#       "Command line too long" errors:

FSRC = \
alocatbl.f alocchrt.f alocctbl.f alocetbl.f alocgtbl.f alocmtbl.f \
alocptbl.f alocstbl.f alocttbl.f applreac.f applumat.f bldcsrc.f \
bldenams.f blkorcmt.f calcrelhum.f checkmem.f chkcpvar.f chkemepi.f \
chkexpscc.f chkexpsic.f chkgrid.f chkint.f chkisiz.f chkmetem.f \
chkptdef.f chkreal.f chksrcno.f convrtll.f convrtxy.f cvtrdtype.f \
cvtvehtype.f dscm3grd.f dscm3lay.f dscsprof.f efsetup.f evalcrit.f \
fillatbl.f fillchrt.f fillctbl.f filletbl.f fillgtbl.f fillmtbl.f \
fillptbl.f fillstbl.f fillttbl.f find1first.f flterr.f fltrneg.f \
fltrxref.f fmtcsrc.f genptcel.f genptvcel.f genuslst.f getbasyr.f \
getcfdsc.f getctgry.f getdysav.f getfline.f getformt.f getidasz.f \
getifdsc.f getiname.f getinvyr.f getm3epi.f getnlist.f getntisz.f \
getrfdsc.f getsinfo.f gettzone.f grd2cnty.f grdfips.f hdrmiss3.f \
ingrid.f initem.f ioapi_grd_size.f lnk2grd.f multunit.f nextime.f \
normtpro.f openphys.f padnzero.f padzero.f parscsrc.f parsline.f \
pdsetup.f polmesg.f prclinrc.f procspro.f progdesc.f rd3mask.f \
rdar2pt.f rdascc.f rdchrscc.f rdcodnam.f rddates.f rdeproc.f \
rdgeocodes.f rdgmat.f rdgref.f rdhdays.f rdinvchr.f rdinvmap.f \
rdinvpol.f rdlines.f rdmactdsc.f rdmapmask.f rdmappol.f rdnaicsdsc.f \
rdorsdsc.f rdpelv.f rdrmat.f rdsccdsc.f rdsccmap.f rdsconv.f \
rdsetmask.f rdsicdsc.f rdsmat.f rdspdprof.f rdspdref.f rdsprof.f \
rdsrcgrps.f rdsref.f rdsrgdesc.f rdsrg.f rdsrghdr.f rdstcy.f \
rdtref.f rdtzone.f rdumat.f rdxclude.f readwr3.f \
rmcommnt.f setscctype.f setsrcdy.f srcgrpcnt.f tagtable.f ungridbv.f \
unitfac.f unitmatch.f useexpgeo.f verchar.f wrchrscc.f wrdaymsg.f \
wrorlout.f wrsrcgrps.f xreftbl.f

OBJ = $(FSRC:.f=.o)

LIB  = libsmoke.a

#      ----------------------   TOP-LEVEL TARGETS:   ------------------

all: lib

clean:
	cd ${OBJDIR}; rm $(LIB)
	cd ${OBJDIR}; rm $(OBJ)

install: $(INSTDIR)
	cd ${OBJDIR}; cp $(LIB) $(INSTDIR)

lib: ${LIB}


#      -----------------------   RULES:   -------------------------

.F.o:
	cd ${OBJDIR}; $(FC) $(FPPFLAGS) $(FFLAGS) -c $(SRCDIR)/$<

.f.o:
	cd ${OBJDIR}; $(FC) $(FFLAGS) -c $(SRCDIR)/$<


#      -----------------------   DEPENDENCIES:------------------------

${LIB}: ${OBJ}
	cd ${OBJDIR}; $(AR) ${ARFLAGS} ${LIB} ${OBJ} ; ranlib ${LIB}

alocatbl.o : $(MODXREF)
alocchrt.o : $(MODXREF)
alocctbl.o : $(MODXREF)
alocetbl.o : $(MODXREF)
alocgtbl.o : $(MODXREF)
alocmtbl.o : $(MODXREF)
alocptbl.o : $(MODXREF)
alocstbl.o : $(MODXREF)
alocttbl.o : $(MODXREF)
applumat.o : $(MODGRID)
bldenams.o : $(MODLISTS)
chkgrid.o : $(MODGRID)
chkptdef.o : $(MODINFO)
dscsprof.o : $(MODINFO)
dscsprof.o : $(MODSPRO)
efsetup.o : $(MODEMFAC)
fillatbl.o : $(MODXREF)
fillchrt.o : $(MODXREF)
fillchrt.o : $(MODINFO)
fillctbl.o : $(MODXREF)
filletbl.o : $(MODXREF)
fillgtbl.o : $(MODXREF)
fillmtbl.o : $(MODXREF)
fillptbl.o : $(MODXREF)
fillstbl.o : $(MODXREF)
fillttbl.o : $(MODXREF)
fltrxref.o : $(MODSOURC)
fltrxref.o : $(MODLISTS)
fltrxref.o : $(MODINFO)
fmtcsrc.o : $(MODINFO)
genptvcel.o : $(MODGRID)
genuslst.o : $(MODSOURC)
genuslst.o : $(MODLISTS)
genuslst.o : $(MODINFO)
getbasyr.o : $(MODSOURC)
getctgry.o : $(MODINFO)
getntisz.o : $(MODLISTS)
getsinfo.o : $(MODFILESET)
getsinfo.o : $(MODINFO)
gettzone.o : $(MODSTCY)
grd2cnty.o : $(MODSTCY)
grd2cnty.o : $(MODSURG)
grd2cnty.o : $(MODGRID)
grdfips.o : $(MODSURG)
grdfips.o : $(MODGRID)
grdfips.o : $(MODMET)
lnk2grd.o : $(MODGRID)
normtpro.o : $(MODTMPRL)
prclinrc.o : $(MODREPRT)
prclinrc.o : $(MODINFO)
procspro.o : $(MODSPRO)
rdar2pt.o : $(MODXREF)
rdar2pt.o : $(MODLISTS)
rdar2pt.o : $(MODAR2PT)
rdcodnam.o : $(MODLISTS)
rddates.o : $(MODTMPRL)
rdeproc.o : $(MODEMFAC)
rdeproc.o : $(MODINFO)
rdgeocodes.o : $(MODSTCY)
rdgref.o : $(MODXREF)
rdgref.o : $(MODINFO)
rdhdays.o : $(MODTMPRL)
rdinvchr.o : $(MODSOURC)
rdinvmap.o : $(MODINFO)
rdinvpol.o : $(MODFILESET)
rdmactdsc.o : $(MODLISTS)
rdmappol.o : $(MODINFO)
rdmappol.o : $(MODFILESET)
rdnaicsdsc.o : $(MODLISTS)
rdorsdsc.o : $(MODSTCY)
rdpelv.o : $(MODELEV)
rdsccdsc.o : $(MODLISTS)
rdsccmap.o : $(MODMOBIL)
rdsconv.o : $(MODLISTS)
rdsconv.o : $(MODSPRO)
rdsconv.o : $(MODINFO)
rdsicdsc.o : $(MODLISTS)
rdsmat.o : $(MODFILESET)
rdspdprof.o : $(MODMOBIL)
rdspdref.o : $(MODXREF)
rdsprof.o : $(MODSPRO)
rdsprof.o : $(MODLISTS)
rdsrcgrps.o : $(MODMERGE)
rdsrcgrps.o : $(MODLISTS)
rdsrcgrps.o : $(MODGRID)
rdsrcgrps.o : $(MODELEV)
rdsrcgrps.o : $(MODSURG)
rdsrcgrps.o : $(MODSOURC)
rdsref.o : $(MODXREF)
rdsref.o : $(MODINFO)
rdsrgdesc.o : $(MODSURG)
rdsrg.o : $(MODSURG)
rdsrg.o : $(MODGRID)
rdstcy.o : $(MODSTCY)
rdtref.o : $(MODXREF)
rdtref.o : $(MODINFO)
rdtref.o : $(MODTMPRL)
rdxclude.o : $(MODXREF)
rdxclude.o : $(MODINFO)
setscctype.o : $(MODINFO)
setsrcdy.o : $(MODINFO)
srcgrpcnt.o : $(MODMERGE)
srcgrpcnt.o : $(MODGRID)
srcgrpcnt.o : $(MODELEV)
tagtable.o : $(MODXREF)
tagtable.o : $(MODSPRO)
tagtable.o : $(MODTAG)
tagtable.o : $(MODINFO)
wrchrscc.o : $(MODLISTS)
wrchrscc.o : $(MODINFO)
wrorlout.o : $(MODSOURC)
wrorlout.o : $(MODSTCY)
wrorlout.o : $(MODLISTS)
wrorlout.o : $(MODINFO)
wrsrcgrps.o : $(MODMERGE)
wrsrcgrps.o : $(MODGRID)
wrsrcgrps.o : $(MODELEV)
xreftbl.o : $(MODXREF)
xreftbl.o : $(MODINFO)
