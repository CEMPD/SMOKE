# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 

# Makefile for SMOKE library

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<
.F.o:
	$(FC) $(FFLAGS) -c -o $@ $<

include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

FILES_ABCDE = \
  alocatbl.o   alocchrt.o   alocctbl.o   alocetbl.o   alocgtbl.o   alocmtbl.o \
  alocptbl.o   alocstbl.o   alocttbl.o   applreac.o   applumat.o   bldcsrc.o  \
  bldenams.o   blkorcmt.o   chkcpvar.o   chkemepi.o   chkisiz.o    chkmetem.o \
  chkptdef.o   chksrcno.o   convrtxy.o   cvtrdtype.o  cvtvehtype.o dscsprof.o \
  efsetup.o    evalcrit.o

FILES_FGIMNP = \
  fillatbl.o   fillchrt.o   fillctbl.o   filletbl.o   fillgtbl.o   fillmtbl.o \
  fillptbl.o   fillstbl.o   fillttbl.o   find1first.o fltrxref.o   fmtcsrc.o  \
  genptcel.o   genuslst.o   getbasyr.o   getctgry.o   getdysav.o   getformt.o \
  getidasz.o   getiname.o   getinvyr.o   getntisz.o   getrfdsc.o   getsinfo.o \
  gettzone.o   getvmix.o    grd2cnty.o   ingrid.o     initem.o     mbsccadj.o \
  normtpro.o   parscsrc.o   pdsetup.o    prclinrc.o   procspro.o   progdesc.o

FILES_R = \
  rdar2pt.o    rdascc.o     rdchrscc.o   rdcodnam.o   rdeproc.o    rdgmat.o   \
  rdgref.o     rdhdays.o    rdinvchr.o   rdinvmap.o   rdinvpol.o   rdlines.o  \
  rdmapmask.o  rdmappol.o   rdmvinfo.o   rdorsdsc.o   rdpelv.o     rdrmat.o   \
  rdsccdsc.o   rdsconv.o    rdsicdsc.o   rdsmat.o     rdspdprof.o  rdspdref.o \
  rdsprof.o    rdsref.o     rdsrghdr.o   rdsrg.o      rdstcy.o     rdtprof.o  \
  rdtref.o     rdtzone.o    rdumat.o     rdvmix.o     rdxclude.o   readwr3.o  \
  rmcommnt.o

FILES_SVWX = \
  setscctype.o setsrcdy.o   verchar.o    wrchrscc.o   wrdaymsg.o   wridaout.o \
  wridapol.o   xreftbl.o

FILES_TMP = 

# NOTE: At some point, insert dependencies on include files

all: lib

lib: $(MODALL) $(MODGRID) $(FILES_ABCDE) $(FILES_FGIMNP) $(FILES_R) $(FILES_SVWX) $(FILES_TMP) Makefile
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_ABCDE)
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_FGIMNP)
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_R)
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_SVWX)
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_TMP)
	ar rv $(SMK_BIN)/libsmoke.a $(MODALL)

debug: $(MODALL) $(MODGRID) $(FILES_ABCDE) $(FILES_FGIMNP) $(FILES_R) $(FILES_SVWX) $(FILES_TMP) Makefile
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_ABCDE)
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_FGIMNP)
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_R)
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_SVWX)
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_TMP)
	ar rv $(SMK_BIN)/libsmoke.debug.a $(MODALL)

#
# Module dependencies
#

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
chkptdef.o : $(MODINFO)
dscsprof.o : $(MODSPRO)
efsetup.o  : $(MODEMFAC)
fillatbl.o : $(MODXREF)
fillchrt.o : $(MODXREF) $(MODINFO)
fillctbl.o : $(MODXREF) 
filletbl.o : $(MODXREF) 
fillgtbl.o : $(MODXREF) 
fillmtbl.o : $(MODXREF) 
fillptbl.o : $(MODXREF)
fillstbl.o : $(MODXREF)
fillttbl.o : $(MODXREF)
fltrxref.o : $(MODSOURC) $(MODLISTS) $(MODINFO)
fmtcsrc.o  : $(MODINFO)
genuslst.o : $(MODSOURC) $(MODLISTS) $(MODINFO)
getbasyr.o : $(MODSOURC)
getctgry.o : $(MODINFO)
getisize.o : $(MODMOBIL)
getntisz.o : $(MODLISTS)
getsinfo.o : $(MODFILESET) $(MODINFO)
gettzone.o : $(MODSTCY)
getvmix.o  : $(MODXREF) $(MODINFO)
grd2cnty.o : $(MODSTCY) $(MODSURG) $(MODGRID)
mbsccadj.o : $(MODMOBIL)
normtpro.o : $(MODTMPRL)
prclinrc.o : $(MODREPRT) $(MODINFO)
procspro.o : $(MODSPRO)
rdar2pt.o  : $(MODXREF) $(MODLISTS) $(MODAR2PT)
rdcodnam.o : $(MODLISTS)
rdeproc.o  : $(MODEMFAC) $(MODINFO)
rdgref.o   : $(MODXREF) $(MODINFO)
rdhdays.o  : $(MODTMPRL)
rdinvchr.o : $(MODSOURC)
rdinvmap.o : $(MODINFO)
rdinvpol.o : $(MODFILESET) 
rdmappol.o : $(MODINFO) $(MODFILESET)
rdmvinfo.o : $(MODMOBIL)
rdorsdsc.o : $(MODSTCY)
rdpelv.o   : $(MODELEV)
rdsccdsc.o : $(MODLISTS)
rdsicdsc.o : $(MODLISTS)
rdsconv.o  : $(MODLISTS) $(MODSPRO) $(MODINFO)
rdsmat.o   : $(MODFILESET)
rdspdprof.o: $(MODMOBIL)
rdspdref.o : $(MODXREF)
rdsprof.o  : $(MODSPRO)
rdsref.o   : $(MODXREF) $(MODINFO)
rdsrg.o    : $(MODSURG) $(MODGRID)
rdstcy.o   : $(MODSTCY) 
rdtprof.o  : $(MODTMPRL)
rdtref.o   : $(MODXREF) $(MODINFO)
rdvmix.o   : $(MODMOBIL) $(MODXREF)
rdxclude.o : $(MODXREF)
setscctype.o: $(MODINFO)
setsrcdy.o : $(MODINFO)
wrchrscc.o : $(MODLISTS) $(MODINFO)
wridaout.o : $(MODSOURC) $(MODSTCY) $(MODINFO)
xreftbl.o  : $(MODXREF) $(MODINFO)

clean:
	/bin/rm -f *.o
	/bin/rm -f $(MODALL)
