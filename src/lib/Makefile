# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 

# Makefile for SMOKE library

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc

FILES_A = \
  applumat.o   evalcrit.o   fltrxref.o   getrfdsc.o   procspro.o   rdsrghdr.o \
  alocchrt.o   chkisiz.o    fmtcsrc.o    getsinfo.o   progdesc.o   rdtprof.o  \
  alocctbl.o   chkmetem.o   gettzone.o   rdtref.o     setsrcdy.o   blkorcmt.o \
  alocetbl.o   chkptdef.o   grd2cnty.o   rdascc.o     rdtzone.o  \
  alocgtbl.o   genuslst.o   initem.o     rdchrscc.o   readwr3.o  \
  alocmtbl.o   filletbl.o   getbasyr.o   rdsrg.o      rdsrghdr.o

FILES_B = \
  alocstbl.o   chksrcno.o   rdstcy.o     rdgmat.o     rdinvchr.o \
  alocttbl.o   convrtxy.o   getctgry.o   mbsccadj.o   rdmvinfo.o     \
  applreac.o   rdeproc.o    rdinvpol.o   rdhdays.o    rdsccdsc.o \
  asgnpsi.o    dscsprof.o   getformt.o   rdlines.o    verchar.o  \
  efsetup.o    getidasz.o   normtpro.o   rdpelv.o     tmprinfo.o \
  getdysav.o   ingrid.o     rdrmat.o    

FILES_C = \
  bldcsrc.o    fillchrt.o   getiname.o   parscsrc.o   rdsconv.o    wrchrscc.o \
  bldenams.o   fillctbl.o   getinvyr.o   rdcodnam.o   wrdaymsg.o   wridaout.o \
  fillgtbl.o   getisize.o   pdsetup.o    rdsmat.o     rdumat.o     wridapol.o \
  fillmtbl.o   rdsprof.o    prclinrc.o   rmcommnt.o   rdorsdsc.o \
  chkcpvar.o   fillstbl.o   getvmix.o    rdgref.o     rdvmix.o     xreftbl.o  \
  fillttbl.o   rdefxref.o   rdsref.o     genptcel.o   chkemepi.o  rdiusrc.o

FILES_TMP =

INCLUDES = $(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

# NOTE: A some point, insert dependencies on include files

all: lib

lib: $(MODALL) $(MODGRID) $(FILES_A) $(FILES_B) $(FILES_C) $(FILES_TMP) Makefile
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_A) 
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_B) 
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_C)
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_TMP)
	ar rv $(SMK_BIN)/libsmoke.a $(MODALL)

debug: $(MODALL) $(MODGRID) $(FILES_A) $(FILES_B) $(FILES_C) $(FILES_TMP) Makefile
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_A) 
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_B) 
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_C) 
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_TMP)
	ar rv $(SMK_BIN)/libsmoke.debug.a $(MODALL)

#
# Module dependencies
#

alocchrt.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ alocchrt.f
alocctbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ alocctbl.f
alocetbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ alocetbl.f
alocgtbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ alocgtbl.f
alocmtbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ alocmtbl.f
alocstbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ alocstbl.f
alocttbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ alocttbl.f
applumat.o: $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ applumat.f
asgnpsi.o  : $(MODSOURC) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ asgnpsi.f
bldenams.o  : $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ bldenams.f
chkptdef.o: $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ chkptdef.f
dscsprof.o: $(MODSPRO) 
	$(FC) $(FFLAGS) -c -o $@ dscsprof.f
efsetup.o : $(MODEMFAC) $(INC)/M5CNST3.EXT
	$(FC) $(FFLAGS) -c -o $@ efsetup.f
fillchrt.o: $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ fillchrt.f
fillctbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ fillctbl.f
filletbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ filletbl.f
fillgtbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ fillgtbl.f
fillmtbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ fillmtbl.f
fillstbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ fillstbl.f
fillttbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ fillttbl.f
fltrxref.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODUNCERT)
	$(FC) $(FFLAGS) -c -o $@ fltrxref.f
fmtcsrc.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ fmtcsrc.f
genuslst.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ genuslst.f
getbasyr.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ getbasyr.f
getctgry.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ getctgry.f
getformt.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ getformt.f
getsinfo.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ getsinfo.f
gettzone.o: $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ gettzone.f
getvmix.o: $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ getvmix.f
grd2cnty.o: $(MODSTCY) $(MODSURG) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ grd2cnty.f
lambert.o: 
	$(FC) -static $(FFLAGS) -c -o $@ lambert.f
mbsccadj.o  : $(MODMOBIL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ mbsccadj.f
normtpro.o: $(MODTMPRL)
	$(FC) $(FFLAGS) -c -o $@ normtpro.f
prclinrc.o: $(MODREPRT) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ prclinrc.f
procspro.o: $(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ procspro.f
rdcodnam.o: $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ rdcodnam.f
rdefxref.o: $(MODMOBIL) $(MODXREF) $(MODEMFAC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdefxref.f
rdeproc.o: $(MODEMFAC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdeproc.f
rdgref.o : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdgref.f
rdinvchr.o: $(MODSOURC) 
	$(FC) $(FFLAGS) -c -o $@ rdinvchr.f
rdiusrc.o: $(MODUNCERT) 
	$(FC) $(FFLAGS) -c -o $@ rdiusrc.f
rdmvinfo.o: $(MODMOBIL)
	$(FC) $(FFLAGS) -c -o $@ rdmvinfo.f
rdorsdsc.o : $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ rdorsdsc.f
rdpelv.o : $(MODINFO) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ rdpelv.f
rdsccdsc.o : $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ rdsccdsc.f
rdsconv.o : $(MODLISTS) $(MODSPRO) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdsconv.f
rdsref.o  : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdsref.f
rdsrg.o   : $(MODSURG) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ rdsrg.f
rdsrghdr.o: $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ rdsrghdr.f
rdstcy.o  : $(MODSTCY) 
	$(FC) $(FFLAGS) -c -o $@ rdstcy.f
rdtprof.o : $(MODTMPRL)
	$(FC) $(FFLAGS) -c -o $@ rdtprof.f
rdtref.o  : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdtref.f
rdvmix.o  : $(MODMOBIL) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdvmix.f
tmprinfo.o  : $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ tmprinfo.f
wrchrscc.o: $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ wrchrscc.f
wridaout.o: $(MODSOURC) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ wridaout.f
xreftbl.o: $(MODXREF) $(MODINFO) $(MODUNCERT)
	$(FC) $(FFLAGS) -c -o $@ xreftbl.f

clean:
	/bin/rm *.o
