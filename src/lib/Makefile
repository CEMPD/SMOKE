# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 

# Makefile for SMOKE library

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc

FILES_A = \
  applumat.o   fltrxref.o   getrfdsc.o   procspro.o   rdsrghdr.o \
  alocchrt.o   chkisiz.o    fmtcsrc.o    getsinfo.o   progdesc.o   rdtprof.o  \
  alocctbl.o   chkmetem.o   gettzone.o   rdtref.o     \
  alocetbl.o   chkptdef.o   grd2cnty.o   rdascc.o     rdtzone.o  \
  alocgtbl.o   genuslst.o   initem.o     rdchrscc.o   readwr3.o  \
  alocmtbl.o   filletbl.o   getbasyr.o   rdsrg.o      rdsrghdr.o

FILES_B = \
  alocstbl.o   chksrcno.o   rdstcy.o     rdgmat.o     rdinvchr.o \
  alocttbl.o   convrtxy.o   getctgry.o   mbsccadj.o   rdmvinfo.o     \
  applreac.o   rdeproc.o    rdinvpol.o   rdhdays.o \
  asgnpsi.o    dscsprof.o   getformt.o   rdlines.o    verchar.o  \
  efsetup.o    getidasz.o   normtpro.o   rdpelv.o     tmprinfo.o \
               ingrid.o     rdrmat.o    

FILES_C = \
  bldcsrc.o    fillchrt.o   getiname.o   parscsrc.o   rdsconv.o    wrchrscc.o \
  bldenams.o   fillctbl.o   getinvyr.o   rdcodnam.o   wrdaymsg.o   wridaout.o \
  fillgtbl.o   getisize.o   pdsetup.o    rdsmat.o     rdumat.o     wridapol.o \
  fillmtbl.o   rdsprof.o    \
  chkcpvar.o   fillstbl.o   getvmix.o    rdgref.o     rdvmix.o     xreftbl.o  \
  fillttbl.o   rdefxref.o   rdsref.o     genptcel.o   chkemepi.o

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
	$(FC) $(FFLAGS) -c -o $@ $< 
alocctbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ $<
alocetbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ $<
alocgtbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ $<
alocmtbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ $<
alocstbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ $< 
alocttbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ $< 
asgnpsi.o  : $(MODSOURC) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
bldenams.o  : $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
chkptdef.o: $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ $<
dscsprof.o: $(MODSPRO) 
	$(FC) $(FFLAGS) -c -o $@ $<     
efsetup.o : $(MODEMFAC) $(INC)/M5CNST3.EXT
	$(FC) $(FFLAGS) -c -o $@ $<
fillchrt.o: $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
fillctbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ $<
filletbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ $<
fillgtbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ $<
fillmtbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ $<
fillstbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ $<
fillttbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ $<
fltrxref.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
fmtcsrc.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
genuslst.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
getbasyr.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
getctgry.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
getformt.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
getsinfo.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
gettzone.o: $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ $<
getvmix.o: $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
grd2cnty.o: $(MODSTCY) $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ $<
lambert.o: 
	$(FC) -static $(FFLAGS) -c -o $@ $<
mbsccadj.o  : $(MODMOBIL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<     
normtpro.o: $(MODTMPRL)
	$(FC) $(FFLAGS) -c -o $@ $<     
procspro.o: $(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdcodnam.o: $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ $<
rdefxref.o: $(MODMOBIL) $(MODXREF) $(MODEMFAC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdeproc.o: $(MODEMFAC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdgref.o : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdinvchr.o: $(MODSOURC) 
	$(FC) $(FFLAGS) -c -o $@ $<
rdmvinfo.o: $(MODMOBIL)
	$(FC) $(FFLAGS) -c -o $@ $<
rdpelv.o : $(MODINFO) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ $<
rdsconv.o : $(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdsref.o  : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdsrg.o   : $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ $<  
rdstcy.o  : $(MODSTCY) 
	$(FC) $(FFLAGS) -c -o $@ $<  
rdtprof.o : $(MODTMPRL)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdtref.o  : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<     
rdvmix.o  : $(MODMOBIL) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<     
tmprinfo.o  : $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
wrchrscc.o: $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
wridaout.o: $(MODSOURC) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
xreftbl.o: $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	/bin/rm *.o
