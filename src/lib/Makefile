# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 

# Makefile for SMOKE library

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

include $(MD_SRC)/modmake.inc

FILES_A = \
  applumat.o   chkint.o     fltrxref.o   getrfdsc.o   procspro.o   rdsrghdr.o \
  alocchrt.o   chkisiz.o    fmtcsrc.o    getsinfo.o   progdesc.o   rdtprof.o  \
  alocctbl.o   chkmetem.o   fltrneg.o    gettzone.o   rd3mask.o    rdtref.o   \
  alocetbl.o   chkptdef.o   grd2cnty.o   hdrmiss3.o   rdascc.o     rdtzone.o  \
  alocgtbl.o   chkreal.o    genuslst.o   initem.o     rdchrscc.o   readwr3.o  \
  alocmtbl.o   filletbl.o   getbasyr.o   lambert.o    rdsrg.o      rdsrghdr.o

FILES_B = \
  alocstbl.o   chksrcno.o   getcfdsc.o   rdstcy.o     rdgmat.o     rdinvchr.o \
  alocttbl.o   convrtxy.o   getctgry.o   mbsccadj.o   rdmvinfo.o     \
  applreac.o   dscm3grd.o   getfline.o   rdeproc.o    rdinvpol.o \
  asgnpsi.o    dscsprof.o   getformt.o   multunit.o   rdlines.o    verchar.o  \
  check3.o     efsetup.o    getidasz.o   normtpro.o   rdpelv.o     tmprinfo.o \
                            getifdsc.o   padzero.o    rdrmat.o     unitfac.o

FILES_C = \
  bldcsrc.o    fillchrt.o   getiname.o   parscsrc.o   rdsconv.o    wrchrscc.o \
  bldenams.o   fillctbl.o   getinvyr.o   parsline.o   rdcodnam.o   wrdaymsg.o \
  chkgrid.o    fillgtbl.o   getisize.o   pdsetup.o    rdsmat.o     rdumat.o   \
  checkmem.o   fillmtbl.o   getm3epi.o   polmesg.o    rdsprof.o    wrinvpol.o \
  chkcpvar.o   fillstbl.o   getvmix.o    rdgref.o     rdvmix.o     xreftbl.o  \
  fillttbl.o   getnlist.o   rdefxref.o   rdsref.o     genptcel.o   chkemepi.o

FILES_TMP = 

INCLUDES = $(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

# NOTE: A some point, insert dependencies on include files

all: lib

lib: $(MODALL) $(FILES_A) $(FILES_B) $(FILES_C) $(FILES_TMP) Makefile
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_A) 
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_B) 
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_C)
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_TMP)
	ar rv $(SMK_BIN)/libsmoke.a $(MODALL)

debug: $(MODALL) $(FILES_A) $(FILES_B) $(FILES_C) $(FILES_TMP) Makefile
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
chkgrid.o: $(MODMERGE) $(INC)/FLTERR.EXT
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
normtpro.o: $(MODTPRO)
	$(FC) $(FFLAGS) -c -o $@ $<     
procspro.o: $(MODSPRO)
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
rdtprof.o : $(MODTPRO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdtref.o  : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<     
rdvmix.o  : $(MODMOBIL) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<     
tmprinfo.o  : $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
wrchrscc.o: $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
xreftbl.o: $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	/bin/rm *.o
