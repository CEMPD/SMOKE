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
               filletbl.o   getbasyr.o   lambert.o    rdsrg.o      rdsrghdr.o

FILES_B = \
  alocstbl.o   chksrcno.o   getcfdsc.o   rdstcy.o     rdgmat.o     rdinvchr.o \
  alocttbl.o   convrtxy.o   getctgry.o   mbsccadj.o   rdmvinfo.o   setvptr.o  \
  applreac.o   dscm3grd.o   getfline.o   mktmat.o     rdeproc.o    rdinvpol.o \
  asgnpsi.o    dscsprof.o   getformt.o   multunit.o   rdlines.o    verchar.o  \
               efsetup.o    getidasz.o   normtpro.o   rdpelv.o     tmprinfo.o \
                            getifdsc.o   padzero.o    rdrmat.o     unitfac.o

FILES_C = \
  bldcsrc.o                              parscsrc.o                wrchrscc.o \
  bldenams.o   fillchrt.o   getiname.o   parsline.o   rdsconv.o    wrdaymsg.o \
  chkgrid.o    fillctbl.o   getinvyr.o   pdsetup.o    rdcodnam.o   rdumat.o   \
  checkmem.o   fillgtbl.o   getisize.o   polmesg.o    rdsmat.o     wrinvpol.o \
  chkcpvar.o   fillstbl.o   getm3epi.o   rdgref.o     rdsprof.o    xreftbl.o  \
  fillttbl.o   getnlist.o   rdefxref.o   rdsref.o     genptcel.o   chkemepi.o

INCLUDES = $(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

# NOTE: A some point, insert dependencies on include files

all: lib

lib: $(MODALL) $(FILES_A) $(FILES_B) $(FILES_C) Makefile
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_A) 
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_B) 
	ar rv $(SMK_BIN)/libsmoke.a $(FILES_C)
	ar rv $(SMK_BIN)/libsmoke.a $(MODALL)

debug: $(MODALL) $(FILES_A) $(FILES_B) $(FILES_C) Makefile
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_A) 
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_B) 
	ar rv $(SMK_BIN)/libsmoke.debug.a $(FILES_C) 
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
alocstbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ $< 
alocttbl.o: $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ $< 
asgnpsi.o  : $(MODSOURC) $(MODXREF) $(MODINFO)
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
fillgtbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ $<
fillctbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ $<
filletbl.o: $(MODXREF) 
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
grd2cnty.o: $(MODSTCY) $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ $<
mbsccadj.o  : $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<     
mktmat.o  : $(MODTPRO)
	$(FC) $(FFLAGS) -c -o $@ $<     
normtpro.o: $(MODTPRO)
	$(FC) $(FFLAGS) -c -o $@ $<     
procspro.o: $(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdefxref.o: $(MODXREF) $(MODEMFAC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdeproc.o: $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ $<
rdgref.o : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdinvchr.o: $(MODSOURC) 
	$(FC) $(FFLAGS) -c -o $@ $<
rdmvinfo.o : $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdpelv.o : $(MODINFO) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ $<
rdsconv.o : $(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdsref.o  : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdsrg.o : $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ $<  
rdstcy.o : $(MODSTCY) 
	$(FC) $(FFLAGS) -c -o $@ $<  
rdtprof.o  : $(MODTPRO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdtref.o  : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<     
tmprinfo.o  : $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
wrchrscc.o: $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ $<
xreftbl.o: $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	/bin/rm *.o
