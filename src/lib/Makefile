# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 

# Makefile for SMOKE library

EDSS_LIB = ./$(EDSS_OS)
FC = f90
FFLAGS = $(DBGFLAGS)

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

include $(MD_SRC)/modmake.inc

FILES_A = \
               chkint.o     fltrxref.o   getsinfo.o   procspro.o   rdsrghdr.o \
  alocchrt.o   chkisiz.o    fmtcsrc.o    gettzone.o   progdesc.o   rdtprof.o  \
  alocctbl.o   chkmetem.o   fltrneg.o    hdrmiss3.o   rd3mask.o    rdtref.o   \
  alocetbl.o   chkptdef.o   grd2cnty.o   initem.o     rdascc.o     rdtzone.o  \
  alocgtbl.o   chkreal.o    genuslst.o   lambert.o    rdchrscc.o   readwr3.o  \
               filletbl.o   getbasyr.o   rdstcy.o     rdsrg.o      rdsrghdr.o

FILES_B = \
  alocstbl.o   chksrcno.o   getcfdsc.o   mktmat.o     rdgmat.o     rdinvchr.o \
  alocttbl.o   convrtxy.o   getctgry.o   normtpro.o   rdmvinfo.o   setvptr.o  \
  applreac.o   dscm3grd.o   getfline.o                rdinvpol.o   verchar.o  \
               dscsprof.o   getformt.o                rdlines.o    wrchrscc.o \
               genptcel.o   getidasz.o   padzero.o    rdpelv.o     wrdaymsg.o \
                            getifdsc.o   parscsrc.o   rdrmat.o

FILES_C = \
  bldcsrc.o    fillchrt.o   getiname.o   parsline.o   rdsconv.o    rdumat.o \
  bldenams.o   fillctbl.o   getinvyr.o   pdsetup.o    rdcodnam.o   wrinvpol.o \
  checkmem.o   fillgtbl.o   getisize.o   polmesg.o    rdsmat.o     xreftbl.o  \
  chkcpvar.o   fillstbl.o   getm3epi.o   rdgref.o     rdsprof.o               \
  chkemepi.o   fillttbl.o   getnlist.o   rdefxref.o   rdsref.o 

INCLUDES = $(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

# NOTE: A some point, insert dependencies on include files

lib: $(MODALL) $(FILES_A) $(FILES_B) $(FILES_C) Makefile
	ar rv $(EDSS_LIB)/libsmoke.a $(FILES_A) 
	ar rv $(EDSS_LIB)/libsmoke.a $(FILES_B) 
	ar rv $(EDSS_LIB)/libsmoke.a $(FILES_C) $(MODALL)

debug: $(MODALL) $(FILES_A) $(FILES_B) $(FILES_C) Makefile
	ar rv $(EDSS_LIB)/libsmoke.debug.a $(FILES_A) 
	ar rv $(EDSS_LIB)/libsmoke.debug.a $(FILES_B) 
	ar rv $(EDSS_LIB)/libsmoke.debug.a $(FILES_C) $(MODALL)

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
chkptdef.o: $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ $<     
dscsprof.o: $(MODSPRO) 
	$(FC) $(FFLAGS) -c -o $@ $<     
fillchrt.o: $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
fillctbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ $<
filletbl.o: $(MODXREF) 
	$(FC) $(FFLAGS) -c -o $@ $<
fillgtbl.o: $(MODXREF) 
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
mktmat.o  : $(MODTPRO)
	$(FC) $(FFLAGS) -c -o $@ $<     
normtpro.o: $(MODTPRO)
	$(FC) $(FFLAGS) -c -o $@ $<     
procspro.o: $(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdefxref.o: $(MODXREF) $(MODEMFAC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdgref.o : $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdinvchr.o: $(MODSOURC) 
	$(FC) $(FFLAGS) -c -o $@ $<     
rdmvinfo.o : $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdpelv.o : $(MODELEV)
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
wrchrscc.o: $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ $<
xreftbl.o: $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<     
