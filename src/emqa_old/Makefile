# Version "@(#)$Id$ $Source$ $Date$ 
#
# Programs listed in alphabetical order
#
all: opts

opts: aoicpop aoicsum areasum asctsum csgarea csgmobil \
      csgpoint getrecs grdarea grdmobil grdpoint mobilsum mvmtsum plantsum \
      pntsum psccsum repbio repnorm repplht reptpro
#
#
LIBS=$(LDIRS) -lemstuff -lm3io -lnetcdf
#
LMOB=$(LDIRS) -lemstuff -lm3io -lnetcdf -lmo5b
#
#
#
AOP=aoicpop.F
#
AOPI=ARDIMS3.EXT CHDIMS3.EXT
#
aoicpop: $(AOP) $(AOPI)
	echo AOICPOP:
	f77 $(FFLAGS)   -o $(SCENDIR)/aoicpop $(AOP) $(LIBS)
#
#
AOS=aoicsum.F
#
AOSI=ARDIMS3.EXT CHDIMS3.EXT
#
aoicsum: $(AOS) $(AOSI)
	echo AOICSUM:
	f77 $(FFLAGS)   -o $(SCENDIR)/aoicsum $(AOS) $(LIBS)
#
ARS=areasum.F
#
ARSI=ARDIMS3.EXT CHDIMS3.EXT
#
areasum: $(ARS) $(ARSI)
	echo AREASUM:
	f77 $(FFLAGS)   -o $(SCENDIR)/areasum $(ARS) $(LIBS)
#
#
ACS=asctsum.F
#
ACSI=ARDIMS3.EXT CHDIMS3.EXT
#
asctsum: $(ACS) $(ACSI)
	echo ASCTSUM:
	f77 $(FFLAGS)   -o $(SCENDIR)/asctsum $(ACS) $(LIBS)
#
CSGS=csgarea.F
#
CSGI=ARDIMS3.EXT GRDIMS3.EXT CHDIMS3.EXT
#
csgarea: $(CSGS) $(CSGI)
	echo CSGAREA:
	f77 $(FFLAGS) -o $(SCENDIR)/csgarea $(CSGS) $(LIBS)
#
#
CSGV=csgmobil.F msmatvec.F csvname.F
#
CSGI=MBDIMS3.EXT CHDGEN3.EXT GRDIMS3.EXT
#
csgmobil: $(CSGV) $(CSGI)
	echo CSGMOBIL:
	f77 $(FFLAGS) -o $(SCENDIR)/csgmobil $(CSGV) $(LIBS)
#
CSG=csgpoint.F 
#
CSGI=PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT
#
csgpoint: $(CSG) $(CSGI)
	echo CSGPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/csgpoint $(CSG) $(LIBS)
#
OBJS = getrecs.F

INCS = PARMS3.EXT FDESC3.EXT GRDIMS3.EXT IODECL3.EXT PTDIMS3.EXT

getrecs: $(OBJS) $(INCS)
	echo GETRECS:
	f77 $(FFLAGS) -o $(SCENDIR)/getrecs $(OBJS) $(LIBS)
#
GRMS=grdarea.F
#
GRMI=ARDIMS3.EXT GRDIMS3.EXT CHDIMS3.EXT TMDIMS3.EXT
#
grdarea: $(GRMS) $(GRMI)
	echo GRDAREA:
	f77 $(FFLAGS) -o $(SCENDIR)/grdarea $(GRMS) $(LIBS)
#
GRDS=grdmobil.F
#
GRDI=MBDIMS3.EXT GRDIMS3.EXT
#
grdmobil: $(GRDS) $(GRDI)
	echo GRDMOBIL:
	f77 $(FFLAGS) -o $(SCENDIR)/grdmobil $(GRDS) $(LIBS)
#
GRDPT=grdpoint.F 
#
GRDI=PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT
#
grdpoint: $(GRDPT) $(GRDI)
	echo GRDPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/grdpoint $(GRDPT) $(LIBS)
#
MOS=mobilsum.F
#
MOSI=MBDIMS3.EXT CHDIMS3.EXT
#
mobilsum: $(MOS) $(MOSI)
	echo MOBILSUM:
	f77 $(FFLAGS)   -o $(SCENDIR)/mobilsum $(MOS) $(LIBS)
#
#
MVS=mvmtsum.F
#
MVSI=MBDIMS3.EXT CHDIMS3.EXT
#
mvmtsum: $(MVS) $(MVSI)
	echo MVMTSUM:
	f77 $(FFLAGS)   -o $(SCENDIR)/mvmtsum $(MVS) $(LIBS)
#
PLS=plantsum.F
#
PLSI=PTDIMS3.EXT CHDIMS3.EXT
#
plantsum: $(PLS) $(PLSI)
	echo PLANTSUM:
	f77 $(FFLAGS)   -o $(SCENDIR)/plantsum $(PLS) $(LIBS)
#
#
PTS=pntsum.F
#
PTSI=PTDIMS3.EXT CHDIMS3.EXT
#
pntsum: $(PTS) $(PTSI)
	echo PNTSUM:
	f77 $(FFLAGS)   -o $(SCENDIR)/pntsum $(PTS) $(LIBS)
#
PSC=psccsum.F
#
PSCI=PTDIMS3.EXT CHDIMS3.EXT
#
psccsum: $(PSC) $(PSCI)
	echo PSCCSUM:
	f77 $(FFLAGS)   -o $(SCENDIR)/psccsum $(PSC) $(LIBS)
#
BXP=repbio.F
#
BXPI=BIODIMS3.EXT ARDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT PARMS3.EXT IODECL3.EXT \
     FDESC3.EXT
#
repbio: $(BXP) $(BXPI)
	echo REPBIO:
	f77 $(FFLAGS)   -o $(SCENDIR)/repbio $(BXP) $(LIBS)
#
BNP=repnorm.F
#
BNPI=BIODIMS3.EXT ARDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT PARMS3.EXT IODECL3.EXT \
     FDESC3.EXT
#
repnorm: $(BNP) $(BNPI)
	echo REPNORM:
	f77 $(FFLAGS)   -o $(SCENDIR)/repnorm $(BNP) $(LIBS)
#
#
RPH=repplht.F
#
RPHI=PTDIMS3.EXT CHDIMS3.EXT TMDIMS3.EXT VGRD3.EXT
#
repplht: $(RPH) $(RPHI)
	echo REPPLHT:
	f77 $(FFLAGS)   -o $(SCENDIR)/repplht $(RPH) $(LIBS)
#
#
#
RTP=reptpro.F
#
RTPI=PTDIMS3.EXT CHDIMS3.EXT TMDIMS3.EXT
#
reptpro: $(RTP) $(RTPI)
	echo REPTPRO:
	f77 $(FFLAGS)   -o $(SCENDIR)/reptpro $(RTP) $(LIBS)
#
#
clean:
	rm *.o
