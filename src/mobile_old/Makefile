# Version "@(#)$Id$ $Source$ $Date$ 
#
all: opts
#
opts: csgmobil emisfac grdmmat grdmobil prediur rawmobil spcmmat tmpmobil
#
LIBS=$(LDIRS) -lemstuff -lm3io -lnetcdf 
#
LMOB=$(LDIRS) -lemstuff -lm3io -lnetcdf -lmo5av1
#
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
#
#
EMIS=emisfac.F efcombo.F moinbld.F mobile.F

EMISI=MBDIMS3.EXT CHDGEN3.EXT

emisfac: $(EMIS) $(EMISI)
	echo EMISFAC:
	f77 $(FFLAGS) -o $(SCENDIR)/emisfac $(EMIS) $(LMOB)
#
#
#
GRMS=grdmmat.F lnk2grd.F
#
GRMI=MBDIMS3.EXT GRDIMS3.EXT
#
grdmmat: $(GRMS) $(GRMI)
	echo GRDMMAT:
	f77 $(FFLAGS) -o $(SCENDIR)/grdmmat $(GRMS) $(LIBS)
#
#
#
GRDS=grdmobil.F
#
GRDI=MBDIMS3.EXT GRDIMS3.EXT
#
grdmobil: $(GRD) $(GRDI)
	echo GRDMOBIL:
	f77 $(FFLAGS) -o $(SCENDIR)/grdmobil $(GRDS) $(LIBS)
#
#
#
PRED=prediur.F rdmplist.F

PREDI=MBDIMS3.EXT GRDIMS3.EXT

prediur: $(PRED) $(PREDI)
	echo PREDIUR:
	f77 $(FFLAGS) -o $(SCENDIR)/prediur $(PRED) $(LIBS)
#
#
#
RAW=rawmobil.F 
#
RAWI=MBDIMS3.EXT TMDIMS3.EXT
#
rawmobil: $(RAW) $(RAWI)
	echo RAWMOBIL:
	f77 $(FFLAGS)  -o $(SCENDIR)/rawmobil $(RAW) $(LIBS)
#
#
#
SPCS=spcmmat.F csvname.F rdsprof.F rdmconv.F
#
SPCI=MBDIMS3.EXT CHDGEN3.EXT
#
spcmmat: $(SPC) $(SPCI)
	echo SPCMMAT:
	f77 $(FFLAGS) -o $(SCENDIR)/spcmmat $(SPCS) $(LIBS)
#
#
#
TMPV=tmpmobil.F mkvmat.F rdmplist.F csvname.F rdvmix.F
#
TMPI=MBDIMS3.EXT CHDGEN3.EXT TMDIMS3.EXT GRDIMS3.EXT
#
tmpmobil: $(TMPV) $(TMPI)
	echo TMPMOBIL:
	f77 $(FFLAGS) -o $(SCENDIR)/tmpmobil $(TMPV) $(LIBS)
#
#
#
clean:
	rm *.o
