# Version "@(#)$Id$ $Source$ $Date$ 
#
all: opts
#
opts: rawpoint emspoint pdpoint elevpoint grdpmat spcpmat \
      ctlpmat tmppoint csgpoint cselev laypoint
#
LDIRS=-L$(em) \
      -L$(io) \
      -L$(ncf)
#
LIBS=$(LDIRS) -lemstuff -lm3io -lnetcdf
#
#
RAW=rawpoint.F
#
RAWI=PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT
#
rawpoint: $(RAW) $(RAWI)
	echo RAWPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/rawpoint $(RAW) $(LIBS)
#
#
#
EMS=emspoint.F
#
EMSI=PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT TMDIMS3.EXT 
#
emspoint: $(EMS) $(EMSI)
	echo EMSPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/emspoint $(EMS) $(LIBS)
#
#
#
DAYP=pdpoint.F dayrec.F
#
pdpoint: $(DAYP) $(RAWI)
	echo PDPOINT:
	f77 $(FFLAGS) -o $(SCENDIR)/pdpoint $(DAYP) $(LIBS)
#
#
#
ELEV=elevpoint.F plumris.F
#
elevpoint: $(ELEV) $(RAWI)
	echo ELEVPOINT:
	f77 $(FFLAGS) -o $(SCENDIR)/elevpoint $(ELEV) $(LIBS)
#
#
#
GRD=grdpmat.F
#
GRDI=PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT
#
grdpmat: $(GRD) $(GRDI)
	echo GRDPMAT:
	f77 $(FFLAGS)   -o $(SCENDIR)/grdpmat $(GRD) $(LIBS)
#
#
#
SPC=spcpmat.F rdsprof.F rdpconv.F
#
SPCI=PTDIMS3.EXT CHDIMS3.EXT
#
spcpmat: $(SPCI) $(SPC)
	echo SPCPMAT:
	f77 $(FFLAGS)   -o $(SCENDIR)/spcpmat $(SPC) $(LIBS)
#
#
#
TMP=tmppoint.F mktmat.F
#
TMPI=PTDIMS3.EXT CHDIMS3.EXT
#
tmppoint: $(TMP) $(TMPI)
	echo TMPPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/tmppoint $(TMP) $(LIBS)
#
#
#
CTLS=ctlpmat.F
#
CTLI=PTDIMS3.EXT CHDIMS3.EXT
#
ctlpmat:  $(CTLS) $(CTLI)
	echo CTLPMAT:
	f77 $(FFLAGS) -o $(SCENDIR)/ctlpmat $(CTLS) $(LIBS)
#
#
#
CSG=csgpoint.F
#
CSGI=PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT
#
csgpoint: $(CSG) $(CSGI)
	echo CSGPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/csgpoint $(CSG) $(LIBS)
#
#
#
CSE=cselev.F
#
CSGJ=PTDIMS3.EXT CHDIMS3.EXT
#
cselev: $(CSE) $(CSGJ)
	echo CSELEV:
	f77 $(FFLAGS)   -o $(SCENDIR)/cselev $(CSE) $(LIBS)
#
#
#
LAY=laypoint.F plmris.F fixstk.F
#
LAYE=PTDIMS3.EXT TMDIMS3.EXT LADIMS3.EXT GRDIMS3.EXT
#
laypoint: $(LAY) $(LAYE)
	echo LAYPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/laypoint $(LAY) $(LIBS)
#
#
clean:
	rm *.o
