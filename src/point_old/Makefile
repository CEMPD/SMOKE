# Version "@(#)$Id$ $Source$ $Date$ 
#
all: opts
#
opts: cselev csgpoint ctlpmat elevpoint emspoint grdpmat \
      grdpoint laypoint pdpoint rawpoint spcpmat tmppoint
#
LIBS=$(LDIRS) -lemstuff -lm3io -lnetcdf
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
ELEV=elevpoint.F plumris.F
#
ELVI=PTDIMS3.EXT
#
elevpoint: $(ELEV) $(ELVI)
	echo ELEVPOINT:
	f77 $(FFLAGS) -o $(SCENDIR)/elevpoint $(ELEV) $(LIBS)
#
#
#
EMS=emspoint.F fixstk.F
#
EMSI=PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT TMDIMS3.EXT 
#
emspoint: $(EMS) $(EMSI)
	echo EMSPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/emspoint $(EMS) $(LIBS)
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
GRP=grdpoint.F
#
GRPI=PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT
#
grdpoint: $(GRP) $(GRPI)
	echo GRDPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/grdpoint $(GRP) $(LIBS)
#
#
#
LAY=laypoint.F plmris.F
#
LAYE=PTDIMS3.EXT TMDIMS3.EXT VGRD3.EXT GRDIMS3.EXT
#
laypoint: $(LAY) $(LAYE)
	echo LAYPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/laypoint $(LAY) $(LIBS)
#
#
#
DAYP=pdpoint.F dayrec.F
#
DAYI=PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT TMDIMS3.EXT
#
pdpoint: $(DAYP) $(DAYI)
	echo PDPOINT:
	f77 $(FFLAGS) -o $(SCENDIR)/pdpoint $(DAYP) $(LIBS)
#
#
#
RAW=rawpoint.F fixstk.F
#
RAWI=PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT TMDIMS3.EXT
#
rawpoint: $(RAW) $(RAWI)
	echo RAWPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/rawpoint $(RAW) $(LIBS)
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
TMPI=PTDIMS3.EXT CHDIMS3.EXT TMDIMS3.EXT
#
tmppoint: $(TMP) $(TMPI)
	echo TMPPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/tmppoint $(TMP) $(LIBS)
#
#
#
clean:
	rm *.o
