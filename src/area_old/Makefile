# Version "@(#)$Id$ $Source$ $Date$ 
#
all: opts
#
opts: csgarea ctlamat emsarea grdamat rawarea spcamat tmparea 
#
LIBS=$(LDIRS) -lemstuff -lm3io -lnetcdf
#
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
#
CTLS=ctlamat.F
#
CTLI=ARDIMS3.EXT CHDIMS3.EXT
#
ctlamat:  $(CTLS) $(CTLI)
	echo CTLAMAT:
	f77 $(FFLAGS) -o $(SCENDIR)/ctlamat $(CTLS) $(LIBS)
#
#
#
EMS=emsarea.F
#
EMSI=ARDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT TMDIMS3.EXT
#
emsarea: $(EMS) $(EMSI)
	echo EMSAREA:
	f77 $(FFLAGS) -o $(SCENDIR)/emsarea $(EMS) $(LIBS)
#
#
#
GRMS=grdamat.F
#
GRMI=ARDIMS3.EXT GRDIMS3.EXT CHDIMS3.EXT
#
grdamat: $(GRMS) $(GRMI)
	echo GRDAMAT:
	f77 $(FFLAGS) -o $(SCENDIR)/grdamat $(GRMS) $(LIBS)
#
#
#
RAWS=rawarea.F
#
RAWI=ARDIMS3.EXT GRDIMS3.EXT CHDIMS3.EXT TMDIMS3.EXT
#
rawarea: $(RAWS) $(RAWI)
	echo RAWAREA:
	f77 $(FFLAGS) -o $(SCENDIR)/rawarea $(RAWS) $(LIBS)
#
#
#
SPCS=spcamat.F rdaconv.F rdsprof.F rdascc.F
#
SPCI=ARDIMS3.EXT CHDIMS3.EXT
#
spcamat:  $(SPCS) $(SPCI)
	echo SPCAMAT:
	f77 $(FFLAGS) -o $(SCENDIR)/spcamat $(SPCS) $(LIBS)
#
#
#
TMPS=tmparea.F mktmatsi.F rdascc.F
#
TMPI=ARDIMS3.EXT CHDIMS3.EXT TMDIMS3.EXT
#
tmparea: $(TMPS) $(TMPI)
	echo TMPAREA:
	f77 $(FFLAGS) -o $(SCENDIR)/tmparea $(TMPS) $(LIBS) 
#
#
#
clean:
	rm *.o
