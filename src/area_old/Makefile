# Version "@(#)$Id$ $Source$ $Date$ 
#
all: opts
#
opts: grdamat spcamat ctlamat rawarea tmparea csgarea
#
LDIRS=-L$(em) \
      -L$(io) \
      -L$(ncf)
#
LIBS=$(LDIRS) -lemstuff -lm3io -lnetcdf
#
#
GRMS=grdamat.F
#
GRMI=ARDIMS3.EXT GRDIMS3.EXT CHDIMS3.EXT
#
grdamat: $(GRMS) $(GRMI)
	echo GRDAMAT:
	f77 $(FFLAGS) -o $(DISTDIR)/grdamat $(GRMS) $(LIBS)
#
#
#
SPCS=spcamat.F
#
SPCI=ARDIMS3.EXT CHDIMS3.EXT
#
spcamat:  $(SPCS) $(SPCI)
	echo SPCAMAT:
	f77 $(FFLAGS) -o $(DISTDIR)/spcamat $(SPCS) $(LIBS)
#
#
#
CTLS=ctlamat.F
#
CTLI=ARDIMS3.EXT CHDIMS3.EXT
#
ctlamat:  $(CTLS) $(CTLI)
	echo CTLAMAT:
	f77 $(FFLAGS) -o $(DISTDIR)/ctlamat $(CTLS) $(LIBS)
#
#
#
RAWS=rawarea.F
#
RAWI=ARDIMS3.EXT GRDIMS3.EXT CHDIMS3.EXT
#
rawarea: $(RAWS) $(RAWI)
	echo RAWAREA:
	f77 $(FFLAGS) -o $(DISTDIR)/rawarea $(RAWS) $(LIBS)
#
#
#
TMPS=tmparea.F mktmat.F
#
TMPI=ARDIMS3.EXT CHDIMS3.EXT  TMDIMS3.EXT
#
tmparea: $(TMPS) $(TMPI)
	echo TMPAREA:
	f77 $(FFLAGS) -o $(DISTDIR)/tmparea $(TMPS) $(LIBS) 
#
#
#
CSGS=csgarea.F
#
CSGI=ARDIMS3.EXT GRDIMS3.EXT CHDIMS3.EXT
#
csgarea: $(CSGS) $(CSGI)
	echo CSGAREA:
	f77 $(FFLAGS) -o $(DISTDIR)/csgarea $(CSGS) $(LIBS)
#
#
#
clean:
	rm *.o
