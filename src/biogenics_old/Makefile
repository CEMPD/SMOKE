# Version "@(#)$Id$ $Source$ $Date$ 
#
all: opts
#
opts: grdbio rawbio tmpbio
#
LIBS=$(LDIRS) -lemlib -lioapi -lnetcdf
#
#
INCS=BIODIMS3.EXT GRDIMS3.EXT TMDIMS3.EXT
#
#
#
GRDS=grdbio.F
#
grdbio:  $(INCS) $(GRDS)
	f77 $(FFLAGS) -o $(SCENDIR)/grdbio $(GRDS) $(LIBS)
#
#
#
RAWS=rawbio.F
#
rawbio:  $(INCS) $(RAWS)
	f77 $(FFLAGS) -o $(SCENDIR)/rawbio $(RAWS) $(LIBS)
#
#
#
TMPS=tmpbio.F hrbio.F hrbios.F cldaten.F cldaten_calc.F czangle.F
#
tmpbio:  $(INCS) $(TMPS)
	f77 $(FFLAGS) -o $(SCENDIR)/tmpbio $(TMPS) $(LIBS)
#
#
#
