# Version "@(#)$Id$ $Source$ $Date$ 
#
all: debug 
#
debug: rawpoint.debug grdpmat.debug
#
FC=f90
FFLAGS = -g 
LDFLAGS = -g
LDIRS = -L$(EDSS_LIB)  \
        -L/env/proj/smoke/models-3/lib/$(EDSS_OS)
LIBS=$(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
#

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

#
RAWP=rawpoint.o rdptinv.o rdemspt.o
#
RAW=ptmodule.o $(RAWP) fixstk.o openpnts.o wpntschr.o wpntspol.o wrptref.o 
#
RAWI=EMCNST3.EXT 
#
$(RAWP):ptmodule.o

$(RAW):$(RAWI)
#
rawpoint.debug: $(RAW) $(RAWI)
	echo RAWPOINT:
	$(FC) $(FFLAGS) -o $(@) $(RAW) $(LIBS)
#
GRD=grdpmat.o genpgmat.o
#
GRDI=EMCNST3.EXT
#
$(GRD):$(GRDI)
#
grdpmat.debug: $(GRD) $(GRDI)
	echo GRDPMAT:
	$(FC) $(FFLAGS) -o $(@) $(GRD) $(LIBS)
#
TST=testprog.o rdsipols.o
#
TSTI=EMCNST3.EXT
#
testprog.debug: $(TST) $(TSTI)
	echo TESTPROG:
	$(FC) $(FFLAGS) -o $(@) $(TST) $(LIBS)
#
#
#
clean:
	rm *.o
