# Version "@(#)$Id$ $Source$ $Date$ 
#
all: debug 
#
debug: rawpoint.debug grdpmat.debug
#
FC=f90
FFLAGS = $(DBGFLAGS)
LDFLAGS = -g
LDIRS = -L$(EDSS_LIB)  \
        -L$(EL_SRC)/$(EDSS_OS)
LIBS=$(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH= $(EL_SRC)/$(EDSS_OS)/libsmoke.debug.a
#

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
RAW=rawpoint.o rdptinv.o rdemspt.o fixstk.o openpnts.o wpntschr.o \
    wpntspol.o wrptref.o
#
rawpoint.debug: $(RAW) $(ALLI) $(LIBPATH)
	echo RAWPOINT:
	$(FC) $(FFLAGS) -o $(@) $(RAW) $(LIBS)
#
GRD=grdpmat.o genpgmat.o rpntschr.o
#
grdpmat.debug: $(GRD) $(ALLI) $(LIBPATH)
	echo GRDPMAT:
	$(FC) $(FFLAGS) -o $(@) $(GRD) $(LIBS)
#
#
TMP=tmppoint.o openptmp.o rpntschr.o
#
tmppoint.debug: $(TMP) $(ALLI) $(LIBPATH)
	echo TMPPOINT:
	$(FC) $(FFLAGS) -o $(@) $(TMP) $(LIBS)
#
TST=testprog.o rdsipols.o
#
TSTI=EMCNST3.EXT
#
testprog.debug: $(TST) $(TSTI)
	echo TESTPROG:
	$(FC) $(FFLAGS) -o $(@) $(TST) $(LIBS)
#
# Modules
#
MODSOURC = $(MD_SRC)/modsourc.o
MODXREF  = $(MD_SRC)/modxref.o
MODTPRO  = $(MD_SRC)/modtpro.o
MODSPRO  = $(MD_SRC)/modspro.o
#
# Module dependencies
#
$(MODSOURC):
	$(FC) $(FFLAGS) -c -o $@ $<
	mv modsourc.M $(MD_SRC)

$(MODXREF):
	$(FC) $(FFLAGS) -c -o $@ $<
	mv modxref.M $(MD_SRC)

$(MODTPRO):
	$(FC) $(FFLAGS) -c -o $@ $<
	mv modtpro.M $(MD_SRC)

$(MODSPRO):
	$(FC) $(FFLAGS) -c -o $@ $<
	mv modspro.M $(MD_SRC)

tmppoint.o: $(MODSOURC) $(MODXREF) $(MODTPRO)
	$(FC) $(FFLAGS) -c -o $@ $<
fixstk.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
grdpmat.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
rawpoint.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
rdemspt.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
rdidapt.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
rdptinv.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
rpntschr.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
wpntschr.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
#
clean:
	rm *.o
