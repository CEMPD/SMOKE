# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT
#
RAW=smkinven.o rdinven.o rdidaar.o rdidamb.o rdidapt.o rdemspt.o fixstk.o \
    procinven.o openinvout.o srcmem.o initinfo.o chklstfl.o wrinvchr.o \
    wrinvpol.o wrptref.o openinvin.o
#    
#
smkinven: $(RAW) $(ALLI) $(LIBPATH)
	echo SMKINVEN:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(RAW) $(LIBS)
#
#
smkinven.debug: $(RAW) $(ALLI) $(LIBPATH_DEBUG)
	echo SMKINVEN \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(RAW) $(LIBS_DEBUG)

include $(MD_SRC)/modmake.inc

fixstk.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
initinfo.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
openinvout.o: $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
procinven.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdemspt.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdidaar.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdidamb.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdidapt.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdinven.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
smkinven.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
srcmem.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
wrinvchr.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
#
clean:
	rm *.o
