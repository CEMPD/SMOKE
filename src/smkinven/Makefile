# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a

all: smkinven

debug: smkinven.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT
#
RAW=chklstfl.o fixstk.o genpdout.o gethdr.o getpdinfo.o initinfo.o openinvin.o \
    openinvout.o openpdout.o procinven.o rdemsar.o rdemsmb.o rdemspd.o \
    rdemspt.o rdepsar.o rdepspd.o rdepspt.o rdidaar.o rdidamb.o rdidapt.o \
    rdinven.o rdlooppd.o smkinven.o srcmem.o wrinvchr.o wrinvemis.o wrinvpol.o \
    wrptref.o wrpdemis.o 
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
genpdout.o: $(MODSOURC) $(MODINFO) $(MODDAYHR) 
	$(FC) $(FFLAGS) -c -o $@ $<
genpdinfo.o: $(MODINFO) $(MODDAYHR)
	$(FC) $(FFLAGS) -c -o $@ $<
gethdr.o: $(MODINFO) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ $<
initinfo.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
openinvout.o: $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
openpdout.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
procinven.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdemsar.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdemsmb.o: $(MODSOURC) $(MODMOBIL) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdemspd.o: $(MODSOURC) $(MODINFO) $(MODDAYHR)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdemspt.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdepspd.o: $(MODSOURC) $(MODINFO) $(MODDAYHR) $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdepspt.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdidaar.o: $(MODSOURC) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdidamb.o: $(MODSOURC) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdidapt.o: $(MODSOURC) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdinven.o: $(MODSOURC) $(MODMOBIL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdlooppd.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
smkinven.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
srcmem.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
wrinvchr.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
wrinvemis.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
wrptref.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
wrpdemis.o: $(MODSOURC) $(MODINFO) $(MODDAYHR) 
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	/bin/rm *.o
