# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -ledsstools -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -ledsstools.debug -lioapi.debug -lnetcdf.debug

LIBPATH       = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a

all: grwinven smkinven

debug: grwinven.debug smkinven.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT
#
GRW=grwinven.o opengrwout.o srcmem.o wrinvchr.o wrinvemis.o wrinvpol.o
#    
#
grwinven: $(GRW) $(ALLI) $(LIBPATH)
	echo GRWINVEN:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(GRW) $(LIBS)
#
#
grwinven.debug: $(GRW) $(ALLI) $(LIBPATH_DEBUG)
	echo GRWINVEN \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(GRW) $(LIBS_DEBUG)

RAW=chklstfl.o fixstk.o genpdout.o gethdr.o getpdinfo.o initinfo.o openinvin.o \
    openinvout.o openpdout.o procinven.o rdemsar.o rdemsmb.o rdemspd.o \
    rdemspt.o rdepsar.o rdepspd.o rdepspt.o rdgrdapi.o rdidaar.o rdidamb.o \
    rdidapt.o rdinven.o rdlooppd.o smkinven.o srcmem.o wrinvchr.o wrinvemis.o \
    wrinvpol.o wrptref.o wrpdemis.o 
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
gethdr.o: $(MODINFO) $(MODLISTS) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ $<
grwinven.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
initinfo.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
opengrwout.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
openinvout.o: $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
openpdout.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
procinven.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdemsar.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdemsmb.o: $(MODSOURC) $(MODMOBIL) $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdemspd.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdemspt.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdepsar.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdepspd.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdepspt.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdgrdapi.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdidaar.o: $(MODSOURC) $(MODLISTS) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdidamb.o: $(MODSOURC) $(MODMOBIL) $(MODLISTS) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdidapt.o: $(MODSOURC) $(MODLISTS) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdinven.o: $(MODSOURC) $(MODLISTS) $(MODMOBIL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdlooppd.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
smkinven.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
srcmem.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
wrinvchr.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
wrinvemis.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
wrptref.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
wrpdemis.o: $(MODSOURC) $(MODINFO) $(MODDAYHR) 
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	/bin/rm *.o
