# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -liotools -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -liotools.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libiotools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libiotools.debug.a

all: smkreport
#
debug: smkreport.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT


#
#GTR=getrecs.o
#
#getrecs: $(GTR) $(ALLI) $(LIBPATH)
#	echo GETRECS:
#	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(GTR) $(LIBS)
#
#getrecs.debug: $(GTR) $(ALLI) $(LIBPATH_DEBUG)
#	echo GETRECS \(debug\):
#	$(FC) $(FFLAGS) -o $(@) $(GTR) $(LIBS_DEBUG)
#
#
SMR=smkreport.o asgnbins.o bldrepidx.o blkorcmt.o genrprt.o openrepin.o   \
    openrepout.o prclinrc.o rdrepin.o rdrprts.o repmrggrd.o repunits.o  \
    rmcommnt.o scanrepc.o selectsrc.o wrrephdr.o wrrepout.o

smkreport: $(SMR) $(ALLI) $(LIBPATH)
	echo SMKREPORT:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(SMR) $(LIBS)

smkreport.debug: $(SMR) $(ALLI) $(LIBPATH_DEBUG)
	echo SMKREPORT \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(SMR) $(LIBS_DEBUG)
#
#
include $(MD_SRC)/modmake.inc

asgnbins.o: $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODELEV) \
            $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
bldrepidx.o: $(MODREPRT) $(MODREPBN) $(MODTMPRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
genrprt.o: $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
getrecs.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
openrepin.o: $(MODREPRT) $(MODREPBN) $(MODTMPRL) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ $<
prclinrc.o: $(MODREPRT) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ $<
rdrepin.o: $(MODSOURC)  $(MODREPRT) $(MODREPBN) $(MODELEV) $(MODLISTS) \
           $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdrprts.o: $(MODREPRT)
	$(FC) $(FFLAGS) -c -o $@ $<
repmrggrd.o: $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ $<
repunits.o: $(MODREPRT) $(MODREPBN)  $(MODTMPRL) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ $<
scanrepc.o: $(MODREPRT) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ $<
selectsrc.o: $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ $<
smkreport.o: $(MODREPRT) $(MODREPBN) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ $<
wrrephdr.o: $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODSTCY) \
            $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
wrrepout.o: $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODSTCY) \
            $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ $<
#
#
clean:
	/bin/rm *.o
