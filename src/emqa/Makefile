# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -ledsstools -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -ledsstools.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a

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
SMR=smkreport.o asgnbins.o bldrepidx.o genrprt.o openrepin.o   \
    openrepout.o rdgrps.o rdrepin.o rdrprts.o repmrggrd.o repunits.o\
    scanrepc.o selectsrc.o wrrephdr.o wrrepout.o

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
            $(MODSTCY) $(MODINFO) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ asgnbins.f
bldrepidx.o: $(MODREPRT) $(MODREPBN) $(MODTMPRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ bldrepidx.f
genrprt.o: $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODINFO) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ genrprt.f
getrecs.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ getrecs.f
openrepin.o: $(MODREPRT) $(MODREPBN) $(MODTMPRL) $(MODINFO) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ openrepin.f
rdgrps.o: $(MODSOURC)  $(MODREPRT) $(MODREPBN) $(MODLISTS) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ rdgrps.f
rdrepin.o: $(MODSOURC)  $(MODREPRT) $(MODREPBN) $(MODELEV) $(MODLISTS) \
           $(MODINFO) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ rdrepin.f
rdrprts.o: $(MODREPRT) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdrprts.f
repmrggrd.o: $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ repmrggrd.f
repunits.o: $(MODREPRT) $(MODREPBN)  $(MODTMPRL) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ repunits.f
scanrepc.o: $(MODREPRT) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ scanrepc.f
selectsrc.o: $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ selectsrc.f
smkreport.o: $(MODREPRT) $(MODREPBN) $(MODINFO) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ smkreport.f
wrrephdr.o: $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODSTCY) \
            $(MODINFO) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ wrrephdr.f
wrrepout.o: $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODSTCY) \
            $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ wrrepout.f
#
#
clean:
	/bin/rm *.o
