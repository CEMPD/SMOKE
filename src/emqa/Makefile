# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)
LIBPATH       = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a

all: smkreport
#
debug: smkreport.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

SMR=smkreport.o asgnbins.o bldrepidx.o genrprt.o openrepin.o   \
    openrepout.o qarepin.o rdgrps.o rdrepin.o rdrprts.o repmrggrd.o repunits.o\
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
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

asgnbins.o: $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODELEV) \
            $(MODSTCY) $(MODINFO) $(MODGRID)
bldrepidx.o: $(MODREPRT) $(MODREPBN) $(MODTMPRL) $(MODINFO) $(MODFILESET) $(MODCNTRL)
genrprt.o: $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODCNTRL) $(MODINFO) $(MODSTCY) 
openrepin.o: $(MODREPRT) $(MODREPBN) $(MODTMPRL) $(MODINFO) $(MODGRID) $(MODFILESET) $(MODCNTRL)
qarepin.o: $(MODREPRT) $(MODSTCY) $(MODINFO)
rdgrps.o: $(MODSOURC)  $(MODREPRT) $(MODREPBN) $(MODLISTS) $(MODGRID)
rdrepin.o: $(MODSOURC)  $(MODREPRT) $(MODREPBN) $(MODELEV) $(MODLISTS) \
           $(MODINFO) $(MODGRID) $(MODFILESET) $(MODCNTRL)
rdrprts.o: $(MODREPRT) $(MODINFO)
repmrggrd.o: $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODINFO) 
repunits.o: $(MODREPRT) $(MODREPBN)  $(MODTMPRL) $(MODINFO) 
scanrepc.o: $(MODREPRT) $(MODINFO) 
selectsrc.o: $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODINFO) 
smkreport.o: $(MODREPRT) $(MODREPBN) $(MODINFO) $(MODGRID)
wrrephdr.o: $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODSTCY) \
            $(MODINFO) $(MODGRID)
wrrepout.o: $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODSTCY) \
            $(MODINFO) 
#
#
clean:
	/bin/rm -f *.o
