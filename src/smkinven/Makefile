# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS)

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
	$(FC) $(FFLAGS) -L$(FSDIR) -o $(SMK_BIN)/$(@) $(GRW) $(LIBS)
#
#
grwinven.debug: $(GRW) $(ALLI) $(LIBPATH_DEBUG)
	echo GRWINVEN \(debug\):
	$(FC) $(FFLAGS) -L$(FSDIR) -o $(@) $(GRW) $(LIBS_DEBUG)

RAW=adjustinv.o asgnar2pt.o asgnnhapx.o chklstfl.o fixstk.o genpdout.o  \
    gethdr.o getpdinfo.o initinfo.o \
    openinvin.o openinvout.o openpdout.o procinven.o rdcempd.o rdemsar.o \
    rdemsmb.o rdemspd.o rdemspt.o rdepsar.o rdepspd.o rdepspt.o rdgrdapi.o \
    rdidaar.o rdidamb.o rdidapt.o rdinven.o rdlooppd.o rdntiar.o rdntimb.o smkinven.o \
    srcmem.o wrinvchr.o wrinvemis.o wrptref.o wrpdemis.o wrepinven.o \
    setnonhap.o
#    
#
smkinven: $(RAW) $(ALLI) $(LIBPATH)
	echo SMKINVEN:
	$(FC) $(FFLAGS) -L$(FSDIR) -o $(SMK_BIN)/$(@) $(RAW) $(LIBS)
#
#
smkinven.debug: $(RAW) $(ALLI) $(LIBPATH_DEBUG)
	echo SMKINVEN \(debug\):
	$(FC) $(FFLAGS) -L$(FSDIR) -o $(@) $(RAW) $(LIBS_DEBUG)

include $(MD_SRC)/modmake.inc

adjustinv.o: $(MODSOURC) $(MODXREF) $(MODINFO) $(MODAR2PT)
	$(FC) $(FFLAGS) -c -o $@ adjustinv.f
asgnar2pt.o: $(MODSOURC) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ asgnar2pt.f
asgnnhapx.o: $(MODSOURC) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ asgnnhapx.f
fixstk.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ fixstk.f
genpdout.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ genpdout.f
gethdr.o: $(MODINFO) $(MODLISTS) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ gethdr.f
getpdinfo.o: $(MODINFO) $(MODDAYHR)
	$(FC) $(FFLAGS) -c -o $@ getpdinfo.f
grwinven.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ grwinven.f
initinfo.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ initinfo.f
opengrwout.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ opengrwout.f
openinvout.o: $(MODLISTS) $(MODINFO) $(MODFILESET)
	$(FC) $(FFLAGS) -I$(FSDIR) -c -o $@ openinvout.f
openpdout.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ openpdout.f
procinven.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ procinven.f
rdcempd.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
	$(FC) $(FFLAGS) -c -o $@ rdcempd.f
rdemsar.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdemsar.f
rdemsmb.o: $(MODSOURC) $(MODMOBIL) $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdemsmb.f
rdemspd.o: $(MODSOURC) $(MODINFO) $(MODDAYHR) 
	$(FC) $(FFLAGS) -c -o $@ rdemspd.f
rdemspt.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdemspt.f
rdepsar.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdepsar.f
rdepspd.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ rdepspd.f
rdepspt.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdepspt.f
rdgrdapi.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdgrdapi.f
rdidaar.o: $(MODSOURC) $(MODLISTS) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdidaar.f
rdidamb.o: $(MODSOURC) $(MODMOBIL) $(MODLISTS) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdidamb.f
rdidapt.o: $(MODSOURC) $(MODLISTS) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdidapt.f
rdinven.o: $(MODSOURC) $(MODLISTS) $(MODMOBIL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdinven.f
rdlooppd.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdlooppd.f
rdntiar.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdntiar.f
rdntimb.o: $(MODSOURC) $(MODMOBIL) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdntimb.f
smkinven.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ smkinven.f
srcmem.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ srcmem.f
wrinvchr.o: $(MODSOURC) $(MODINFO) $(MODFILESET)
	$(FC) $(FFLAGS) -I$(FSDIR) -c -o $@ wrinvchr.f
wrinvemis.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ wrinvemis.f
wrinvpol.o: $(MODFILESET)
	$(FC) $(FFLAGS) -I$(FSDIR) -c -o $@ wrinvpol.f
wrpdemis.o: $(MODSOURC) $(MODINFO) $(MODDAYHR) 
	$(FC) $(FFLAGS) -c -o $@ wrpdemis.f
wrptref.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ wrptref.f
wrepinven.o: $(MODINFO) $(MODLISTS) $(MODAR2PT)
	$(FC) $(FFLAGS) -c -o $@ wrepinven.f

clean:
	/bin/rm *.o
