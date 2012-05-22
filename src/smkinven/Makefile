# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)

LIBPATH       = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a

all: grwinven smkinven

debug: grwinven.debug smkinven.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<
.F.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT
#
GRW=grwinven.o opengrwout.o srcmem.o wrinvchr.o wrinvpol.o
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

RAW=    adjustinv.o asgnar2pt.o asgnnhapx.o chklstfl.o fixstk.o genpdout.o \
        gethdr.o getpdinfo.o initinfo.o formlist.o \
        openinvin.o openinvout.o openpdout.o procinven.o procinvsrcs.o rdcempd.o \
        rdemspd.o rdorlfr.o rdepspd.o rdgrdapi.o \
        rdinvsrcs.o rdlooppd.o smkinven.o \
        srcmem.o wrinvchr.o wrinvemis.o wrinvpol.o wrptref.o wrpdemis.o wrepinven.o \
        setnonhap.o rdsrcidaar.o rdsrcntiar.o rddataidaar.o rdinvdata.o rddatantiar.o \
        procar2pt.o rdsrcidamb.o rdsrcntimb.o rddataidamb.o rddatantimb.o \
        rdsrcff10ar.o rddataff10ar.o rdsrcff10mb.o rddataff10mb.o rdsrcff10pt.o rddataff10pt.o\
        rdsrcidapt.o rddataidapt.o rdsrcemsmb.o rddataemsmb.o rdsrcemsar.o \
        rddataemsar.o rdsrcemspt.o rddataemspt.o procemspt.o rdff10pd.o \
        rddevemspt.o rdfacemspt.o rdprocemspt.o rdstkemspt.o rdsrcntifr.o\
        rdsrcntipt.o rdsrcntinp.o rddatantipt.o rddatantinp.o rddatantifr.o rdcemsum.o
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
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

asgnar2pt.o: $(MODSOURC) $(MODXREF) $(MODINFO)
asgnnhapx.o: $(MODSOURC) $(MODXREF) $(MODINFO)
fixstk.o: $(MODSOURC)
formlist.o: $(MODINFO) $(MODLISTS)
genpdout.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODSTCY)
gethdr.o: $(MODLISTS) $(MODSTCY) $(MODINFO)
getpdinfo.o: $(MODINFO) $(MODDAYHR)
grwinven.o: $(MODSOURC) $(MODINFO) $(MODFILESET)
initinfo.o: $(MODINFO)
opengrwout.o: $(MODINFO)
openinvout.o: $(MODLISTS) $(MODINFO) $(MODFILESET)
openpdout.o: $(MODINFO)
procar2pt.o: $(MODSOURC) $(MODXREF) $(MODINFO) $(MODAR2PT)
procemspt.o: $(MODINFO) $(MODSOURC) $(MODLISTS)
procinven.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
procinvsrcs.o: $(MODSOURC) $(MODINFO) $(MODMOBIL)
rdcempd.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
rdcemsum.o: $(MODLISTS) $(MODDAYHR)
rddataemsar.o: $(MODINFO) $(MODLISTS)
rddataemsmb.o: $(MODINFO)
rddataemspt.o: $(MODINFO) $(MODLISTS)
rddataidaar.o: $(MODINFO)
rddataidamb.o: $(MODINFO)
rddataidapt.o: $(MODINFO)
rddatantiar.o: $(MODINFO)
rddataff10ar.o: $(MODINFO)
rddatantimb.o: $(MODINFO)
rddataff10mb.o: $(MODINFO)
rddatantinp.o: $(MODINFO)
rddatantipt.o: $(MODINFO)
rddataff10pt.o: $(MODINFO)
rddatantifr.o: $(MODINFO)
rddevemspt.o: $(MODINFO) $(MODSOURC)
rdemspd.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) 
rdff10pd.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) 
rdepspd.o: $(MODSOURC) $(MODLISTS) $(MODSTCY) $(MODINFO) $(MODDAYHR)
rdfacemspt.o: $(MODINFO) $(MODSOURC)
rdgrdapi.o: $(MODSOURC) $(MODINFO)
rdinvdata.o: $(MODSOURC) $(MODINFO) $(MODLISTS) $(MODMOBIL)
rdinvsrcs.o: $(MODSOURC) $(MODINFO) $(MODLISTS) $(MODMOBIL)
rdlooppd.o: $(MODINFO)
rdorlfr.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) 
rdprocemspt.o: $(MODINFO) $(MODSOURC)
rdsrcemsar.o: $(MODLISTS)
rdsrcemsmb.o: $(MODINFO)
rdsrcemspt.o: $(MODLISTS)
rdsrcidaar.o:
rdsrcidamb.o: $(MODINFO)
rdsrcidapt.o:
rdsrcntiar.o: $(MODLISTS)
rdsrcff10ar.o: $(MODLISTS)
rdsrcntimb.o: $(MODLISTS)
rdsrcff10mb.o: $(MODLISTS)
rdsrcntifr.o: $(MODINFO)
rdsrcntinp.o: $(MODLISTS)
rdsrcntipt.o: $(MODLISTS)
rdsrcff10pt.o: $(MODLISTS)
rdstkemspt.o: $(MODINFO) $(MODSOURC)
setnonhap.o: $(MODSOURC) $(MODLISTS) $(MODXREF) $(MODINFO)
smkinven.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
srcmem.o: $(MODSOURC)
wrepinven.o: $(MODSOURC) $(MODINFO) $(MODLISTS) $(MODAR2PT)
wrinvchr.o: $(MODSOURC) $(MODINFO)
wrinvemis.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODFILESET)
wrinvpol.o: $(MODFILESET)
wrpdemis.o: $(MODSOURC) $(MODINFO) $(MODDAYHR) 
wrptref.o: $(MODINFO)

clean:
	/bin/rm -f *.o 
