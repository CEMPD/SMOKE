# Version @(#)$Id$ 
# Path	  $Source$
# Date	  $Date$ 
#
LIBS		  = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG	  = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)

LIBPATH		  = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libedsstools.a
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

RAW=    adjustinv.o asgnar2pt.o asgnnhapx.o chklstfl.o fixstk.o genpdout.o \
        gethdr.o getpdinfo.o initinfo.o \
        openinvin.o openinvout.o openpdout.o procinven.o rdcempd.o \
        rdemspd.o rdepspd.o rdgrdapi.o \
        rdinvsrcs.o rdlooppd.o smkinven.o \
        srcmem.o wrinvchr.o wrinvemis.o wrinvpol.o wrptref.o wrpdemis.o wrepinven.o \
        setnonhap.o rdsrcidaar.o rdsrcntiar.o rddataidaar.o rdinvdata.o rddatantiar.o \
        procar2pt.o rdsrcidamb.o rdsrcntimb.o rddataidamb.o rddatantimb.o \
        rdsrcidapt.o rddataidapt.o rdsrcemsmb.o rddataemsmb.o rdsrcemsar.o \
        rddataemsar.o
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
grwinven.o: $(MODSOURC) $(MODINFO) $(MODFILESET)
	$(FC) $(FFLAGS) -c -o $@ grwinven.F
initinfo.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ initinfo.f
opengrwout.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ opengrwout.F
openinvout.o: $(MODLISTS) $(MODINFO) $(MODFILESET)
	$(FC) $(FFLAGS) -c -o $@ openinvout.F
openpdout.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ openpdout.f
procar2pt.o: $(MODSOURC) $(MODXREF) $(MODINFO) $(MODAR2PT)
	$(FC) $(FFLAGS) -c -o $@ procar2pt.f
procinven.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODMOBIL)
	$(FC) $(FFLAGS) -c -o $@ procinven.f
rdcempd.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
	$(FC) $(FFLAGS) -c -o $@ rdcempd.f
rddataemsmb.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rddataemsmb.f
rddataidaar.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rddataidaar.f
rddataidamb.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rddataidamb.f
rddataidapt.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rddataidapt.f
rddatantiar.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rddatantiar.f
rddatantimb.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rddatantimb.f
rdemsar.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdemsar.f
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
rdinvdata.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODMOBIL)
	$(FC) $(FFLAGS) -c -o $@ rdinvdata.f
rdinvsrcs.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODMOBIL)
	$(FC) $(FFLAGS) -c -o $@ rdinvsrcs.F
rdlooppd.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdlooppd.f
rdsrcemsmb.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdsrcemsmb.f
rdsrcidaar.o:
	$(FC) $(FFLAGS) -c -o $@ rdsrcidaar.f
rdsrcidamb.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ rdsrcidamb.f
rdsrcidapt.o:
	$(FC) $(FFLAGS) -c -o $@ rdsrcidapt.f
rdsrcntiar.o: $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ rdsrcntiar.f
rdsrcntimb.o: $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ rdsrcntimb.f
setnonhap.o: $(MODSOURC) $(MODLISTS) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ setnonhap.f
smkinven.o: $(MODSOURC) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ smkinven.f
srcmem.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ srcmem.f
wrepinven.o: $(MODSOURC) $(MODINFO) $(MODLISTS) $(MODAR2PT)
	$(FC) $(FFLAGS) -c -o $@ wrepinven.f
wrinvchr.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ wrinvchr.f
wrinvemis.o: $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODFILESET)
	$(FC) $(FFLAGS) -c -o $@ wrinvemis.F
wrinvpol.o: $(MODFILESET)
	$(FC) $(FFLAGS) -c -o $@ wrinvpol.f
wrpdemis.o: $(MODSOURC) $(MODINFO) $(MODDAYHR) 
	$(FC) $(FFLAGS) -c -o $@ wrpdemis.f
wrptref.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ wrptref.f

clean:
	/bin/rm -f *.o 
ifcclean:
	/bin/rm -f *.o  
	/bin/rm -f *.d  
	/bin/rm -f work.pc*
	/bin/rm -f $(MODALL)
	/bin/rm -f $(IOMODALL)
	/bin/rm -f $(FSMODALL)
