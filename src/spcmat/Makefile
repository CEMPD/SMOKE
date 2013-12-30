# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)
LIBPATH       = $(SMK_BIN)/libsmoke.a  $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a
#

all: spcmat

debug: spcmat.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
SPC= spcmat.o opensmat.o asgnspro.o asgntag.o rdtag.o chknonhap.o rdcombo.o
     
#
spcmat: $(SPC) $(ALLI) $(LIBPATH)
	echo SPCMAT:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(SPC) $(LIBS)
#
#
spcmat.debug: $(SPC) $(ALLI) $(LIBPATH_DEBUG)
	echo SPCMAT \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(SPC) $(LIBS_DEBUG)

include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

chknonhap.o: $(MODLISTS) $(MODSPRO)
asgnspro.o: $(MODSOURC) $(MODXREF) $(MODSPRO) $(MODINFO)
asgntag.o: $(MODTAG) $(MODSOURC) $(MODINFO)
opensmat.o: $(MODINFO) $(MODFILESET) $(MODTAG)
rdcombo.o: $(MODLISTS) $(MODSPRO)
rdtag.o: $(MODLISTS) $(MODXREF) $(MODSPRO) $(MODTAG) $(MODINFO)
spcmat.o: $(MODSOURC) $(MODSPRO) $(MODINFO) $(MODEMFAC) $(MODFILESET) $(MODLISTS) $(MODTAG)
#
clean:
	/bin/rm -f *.o
