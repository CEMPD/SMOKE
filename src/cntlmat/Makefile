# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)
LIBPATH       = $(SMK_BIN)/libsmoke.a  $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a  $(TOOLS_BIN)/libedsstools.debug.a
#

all: cntlmat

debug: cntlmat.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT 

#
CTL=cntlmat.o alocpkts.o asgncntl.o errpkts.o \
    fillcdat.o fillcntl.o  genmultc.o genproj.o genreact.o \
    openpmat.o opencmat.o openctmp.o openrmat.o pktloop.o \
    procpkts.o rdpacket.o wrctmp.o wrrmat.o wcntlrep.o
#
cntlmat: $(CTL) $(ALLI) $(LIBPATH)
	echo CNTLMAT:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(CTL) $(LIBS)
#
cntlmat.debug: $(CTL) $(ALLI) $(LIBPATH_DEBUG)
	echo CNTLMAT \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(CTL) $(LIBS_DEBUG)
#
# Module dependencies
#
include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

alocpkts.o: $(MODCNTRL) $(MODLISTS) $(MODINFO)
asgncntl.o: $(MODSOURC) $(MODXREF) $(MODINFO)
cntlmat.o: $(MODSOURC) $(MODCNTRL) $(MODSPRO) $(MODINFO) $(MODFILESET)
fillcdat.o: $(MODCNTRL) $(CL_SRC)/CPKTDAT.EXT
fillcntl.o: $(MODXREF) $(MODLISTS) $(MODINFO) $(CL_SRC)/CPKTDAT.EXT
genmultc.o: $(MODSOURC) $(MODCNTRL) $(MODINFO)
genproj.o:  $(MODSOURC) $(MODCNTRL) $(MODINFO)
genreact.o: $(MODSOURC) $(MODCNTRL) $(MODSPRO) $(MODINFO)
opencmat.o: $(MODCNTRL) $(MODINFO) $(MODFILESET)
openpmat.o: $(MODINFO) $(MODCNTRL) $(MODFILESET)
openrmat.o: $(MODINFO) $(MODFILESET)
pktloop.o: $(MODXREF) $(MODLISTS) $(MODINFO) $(CL_SRC)/CPKTDAT.EXT
procpkts.o: $(MODCNTRL) $(MODLISTS) $(MODINFO) $(MODXREF)
rdpacket.o: $(MODCNTRL) $(MODINFO) $(CL_SRC)/CPKTDAT.EXT
wrctmp.o: $(MODINFO)
wcntlrep.o: $(MODSOURC) $(MODCNTRL) $(MODINFO)

#
clean:
	/bin/rm -f *.o
