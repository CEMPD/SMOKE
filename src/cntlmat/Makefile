# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS)
LIBPATH       = $(SMK_BIN)/libsmoke.a  $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a  $(TOOLS_BIN)/libedsstools.debug.a
#

all: cntlmat

debug: cntlmat.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT 

#
CTL=cntlmat.o aloccmat.o alocpkts.o asgncntl.o errpkts.o \
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
include $(FS_ROOT)/modmake.inc

aloccmat.o: $(MODXREF) $(MODCNTRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ aloccmat.f
alocpkts.o: $(MODCNTRL) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ alocpkts.f
asgncntl.o: $(MODSOURC) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ asgncntl.f
cntlmat.o: $(MODSOURC) $(MODCNTRL) $(MODSPRO) $(MODINFO) $(MODFILESET)
	$(FC) $(FFLAGS) -c -o $@ cntlmat.f
fillcdat.o: $(MODCNTRL) $(CL_SRC)/CPKTDAT.EXT
	$(FC) $(FFLAGS) -c -o $@ fillcdat.f
fillcntl.o: $(MODXREF) $(MODLISTS) $(MODINFO) $(CL_SRC)/CPKTDAT.EXT
	$(FC) $(FFLAGS) -c -o $@ fillcntl.f
genmultc.o: $(MODSOURC) $(MODCNTRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ genmultc.f
genproj.o:  $(MODSOURC) $(MODCNTRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ genproj.f
genreact.o: $(MODSOURC) $(MODCNTRL) $(MODSPRO) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ genreact.f
opencmat.o: $(MODCNTRL) $(MODINFO) $(MODFILESET)
	$(FC) $(FFLAGS) -c -o $@ opencmat.f
openpmat.o: $(MODINFO) $(MODFILESET)
	$(FC) $(FFLAGS) -c -o $@ openpmat.f
openrmat.o: $(MODINFO) $(MODFILESET)
	$(FC) $(FFLAGS) -c -o $@ openrmat.f
pktloop.o: $(MODXREF) $(MODLISTS) $(MODINFO) $(CL_SRC)/CPKTDAT.EXT
	$(FC) $(FFLAGS) -c -o $@ pktloop.f
procpkts.o: $(MODCNTRL) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ procpkts.f
rdpacket.o: $(MODCNTRL) $(MODINFO) $(CL_SRC)/CPKTDAT.EXT
	$(FC) $(FFLAGS) -c -o $@ rdpacket.f
wrctmp.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ wrctmp.f
wcntlrep.o: $(MODSOURC) $(MODCNTRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ wcntlrep.f


#
clean:
	/bin/rm -f *.o
