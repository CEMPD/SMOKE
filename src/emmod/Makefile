# Version @(#)$Id$
# Path    $Source$
# Date    $Date$ 

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<     

FILES =  modsourc.o modmobil.o modxref.o modtmprl.o modspro.o modcntrl.o \
         modlists.o modinfo.o moddayhr.o modsurg.o modemfac.o modstcy.o \
         modmerge.o modelev.o modmet.o modbiog.o modreprt.o modrepbn.o \
         modmbset.o modbeis3.o modar2pt.o
         

INCLUDES = $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT

include $(MD_SRC)/modmake.inc

clean:
	/bin/rm -f $(MD_OUT)/*.o $(MD_OUT)/*.mod $(MD_OUT)/*.M
