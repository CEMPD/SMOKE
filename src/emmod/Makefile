# Version @(#)$Id$
# Path    $Source$
# Date    $Date$ 

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<     

FILES =  modsourc.o modxref.o modtpro.o modspro.o modcntrl.o modlists.o \
         modinfo.o

INCLUDES = $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT

include $(MD_SRC)/modmake.inc

clean:
	/bin/rm $(MD_OUT)/*.o $(MD_OUT)/*.mod
