FC = f90
FFLAGS   = -g -I$(INC)
DBGFLAGS = -g

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<     

FILES =  modsourc.o modxref.o modtpro.o modspro.o modcntrl.o modlists.o \
         modinfo.o

INCLUDES = $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT

debug: $(FILES)

modsourc.o: $(INCLUDES)
	mv $@ $(EDSS_EXE)
	if ( test -f MODSOURC.mod ) mv MODSOURC.mod $(EDSS_EXE)

modxref.o : $(INCLUDES)
	mv $@ $(EDSS_EXE)
	if ( test -f MODXREF.mod ) mv MODXREF.mod $(EDSS_EXE)

modspro.o : $(INCLUDES)
	mv $@ $(EDSS_EXE)
	if ( test -f MODSPRO.mod ) mv MODSPRO.mod $(EDSS_EXE)

modtpro.o : $(INCLUDES)
	mv $@ $(EDSS_EXE)
	if ( test -f MODTPRO.mod ) mv MODTPRO.mod $(EDSS_EXE)

modcntrl.o : $(INCLUDES)
	mv $@ $(EDSS_EXE)
	if ( test -f MODCNTRL.mod ) mv MODCNTRL.mod $(EDSS_EXE)

modlists.o : $(INCLUDES)
	mv $@ $(EDSS_EXE)
	if ( test -f MODLISTS.mod ) mv MODLISTS.mod $(EDSS_EXE)

modinfo.o : $(INCLUDES)
	mv $@ $(EDSS_EXE)
	if ( test -f MODINFO.mod ) mv MODINFO.mod $(EDSS_EXE)

clean:
	rm $SMK_EXE/*.o $SMK_EXE/*.mod
