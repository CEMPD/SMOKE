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
	mv $@ $(MD_OUT)
	if ( test -f MODSOURC.mod ) mv MODSOURC.mod $(MD_OUT)

modxref.o : $(INCLUDES)
	mv $@ $(MD_OUT)
	if ( test -f MODXREF.mod ) mv MODXREF.mod $(MD_OUT)

modspro.o : $(INCLUDES)
	mv $@ $(MD_OUT)
	if ( test -f MODSPRO.mod ) mv MODSPRO.mod $(MD_OUT)

modtpro.o : $(INCLUDES)
	mv $@ $(MD_OUT)
	if ( test -f MODTPRO.mod ) mv MODTPRO.mod $(MD_OUT)

modcntrl.o : $(INCLUDES)
	mv $@ $(MD_OUT)
	if ( test -f MODCNTRL.mod ) mv MODCNTRL.mod $(MD_OUT)

modlists.o : $(INCLUDES)
	mv $@ $(MD_OUT)
	if ( test -f MODLISTS.mod ) mv MODLISTS.mod $(MD_OUT)

modinfo.o : $(INCLUDES)
	mv $@ $(MD_OUT)
	if ( test -f MODINFO.mod ) mv MODINFO.mod $(MD_OUT)

clean:
	rm $MD_OUT/*.o $MD_OUT/*.mod
