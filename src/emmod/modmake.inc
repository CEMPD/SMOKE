# Version "@(#)$Id$ $Source$ $Date$ 
#
# Include file for making modules and doing needed copies, etc.
#
MD_OUT = $(MD_SRC)/$(EDSS_EXE)

MODSOURC = $(MD_OUT)/modsourc.o
MODXREF  = $(MD_OUT)/modxref.o
MODTPRO  = $(MD_OUT)/modtpro.o
MODSPRO  = $(MD_OUT)/modspro.o
MODLISTS = $(MD_OUT)/modlists.o
MODINFO  = $(MD_OUT)/modinfo.o
MODSURG  = $(MD_OUT)/modsurg.o
#
# Module dependencies...
# The copies from $EL_SRC to here are needed for the debug version, because
# for some reason when there are 6 or more modules in a directory used to 
# link the modules, the debugger will not be able to assign breakpoints
# at specific lines in the code
#
$(MODSOURC): $(MD_SRC)/modsourc.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modsourc.f
	if ( test -f modsourc.mod ) ; then mv modsourc.mod $(MD_OUT) ; fi
	if ( test -f MODSOURC.mod ) ; then mv MODSOURC.mod $(MD_OUT) ; fi

$(MODXREF): $(MD_SRC)/modxref.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modxref.f
	if ( test -f modxref.mod ) ; then mv modxref.mod $(MD_OUT) ; fi
	if ( test -f MODXREF.mod ) ; then mv MODXREF.mod $(MD_OUT) ; fi

$(MODTPRO):  $(MD_SRC)/modtpro.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modtpro.f
	if ( test -f modtpro.mod ) ; then mv modtpro.mod $(MD_OUT) ; fi
	if ( test -f MODTPRO.mod ) ; then mv MODTPRO.mod $(MD_OUT) ; fi

$(MODSPRO): $(MD_SRC)/modspro.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modspro.f
	if ( test -f modspro.mod ) ; then mv modspro.mod $(MD_OUT) ; fi
	if ( test -f MODSPRO.mod ) ; then mv MODSPRO.mod $(MD_OUT) ; fi

$(MODLISTS): $(MD_SRC)/modlists.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modlists.f
	if ( test -f modlists.mod ) ; then mv modlists.mod $(MD_OUT) ; fi
	if ( test -f MODLISTS.mod ) ; then mv MODLISTS.mod $(MD_OUT) ; fi

$(MODINFO): $(MD_SRC)/modinfo.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modinfo.f
	if ( test -f modinfo.mod ) ; then mv modinfo.mod $(MD_OUT) ; fi
	if ( test -f MODINFO.mod ) ; then mv MODINFO.mod $(MD_OUT) ; fi

$(MODSURG): $(MD_SRC)/modsurg.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modsurg.f
	if ( test -f modsurg.mod ) ; then mv modsurg.mod $(MD_OUT) ; fi
	if ( test -f MODSURG.mod ) ; then mv MODSURG.mod $(MD_OUT) ; fi

