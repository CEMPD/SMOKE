# Version "@(#)$Id$ $Source$ $Date$ 
#
# Include file for making modules and doing needed copies, etc.
#
MD_OUT = $(MD_SRC)/$(SMOKE_EXE)

MODSOURC = $(MD_OUT)/modsourc.o
MODXREF  = $(MD_OUT)/modxref.o
MODTPRO  = $(MD_OUT)/modtpro.o
MODSPRO  = $(MD_OUT)/modspro.o
MODLISTS = $(MD_OUT)/modlists.o
MODINFO  = $(MD_OUT)/modinfo.o
MODDAYHR = $(MD_OUT)/moddayhr.o
MODCNTRL = $(MD_OUT)/modcntrl.o
MODSURG  = $(MD_OUT)/modsurg.o
MODEMFAC = $(MD_OUT)/modemfac.o
MODSTCY  = $(MD_OUT)/modstcy.o
MODMERGE = $(MD_OUT)/modmerge.o
MODELEV  = $(MD_OUT)/modelev.o
MODMET   = $(MD_OUT)/modmet.o
#
MODALL   = $(MODSOURC) $(MODXREF) $(MODTPRO)  $(MODSPRO)  $(MODCNTRL) \
           $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODMERGE) $(MODSTCY)  \
           $(MODSURG)  $(MODEMFAC) $(MODELEV) $(MODMET)
#
# Module dependencies...
# The copies from $EL_SRC to here are needed for the debug version, because
# for some reason when there are 6 or more modules in a directory used to 
# link the modules, the debugger will not be able to assign breakpoints
# at specific lines in the code
#
$(MODSOURC): $(MD_SRC)/modsourc.f  $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modsourc.f
	if ( test -f modsourc.mod ) ; then mv modsourc.mod $(MD_OUT) ; fi
	if ( test -f modsourc.M   ) ; then mv modsourc.M   $(MD_OUT) ; fi
	if ( test -f MODSOURC.mod ) ; then mv MODSOURC.mod $(MD_OUT) ; fi

$(MODXREF): $(MD_SRC)/modxref.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modxref.f
	if ( test -f modxref.mod ) ; then mv modxref.mod $(MD_OUT) ; fi
	if ( test -f modxref.M   ) ; then mv modxref.M   $(MD_OUT) ; fi
	if ( test -f MODXREF.mod ) ; then mv MODXREF.mod $(MD_OUT) ; fi

$(MODTPRO):  $(MD_SRC)/modtpro.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modtpro.f
	if ( test -f modtpro.mod ) ; then mv modtpro.mod $(MD_OUT) ; fi
	if ( test -f modtpro.M   ) ; then mv modtpro.M   $(MD_OUT) ; fi
	if ( test -f MODTPRO.mod ) ; then mv MODTPRO.mod $(MD_OUT) ; fi

$(MODSPRO): $(MD_SRC)/modspro.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modspro.f
	if ( test -f modspro.mod ) ; then mv modspro.mod $(MD_OUT) ; fi
	if ( test -f modspro.M   ) ; then mv modspro.M   $(MD_OUT) ; fi
	if ( test -f MODSPRO.mod ) ; then mv MODSPRO.mod $(MD_OUT) ; fi

$(MODLISTS): $(MD_SRC)/modlists.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modlists.f
	if ( test -f modlists.mod ) ; then mv modlists.mod $(MD_OUT) ; fi
	if ( test -f modlists.M   ) ; then mv modlists.M   $(MD_OUT) ; fi
	if ( test -f MODLISTS.mod ) ; then mv MODLISTS.mod $(MD_OUT) ; fi

$(MODINFO): $(MD_SRC)/modinfo.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modinfo.f
	if ( test -f modinfo.mod ) ; then mv modinfo.mod $(MD_OUT) ; fi
	if ( test -f modinfo.M   ) ; then mv modinfo.M   $(MD_OUT) ; fi
	if ( test -f MODINFO.mod ) ; then mv MODINFO.mod $(MD_OUT) ; fi

$(MODDAYHR): $(MD_SRC)/moddayhr.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/moddayhr.f
	if ( test -f moddayhr.mod ) ; then mv moddayhr.mod $(MD_OUT) ; fi
	if ( test -f moddayhr.M   ) ; then mv moddayhr.M   $(MD_OUT) ; fi
	if ( test -f MODDAYHR.mod ) ; then mv MODDAYHR.mod $(MD_OUT) ; fi

$(MODCNTRL): $(MD_SRC)/modcntrl.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modcntrl.f
	if ( test -f modcntrl.mod ) ; then mv modcntrl.mod $(MD_OUT) ; fi
	if ( test -f modcntrl.M   ) ; then mv modcntrl.M   $(MD_OUT) ; fi
	if ( test -f MODCNTRL.mod ) ; then mv MODCNTRL.mod $(MD_OUT) ; fi

$(MODSURG): $(MD_SRC)/modsurg.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modsurg.f
	if ( test -f modsurg.mod ) ; then mv modsurg.mod $(MD_OUT) ; fi
	if ( test -f modsurg.M   ) ; then mv modsurg.M   $(MD_OUT) ; fi
	if ( test -f MODSURG.mod ) ; then mv MODSURG.mod $(MD_OUT) ; fi

$(MODEMFAC): $(MD_SRC)/modemfac.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modemfac.f
	if ( test -f modemfac.mod ) ; then mv modemfac.mod $(MD_OUT) ; fi
	if ( test -f modemfac.M   ) ; then mv modemfac.M   $(MD_OUT) ; fi
	if ( test -f MODEMFAC.mod ) ; then mv MODEMFAC.mod $(MD_OUT) ; fi

$(MODSTCY): $(MD_SRC)/modstcy.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modstcy.f
	if ( test -f modstcy.mod ) ; then mv modstcy.mod $(MD_OUT) ; fi
	if ( test -f modstcy.M   ) ; then mv modstcy.M   $(MD_OUT) ; fi
	if ( test -f MODSTCY.mod ) ; then mv MODSTCY.mod $(MD_OUT) ; fi

$(MODMERGE): $(MD_SRC)/modmerge.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modmerge.f
	if ( test -f modmerge.mod ) ; then mv modmerge.mod $(MD_OUT) ; fi
	if ( test -f modmerge.M   ) ; then mv modmerge.M   $(MD_OUT) ; fi
	if ( test -f MODMERGE.mod ) ; then mv MODMERGE.mod $(MD_OUT) ; fi

$(MODELEV): $(MD_SRC)/modelev.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modelev.f
	if ( test -f modelev.mod ) ; then mv modelev.mod $(MD_OUT) ; fi
	if ( test -f modelev.M   ) ; then mv modelev.M   $(MD_OUT) ; fi
	if ( test -f MODELEV.mod ) ; then mv MODELEV.mod $(MD_OUT) ; fi

$(MODMET): $(MD_SRC)/modmet.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modmet.f
	if ( test -f modmet.mod ) ; then mv modmet.mod $(MD_OUT) ; fi
	if ( test -f modmet.M   ) ; then mv modmet.M   $(MD_OUT) ; fi
	if ( test -f MODMET.mod ) ; then mv MODMET.mod $(MD_OUT) ; fi

