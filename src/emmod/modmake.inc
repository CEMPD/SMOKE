# Version "@(#)$Id$ $Source$ $Date$ 
#
# Include file for making modules and doing needed copies, etc.
#
MD_OUT = $(MD_SRC)/$(EDSS_EXE)

MODSOURC = $(MD_OUT)/modsourc.o
MODMOBIL = $(MD_OUT)/modmobil.o
MODXREF  = $(MD_OUT)/modxref.o
MODTMPRL = $(MD_OUT)/modtmprl.o
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
MODBIOG  = $(MD_OUT)/modbiog.o
MODBEIS3 = $(MD_OUT)/modbeis3.o
MODREPRT = $(MD_OUT)/modreprt.o
MODREPBN = $(MD_OUT)/modrepbn.o

#
MODALL   = $(MODSOURC) $(MODMOBIL) $(MODXREF)  $(MODTMPRL) $(MODSPRO)  \
           $(MODCNTRL) $(MODLISTS) $(MODINFO)  $(MODDAYHR) \
           $(MODMERGE) $(MODSTCY)  $(MODSURG)  $(MODEMFAC) $(MODELEV)  \
           $(MODMET)   $(MODBIOG)  $(MODBEIS3) $(MODREPRT) $(MODREPBN)
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

$(MODMOBIL): $(MD_SRC)/modmobil.f  $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modmobil.f
	if ( test -f modmobil.mod ) ; then mv modmobil.mod $(MD_OUT) ; fi
	if ( test -f modmobil.M   ) ; then mv modmobil.M   $(MD_OUT) ; fi
	if ( test -f MODMOBIL.mod ) ; then mv MODMOBIL.mod $(MD_OUT) ; fi

$(MODXREF): $(MD_SRC)/modxref.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modxref.f
	if ( test -f modxref.mod ) ; then mv modxref.mod $(MD_OUT) ; fi
	if ( test -f modxref.M   ) ; then mv modxref.M   $(MD_OUT) ; fi
	if ( test -f MODXREF.mod ) ; then mv MODXREF.mod $(MD_OUT) ; fi

$(MODTMPRL):  $(MD_SRC)/modtmprl.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modtmprl.f
	if ( test -f modtmprl.mod ) ; then mv modtmprl.mod $(MD_OUT) ; fi
	if ( test -f modtmprl.M   ) ; then mv modtmprl.M   $(MD_OUT) ; fi
	if ( test -f MODTMPRL.mod ) ; then mv MODTMPRL.mod $(MD_OUT) ; fi

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

$(MODCNTRL): $(MD_SRC)/modcntrl.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modcntrl.f
	if ( test -f modcntrl.mod ) ; then mv modcntrl.mod $(MD_OUT) ; fi
	if ( test -f modcntrl.M   ) ; then mv modcntrl.M   $(MD_OUT) ; fi
	if ( test -f MODCNTRL.mod ) ; then mv MODCNTRL.mod $(MD_OUT) ; fi

$(MODSURG): $(MD_SRC)/modsurg.f
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modsurg.f
	if ( test -f modsurg.mod ) ; then mv modsurg.mod $(MD_OUT) ; fi
	if ( test -f modsurg.M   ) ; then mv modsurg.M   $(MD_OUT) ; fi
	if ( test -f MODSURG.mod ) ; then mv MODSURG.mod $(MD_OUT) ; fi

$(MODEMFAC): $(MD_SRC)/modemfac.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modemfac.f
	if ( test -f modemfac.mod ) ; then mv modemfac.mod $(MD_OUT) ; fi
	if ( test -f modemfac.M   ) ; then mv modemfac.M   $(MD_OUT) ; fi
	if ( test -f MODEMFAC.mod ) ; then mv MODEMFAC.mod $(MD_OUT) ; fi

$(MODSTCY): $(MD_SRC)/modstcy.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modstcy.f
	if ( test -f modstcy.mod ) ; then mv modstcy.mod $(MD_OUT) ; fi
	if ( test -f modstcy.M   ) ; then mv modstcy.M   $(MD_OUT) ; fi
	if ( test -f MODSTCY.mod ) ; then mv MODSTCY.mod $(MD_OUT) ; fi

$(MODMERGE): $(MD_SRC)/modmerge.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modmerge.f
	if ( test -f modmerge.mod ) ; then mv modmerge.mod $(MD_OUT) ; fi
	if ( test -f modmerge.M   ) ; then mv modmerge.M   $(MD_OUT) ; fi
	if ( test -f MODMERGE.mod ) ; then mv MODMERGE.mod $(MD_OUT) ; fi

$(MODELEV): $(MD_SRC)/modelev.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modelev.f
	if ( test -f modelev.mod ) ; then mv modelev.mod $(MD_OUT) ; fi
	if ( test -f modelev.M   ) ; then mv modelev.M   $(MD_OUT) ; fi
	if ( test -f MODELEV.mod ) ; then mv MODELEV.mod $(MD_OUT) ; fi

$(MODMET): $(MD_SRC)/modmet.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modmet.f
	if ( test -f modmet.mod ) ; then mv modmet.mod $(MD_OUT) ; fi
	if ( test -f modmet.M   ) ; then mv modmet.M   $(MD_OUT) ; fi
	if ( test -f MODMET.mod ) ; then mv MODMET.mod $(MD_OUT) ; fi

$(MODBIOG): $(MD_SRC)/modbiog.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modbiog.f
	if ( test -f modbiog.mod ) ; then mv modbiog.mod $(MD_OUT) ; fi
	if ( test -f modbiog.M   ) ; then mv modbiog.M   $(MD_OUT) ; fi
	if ( test -f MODBIOG.mod ) ; then mv MODBIOG.mod $(MD_OUT) ; fi

$(MODBEIS3): $(MD_SRC)/modbeis3.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modbeis3.f
	if ( test -f modbeis3.mod ) ; then mv modbeis3.mod $(MD_OUT) ; fi
	if ( test -f modbeis3.M   ) ; then mv modbeis3.M   $(MD_OUT) ; fi
	if ( test -f MODBEIS3.mod ) ; then mv MODBEIS3.mod $(MD_OUT) ; fi

$(MODREPRT): $(MD_SRC)/modreprt.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modreprt.f
	if ( test -f modreprt.mod ) ; then mv modreprt.mod $(MD_OUT) ; fi
	if ( test -f modreprt.M   ) ; then mv modreprt.M   $(MD_OUT) ; fi
	if ( test -f MODREPRT.mod ) ; then mv MODREPRT.mod $(MD_OUT) ; fi

$(MODREPBN): $(MD_SRC)/modrepbn.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ $(MD_SRC)/modrepbn.f
	if ( test -f modrepbn.mod ) ; then mv modrepbn.mod $(MD_OUT) ; fi
	if ( test -f modrepbn.M   ) ; then mv modrepbn.M   $(MD_OUT) ; fi
	if ( test -f MODREPBN.mod ) ; then mv MODREPBN.mod $(MD_OUT) ; fi
