# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

#
# Include file for making modules and doing needed copies, etc.
#
MODBEIS3 = modbeis3.o
#
MODALL   = $(MODBEIS3) 
#
# Module dependencies...
# The copies from $EL_SRC to here are needed for the debug version, because
# for some reason when there are 6 or more modules in a directory used to 
# link the modules, the debugger will not be able to assign breakpoints
# at specific lines in the code
#

$(MODBEIS3): modbeis3.f $(INC)/EMPRVT3.EXT $(INC)/EMSTRG3.EXT
	$(FC) $(FFLAGS) -c -o $@ modbeis3.f
#	if ( test -f modbeis3.mod ) ; then mv modbeis3.mod $(MD_OUT) ; fi
#	if ( test -f modbeis3.M   ) ; then mv modbeis3.M   $(MD_OUT) ; fi
#	if ( test -f MODBEIS3.mod ) ; then mv MODBEIS3.mod $(MD_OUT) ; fi

