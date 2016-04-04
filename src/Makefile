#.........................................................................
# Version "@(#)$Header$"
# Copyright (C) 2009 Baron Advanced Meteorological Systems, LLC., and
# (C) 2016 UNC Institute for the Environment
# Distributed under the GNU GENERAL PUBLIC LICENSE version 2
# See file "GPL.txt" for conditions of use.
#.........................................................................
#  Environment Variables:
#       BIN     machine/OS/compiler/mode type. Shows up as suffix
#               for "Makeinclude.${BIN}" to determine compilation
#               flags, and in ${OBJDIR} and $(INSTALL) to determine
#               binary directories
#       INSTALL installation-directory root, used for "make install":
#               SMOKE executables will be installed in $(INSTALL)/${BIN}
#.........................................................................
#  Directories:
#       ${BASEDIR}  is the root directory for the SMOKE source and
#                   the (machine/compiler/flag-specific) binary
#                   object/library/executable directories.

include Makeinclude

all:
	if [ ! -d ${OBJDIR} ]; then mkdir ${OBJDIR}; fi;
	cd ${BASEDIR}/filesetapi; make all
	cd ${BASEDIR}/emmod     ; make all
	cd ${BASEDIR}/lib       ; make all
	cd ${BASEDIR}/emqa      ; make all
	cd ${BASEDIR}/emutil    ; make all
	cd ${BASEDIR}/biog      ; make all
	cd ${BASEDIR}/cntlmat   ; make all
	cd ${BASEDIR}/grdmat    ; make all
	cd ${BASEDIR}/spcmat    ; make all
	cd ${BASEDIR}/point     ; make all
	cd ${BASEDIR}/movesmrg  ; make all
	cd ${BASEDIR}/smkinven  ; make all
	cd ${BASEDIR}/smkmerge  ; make all
	cd ${BASEDIR}/temporal  ; make all
	echo -e "SMOKE:  completed 'make all'"

clean:
	cd ${BASEDIR}/filesetapi; make -i clean
	cd ${BASEDIR}/emmod     ; make -i clean
	cd ${BASEDIR}/lib       ; make -i clean
	cd ${BASEDIR}/emqa      ; make -i clean
	cd ${BASEDIR}/emutil    ; make -i clean
	cd ${BASEDIR}/biog      ; make -i clean
	cd ${BASEDIR}/cntlmat   ; make -i clean
	cd ${BASEDIR}/grdmat    ; make -i clean
	cd ${BASEDIR}/spcmat    ; make -i clean
	cd ${BASEDIR}/point     ; make -i clean
	cd ${BASEDIR}/movesmrg  ; make -i clean
	cd ${BASEDIR}/smkinven  ; make -i clean
	cd ${BASEDIR}/smkmerge  ; make -i clean
	cd ${BASEDIR}/temporal  ; make -i clean
	echo -e "SMOKE:  completed 'make clean'"

install:
	echo -e "\nSMOKE:  making install in ${INSTDIR}"
	if [ ! -d ${INSTDIR} ] ; then mkdir ${INSTDIR} ; fi ;
	cd ${BASEDIR}/biog      ; make install
	cd ${BASEDIR}/cntlmat   ; make install
	cd ${BASEDIR}/emqa      ; make install
	cd ${BASEDIR}/emutil    ; make install
	cd ${BASEDIR}/grdmat    ; make install
	cd ${BASEDIR}/spcmat    ; make install
	cd ${BASEDIR}/emutil    ; make install
	cd ${BASEDIR}/point     ; make install
	cd ${BASEDIR}/movesmrg  ; make install
	cd ${BASEDIR}/smkinven  ; make install
	cd ${BASEDIR}/smkmerge  ; make install
	cd ${BASEDIR}/temporal  ; make install
	echo -e "SMOKE:  completed 'make install' in ${INSTDIR}"

relink:
	cd ${BASEDIR}/biog      ; make relink
	cd ${BASEDIR}/cntlmat   ; make relink
	cd ${BASEDIR}/emqa      ; make relink
	cd ${BASEDIR}/emutil    ; make relink
	cd ${BASEDIR}/grdmat    ; make relink
	cd ${BASEDIR}/spcmat    ; make relink
	cd ${BASEDIR}/emutil    ; make relink
	cd ${BASEDIR}/point     ; make relink
	cd ${BASEDIR}/movesmrg  ; make relink
	cd ${BASEDIR}/smkinven  ; make relink
	cd ${BASEDIR}/smkmerge  ; make relink
	cd ${BASEDIR}/temporal  ; make relink
	echo -e "SMOKE:  completed 'make relink'"

