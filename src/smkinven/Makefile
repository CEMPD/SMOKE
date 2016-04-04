#.........................................................................
# Version "@(#)$Header$"
# Copyright (C) 2009 Baron Advanced Meteorological Systems, LLC.
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
#       ${SRCDIR}   is the source directory for this Makefile
#       ${INCDIR}   is the source directory for SMOKE INCLUDE-files
#       ${IODIR}    is the source directory for the I/O API library
#       ${OBJDIR}   is the current machine/compiler/flag-specific
#                   build-directory
#.........................................................................
#
#       ---------------     Definitions:   -------------------------

.SUFFIXES: .m4 .c .F .f .mod

include ../Makeinclude

SRCDIR = ${BASEDIR}/smkinven

FSRC =  \
adjustinv.f asgnar2pt.f asgnnhapx.f chklstfl.f fixstk.f formlist.f \
genmedsout.f genpdout.f gethdr.f getpdinfo.f initinfo.f \
opengrwout.f openinvin.f openinvout.f openpdout.f procar2pt.f procinven.f \
procinvsrcs.f rdcempd.f rdcemsum.f rddataff10ar.f rddataff10mb.f rddataff10pt.f \
rddatamedspt.f rddatantiar.f rddatantifr.f rddatantimb.f rddatantinp.f rddatantipt.f \
rdemspd.f rdff10pd.f rdgrdapi.f rdgrdncf.f rdinvdata.f rdinvsrcs.f \
rdlooppd.f rdmedsinfo.f rdmedspd.f rdorlfr.f rdsrcff10ar.f rdsrcff10mb.f \
rdsrcff10pt.f rdsrcmedspt.f rdsrcntiar.f rdsrcntifr.f rdsrcntimb.f rdsrcntinp.f \
rdsrcntipt.f setnonhap.f smkinven.f srcmem.f wrepinven.f wrinvchr.f \
wrinvemis.f wrinvpol.f wrpdemis.f wrptref.f grwinven.f

OBJ  =  $(FSRC:.f=.o)

EXE  = smkinven grwinven


#      ----------------------   TOP-LEVEL TARGETS:   ------------------

all: ${OBJDIR} ${EXE}

clean:
	cd ${OBJDIR}; /bin/rm ${OBJ} ${EXE}

rmexe:
	cd ${OBJDIR}; rm ${EXE}

relink:
	make -i rmexe ; make all

install:  all
	cd ${OBJDIR}; cp ${EXE} ${INSTDIR}


#      -----------------------   RULES:   -----------------------------

%.o : %.mod        #  Disable "gmake"s obnoxious implicit Modula-2 rule !!

.f.o:
	cd ${OBJDIR}; $(FC) $(FFLAGS) -c $(SRCDIR)/$<


#      -----------------------   DEPENDENCIES:------------------------

asgnar2pt.o : $(MODSOURC)
asgnar2pt.o : $(MODXREF)
asgnar2pt.o : $(MODINFO)
asgnnhapx.o : $(MODSOURC)
asgnnhapx.o : $(MODXREF)
asgnnhapx.o : $(MODINFO)
fixstk.o : $(MODSOURC)
formlist.o : $(MODLISTS)
formlist.o : $(MODINFO)
genmedsout.o : $(MODINFO)
genmedsout.o : $(MODDAYHR)
genpdout.o : $(MODSOURC)
genpdout.o : $(MODLISTS)
genpdout.o : $(MODINFO)
genpdout.o : $(MODDAYHR)
genpdout.o : $(MODSTCY)
gethdr.o : $(MODLISTS)
gethdr.o : $(MODSTCY)
gethdr.o : $(MODINFO)
getpdinfo.o : $(MODSOURC)
getpdinfo.o : $(MODINFO)
getpdinfo.o : $(MODDAYHR)
getpdinfo.o : $(MODLISTS)
grwinven.o : $(MODSOURC)
grwinven.o : $(MODINFO)
grwinven.o : $(MODFILESET)
initinfo.o : $(MODINFO)
opengrwout.o : $(MODINFO)
openinvin.o : $(MODINFO)
openinvin.o : $(MODLISTS)
openinvin.o : $(MODXREF)
openinvin.o : $(MODDAYHR)
openinvin.o : $(MODMOBIL)
openinvout.o : $(MODINFO)
openinvout.o : $(MODFILESET)
openpdout.o : $(MODINFO)
procar2pt.o : $(MODSOURC)
procar2pt.o : $(MODXREF)
procar2pt.o : $(MODINFO)
procar2pt.o : $(MODAR2PT)
procinven.o : $(MODSOURC)
procinven.o : $(MODLISTS)
procinven.o : $(MODINFO)
procinvsrcs.o : $(MODSOURC)
procinvsrcs.o : $(MODINFO)
rdcempd.o : $(MODSOURC)
rdcempd.o : $(MODLISTS)
rdcempd.o : $(MODINFO)
rdcempd.o : $(MODDAYHR)
rdcemsum.o : $(MODDAYHR)
rddataff10ar.o : $(MODINFO)
rddataff10ar.o : $(MODDAYHR)
rddataff10mb.o : $(MODINFO)
rddataff10mb.o : $(MODLISTS)
rddataff10pt.o : $(MODINFO)
rddataff10pt.o : $(MODDAYHR)
rddatamedspt.o : $(MODINFO)
rddatamedspt.o : $(MODSOURC)
rddatantiar.o : $(MODINFO)
rddatantifr.o : $(MODINFO)
rddatantimb.o : $(MODINFO)
rddatantinp.o : $(MODINFO)
rddatantipt.o : $(MODINFO)
rdemspd.o : $(MODSOURC)
rdemspd.o : $(MODLISTS)
rdemspd.o : $(MODINFO)
rdemspd.o : $(MODDAYHR)
rdff10pd.o : $(MODSOURC)
rdff10pd.o : $(MODLISTS)
rdff10pd.o : $(MODINFO)
rdff10pd.o : $(MODDAYHR)
rdff10pd.o : $(MODSTCY)
rdgrdapi.o : $(MODSOURC)
rdgrdapi.o : $(MODINFO)
rdgrdncf.o : $(MODSOURC)
rdgrdncf.o : $(MODINFO)
rdinvdata.o : $(MODSOURC)
rdinvdata.o : $(MODINFO)
rdinvdata.o : $(MODLISTS)
rdinvdata.o : $(MODMOBIL)
rdinvsrcs.o : $(MODSOURC)
rdinvsrcs.o : $(MODINFO)
rdinvsrcs.o : $(MODLISTS)
rdinvsrcs.o : $(MODMOBIL)
rdlooppd.o : $(MODINFO)
rdlooppd.o : $(MODLISTS)
rdmedsinfo.o : $(MODSOURC)
rdmedspd.o : $(MODSOURC)
rdmedspd.o : $(MODLISTS)
rdmedspd.o : $(MODINFO)
rdmedspd.o : $(MODDAYHR)
rdmedspd.o : $(MODSTCY)
rdorlfr.o : $(MODSOURC)
rdorlfr.o : $(MODLISTS)
rdorlfr.o : $(MODINFO)
rdorlfr.o : $(MODDAYHR)
rdsrcff10ar.o : $(MODLISTS)
rdsrcff10ar.o : $(MODDAYHR)
rdsrcff10mb.o : $(MODLISTS)
rdsrcff10pt.o : $(MODLISTS)
rdsrcff10pt.o : $(MODDAYHR)
rdsrcmedspt.o : $(MODSOURC)
rdsrcntiar.o : $(MODLISTS)
rdsrcntifr.o : $(MODINFO)
rdsrcntifr.o : $(MODLISTS)
rdsrcntimb.o : $(MODLISTS)
rdsrcntinp.o : $(MODLISTS)
rdsrcntipt.o : $(MODLISTS)
setnonhap.o : $(MODSOURC)
setnonhap.o : $(MODLISTS)
setnonhap.o : $(MODXREF)
setnonhap.o : $(MODINFO)
smkinven.o : $(MODSOURC)
smkinven.o : $(MODLISTS)
smkinven.o : $(MODINFO)
smkinven.o : $(MODDAYHR)
srcmem.o : $(MODSOURC)
wrepinven.o : $(MODSOURC)
wrepinven.o : $(MODINFO)
wrepinven.o : $(MODLISTS)
wrepinven.o : $(MODAR2PT)
wrinvchr.o : $(MODSOURC)
wrinvchr.o : $(MODINFO)
wrinvchr.o : $(MODLISTS)
wrinvemis.o : $(MODSOURC)
wrinvemis.o : $(MODLISTS)
wrinvemis.o : $(MODINFO)
wrinvemis.o : $(MODFILESET)
wrinvpol.o : $(MODFILESET)
wrpdemis.o : $(MODSOURC)
wrpdemis.o : $(MODLISTS)
wrpdemis.o : $(MODINFO)
wrpdemis.o : $(MODDAYHR)
wrpdemis.o : $(MODLISTS)
wrptref.o : $(MODINFO)

#  ---------------------------  ${EXE} Program builds:  -----------------

${OBJDIR}:
	mkdir -p ${OBJDIR}
smkinven: adjustinv.o asgnar2pt.o asgnnhapx.o chklstfl.o fixstk.o formlist.o \
          genmedsout.o genpdout.o gethdr.o getpdinfo.o initinfo.o \
          openinvin.o openinvout.o openpdout.o procar2pt.o procinven.o \
          procinvsrcs.o rdcempd.o rdcemsum.o rddataff10ar.o rddataff10mb.o rddataff10pt.o \
          rddatamedspt.o rddatantiar.o rddatantifr.o rddatantimb.o rddatantinp.o rddatantipt.o \
          rdemspd.o rdff10pd.o rdgrdapi.o rdgrdncf.o rdinvdata.o rdinvsrcs.o \
          rdlooppd.o rdmedsinfo.o rdmedspd.o rdorlfr.o rdsrcff10ar.o rdsrcff10mb.o \
          rdsrcff10pt.o rdsrcmedspt.o rdsrcntiar.o rdsrcntifr.o rdsrcntimb.o rdsrcntinp.o \
          rdsrcntipt.o setnonhap.o smkinven.o srcmem.o wrepinven.o wrinvchr.o \
          wrinvemis.o wrinvpol.o wrpdemis.o wrptref.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

grwinven: opengrwout.o srcmem.o wrinvchr.o wrinvpol.o grwinven.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

