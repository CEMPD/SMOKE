#.........................................................................
# Version "$Id$"
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
#.........................................................................
#
#       ---------------     Definitions:   -------------------------

.SUFFIXES: .m4 .c .F .f .f90 .F90 .mod

include Makeinclude

# find source, object, and module files here:

VPATH = ${OBJDIR}:\
${BASEDIR}/emmod:${BASEDIR}/lib:${BASEDIR}/filesetapi:${BASEDIR}/biog:\
${BASEDIR}/cntlmat:${BASEDIR}/emqa:${BASEDIR}/emutil:${BASEDIR}/grdmat:\
${BASEDIR}/movesmrg:${BASEDIR}/point:${BASEDIR}/smkinven:${BASEDIR}/smkmerge:\
${BASEDIR}/spcmat:${BASEDIR}/temporal

#  emmod:

MODSRC = \
 modar2pt.f   modbeis3.f    modbiog.f         \
 modcntrl.f   moddayhr.f    modelev.f         \
 modemfac.f   modgrid.f     modinfo.f         \
 modlists.f   modmbset.f    modmerge.f        \
 modmet.f     modmobil.f    modmvsmrg.f       \
 modrepbn.f   modreprt.f    modsourc.f        \
 modspro.f    modstcy.f     modsurg.f         \
 modtag.f     modtmprl.f    modxref.f         \
 modgrdlib.f 

MOD  = $(MODSRC:.f=.mod)

MODOBJ = $(MODSRC:.f=.o)

#  lib:

LIBSRC = \
 alocatbl.f     alocchrt.f      alocctbl.f          \
 alocetbl.f     alocgtbl.f      alocmtbl.f          \
 alocptbl.f     alocstbl.f      alocttbl.f          \
 applreac.f     applumat.f      bldcsrc.f           \
 bldenams.f     blkorcmt.f      calcrelhum.f        \
 checkmem.f     chkcpvar.f      chkemepi.f          \
 chkexpscc.f    chkexpsic.f     chkgrid.f           \
 chkint.f       chkisiz.f       chkmetem.f          \
 chkptdef.f     chkreal.f       chksrcno.f          \
                                cvtrdtype.f         \
 cvtvehtype.f   dscm3grd.f      dscm3lay.f          \
 dscsprof.f     efsetup.f       evalcrit.f          \
 fillatbl.f     fillchrt.f      fillctbl.f          \
 filletbl.f     fillgtbl.f      fillmtbl.f          \
 fillptbl.f     fillstbl.f      fillttbl.f          \
 find1first.f   flterr.f        fltrneg.f           \
 fltrxref.f     fmtcsrc.f       genptcel.f          \
 genptvcel.f    genuslst.f      getbasyr.f          \
 getcfdsc.f     getctgry.f      getdysav.f          \
 getfline.f     getformt.f      getidasz.f          \
 getifdsc.f     getiname.f      getinvyr.f          \
 getm3epi.f     getnlist.f      getntisz.f          \
 getrfdsc.f     getsinfo.f      gettzone.f          \
 grd2cnty.f     grdfips.f       hdrmiss3.f          \
 ingrid.f       initem.f        ioapi_grd_size.f    \
 rdsprofdsc.f   multunit.f      nextime.f           \
 normtpro.f     openphys.f      padnzero.f          \
 padzero.f      parscsrc.f      parsline.f          \
 pdsetup.f      polmesg.f       prclinrc.f          \
 procspro.f     progdesc.f      rd3mask.f           \
 rdar2pt.f      rdascc.f        rdchrscc.f          \
 rdcodnam.f     rddates.f       rdeproc.f           \
 rdgeocodes.f   rdgmat.f        rdgref.f            \
 rdhdays.f      rdinvchr.f      rdinvmap.f          \
 rdinvpol.f     rdlines.f       rdmactdsc.f         \
 rdmapmask.f    rdmappol.f      rdnaicsdsc.f        \
 rdorsdsc.f     rdpelv.f        rdrmat.f            \
 rdsccdsc.f     rdsccmap.f      rdsconv.f           \
 rdsetmask.f    rdsicdsc.f      rdsmat.f            \
 rdspdprof.f    rdspdref.f      rdsprof.f           \
 rdsrcgrps.f    rdsref.f        rdsrg.f             \
 rdsrgdesc.f    rdsrghdr.f      rdstcy.f            \
 rdtref.f       rdtzone.f       rdumat.f            \
 rdxclude.f     readwr3.f       rmcommnt.f          \
 setscctype.f   setsrcdy.f      srcgrpcnt.f         \
 tagtable.f                     unitfac.f           \
 unitmatch.f    useexpgeo.f     verchar.f           \
 wrchrscc.f     wrdaymsg.f      wrorlout.f          \
 wrsrcgrps.f    xfractbl.f      xreftbl.f

LIBOBJ = $(LIBSRC:.f=.o)

#  filesetapi:  both ".f" and ".F"
#  special case because of both .f and .F in ${FILSRC}

FILSRC = \
 appendname.f        chkfileset.f         \
 chksetdesc.f        cleanup.f            \
 closeset.f          createset.f          \
 descset.f           modfileset.f         \
 openset.f           promptset.f          \
 readset.F           writeset.F

FILOBJ = \
 appendname.o        chkfileset.o         \
 chksetdesc.o        cleanup.o            \
 closeset.o          createset.o          \
 descset.o           modfileset.o         \
 openset.o           promptset.o          \
 readset.o           writeset.o

FIL  = $(FILSRC:.f=.mod)

#  biog:

BIOSRC = \
 czangle.f     getpar.f      getparb.f           \
 hrbeis3.f     hrbeis_361.f  hrno.f              \
 normbeis3.f   normbeis314.f normbeis361.f       \
 prebmet.f     rdb3fac.f     rdb4fac.f           \
 rdb4fac_csv.f rdbpro.f      tmpbeis3.f          \
 tmpbeis314.f  tmpbeis361.f
 
BIOOBJ =  $(BIOSRC:.f=.o)

#  cntlmat:

CTLSRC = \
 alocpkts.f     asgncntl.f      cntlmat.f       \
 errpkts.f      fillcdat.f      fillcntl.f      \
 genmultc.f     genproj.f       genreact.f      \
 opencmat.f     openctmp.f      openpmat.f      \
 openrmat.f     pktloop.f       procpkts.f      \
 rdpacket.f     wcntlrep.f      wrctmp.f        \
 wrrmat.f
 
CTLOBJ =  $(CTLSRC:.f=.o)

#  emqa (smkreport):

RPTSRC = \
 asgnbins.f    bldrepidx.f    genrprt.f  \
 openrepin.f   openrepout.f   qarepin.f  \
 rdgrps.f      rdrepin.f      rdrprts.f  \
 repmrggrd.f   repunits.f     scanrepc.f \
 selectsrc.f   smkreport.f    wrrephdr.f \
 wrrepout.f    rdssup.f
 
RPTOBJ =  $(RPTSRC:.f=.o)

#  emutil:

UTISRC = \
 aggwndw.f   beld3to2.f       bluesky2inv.f        \
 cemscan.f   extractida.f     gcntl4carb.f         \
 gentpro.f   geofac.f         invsplit.f           \
 layalloc.f  metcombine.f     metscan.f            \
 pktreduc.f  smk2emis.f       surgtool.f           \
 uam2ncf.f   inlineto2d.f     saregroup.f
 
UTIOBJ =  $(UTISRC:.f=.o)

#  grdmat:

GRDSRC = \
 asgnsurg.f  genagmat.f       genggmat.f    \
 genlgmat.f  genmgmat.f       genpgmat.f    \
 grdmat.f    opengmat.f       rdsrg4grd.f   \
 setfrac.f   sizgmat.f
 
GRDOBJ =  $(GRDSRC:.f=.o)

#  movesmrg:

MOVSRC = \
 bldprocidx.f     bldsrccell.f      hourmet.f     \
 mbldmrgidx.f     met4moves.f       mgetmrgev.f   \
 minitstcy.f      mmrgonams.f       mmrgunits.f   \
 mmrgvnams.f      mopenmrgin.f      mopenmrgout.f \
 movesmrg.f       mwrmrgrep.f       rdcfpro.f     \
 rdfmref.f        rdmrclist.f       rdmxref.f     \
 rdrpdemfacs.f    rdrphemfacs.f     rdrppemfacs.f \
 rdrpvemfacs.f    rdspdpro.f        setrefcnty.f  \
 wrtemprof.f      rdspdist.f
 
MOVOBJ =  $(MOVSRC:.f=.o)

#  point:

PNTSRC = \
 asgngrps.f       elevpoint.f       fire_plmris.f     \
 fire_postplm.f   fire_preplm.f     laypoint.f        \
 mxgrpemis.f      openeout.f        openlayout.f      \
 plmris.f         plsprd.f          plumesti.f        \
 postplm.f        preplm.f          rpelvcfg.f        \
 wpingstk.f
 
PNTOBJ =  $(PNTSRC:.f=.o)

#  smkinven:

INVSRC = \
 adjustinv.f    asgnar2pt.f    asgnnhapx.f    chklstfl.f     \
 fixstk.f       formlist.f     genmedsout.f   \
 genpdout.f     gethdr.f       getpdinfo.f    \
 grwinven.f     initinfo.f     opengrwout.f   \
 openinvin.f    openinvout.f   openpdout.f    \
 procar2pt.f    procinven.f    procinvsrcs.f  \
 rdcempd.f      rdcemsum.f     rddataff10ar.f \
 rddataff10mb.f rddataff10pt.f rddatamedspt.f \
 rddatantiar.f  rddatantifr.f  rddatantimb.f  \
 rddatantinp.f  rddatantipt.f  rdemspd.f      \
 rdff10pd.f     rdgrdapi.f     rdgrdncf.f     \
 rdinvdata.f    rdinvsrcs.f    rdlooppd.f     \
 rdmedsinfo.f   rdmedspd.f     rdorlfr.f      \
 rdsrcff10ar.f  rdsrcff10mb.f  rdsrcff10pt.f  \
 rdsrcmedspt.f  rdsrcntiar.f   rdsrcntifr.f   \
 rdsrcntimb.f   rdsrcntinp.f   rdsrcntipt.f   \
 setnonhap.f    smkinven.f     srcmem.f       \
 wrepinven.f    wrinvchr.f     wrinvemis.f    \
 wrinvpol.f     wrpdemis.f     wrptref.f

INVOBJ =  $(INVSRC:.f=.o)

#  smkmerge:

MRGSRC = \
 allocmrg.f    bldmrgidx.f    getmrgev.f     \
 initstcy.f    mrgasciielev.f mrgbio.f       \
 mrgelev.f     mrggrid.f      mrgmult.f      \
 mrgonams.f    mrgunits.f     mrgvnams.f     \
 openmrgin.f   openmrgout.f   rdmrginv.f     \
 setoutdate.f  smkmerge.f     wmrgelev.f     \
 wmrgemis.f    wrmrgrep.f     mrgpt.f

MRGOBJ =  $(MRGSRC:.f=.o)

#  spcmat:

SPCSRC = \
 asgnspro.f       asgntag.f        chknonhap.f      \
 opensmat.f       rdcombo.f        rdtag.f          \
 spcmat.f

SPCOBJ =  $(SPCSRC:.f=.o)

#  temporal:

TMPSRC = \
 asgntpro.f    genhemis.f     mktmat.f       \
 opentmp.f     opentmpin.f    proctpro.f     \
 setdaylt.f    temporal.f     tmnamunt.f     \
 updtmat.f     wrtsup.f

TMPOBJ =  $(TMPSRC:.f=.o)

OBJ  =  $(MODOBJ) $(FILOBJ) $(LIBOBJ) $(BIOOBJ) $(CTLOBJ) $(RPTOBJ) $(UTIOBJ) \
        $(GRDOBJ) $(MOVOBJ) $(PNTOBJ) $(INVOBJ) $(MRGOBJ) $(SPCOBJ) $(TMPOBJ)

MLIB  = libemmod.a

SLIB  = libsmoke.a

FLIB  = libfileset.a

UTIL  = aggwndw     beld3to2    bluesky2inv cemscan     extractida  \
        geofac      invsplit    metcombine  metscan     pktreduc    \
        smk2emis    surgtool    uam2ncf     gentpro     gcntl4carb  \
        layalloc    inlineto2d  saregroup

EXE   = normbeis3   tmpbeis3    cntlmat     smkreport   grdmat      \
        movesmrg    met4moves   laypoint    elevpoint   smkinven    \
        grwinven    smkmerge    mrggrid     mrgelev     spcmat      \
        temporal    mrgpt                                           \
        ${UTIL}


######################################################################
#                           High-level targets:

all:    lib exe

lib:    ${MLIB} ${SLIB} ${FLIB}

exe:    ${EXE}

dir:    ${OBJDIR}

emmod:  ${MOD} ${MLIB}

smklib: ${SLIB}

biog:   normbeis3   tmpbeis3

cntl:   cntlmat

report: smkreport

grd  :  grdmat

moves:  movesmrg    met4moves

point:  laypoint    elevpoint

inven:  smkinven    grwinven

merge:  smkmerge    mrggrid     mrgelev      mrgpt

spc  :  spcmat

tmprl:  temporal

util :  ${UTIL}

clean:
	cd ${OBJDIR}; rm -f ${OBJ} ${EXE} ${MOD} ${FIL} ${MLIB} ${SLIB} ${FLIB}

rmexe:
	cd ${OBJDIR}; rm ${EXE}

relink:
	make BIN=${BIN} -i rmexe ; make BIN=${BIN} all

install:  ${INSTDIR} ${EXE}
	echo "Installing SMOKE in ${INSTDIR}"
	cd ${OBJDIR};  cp ${EXE} ${INSTDIR}


#      -----------------------   RULES:   -------------------------

%.o : %.mod        #  Disable "gmake"s obnoxious implicit Modula-2 rule !!

.f.o:
	cd ${OBJDIR};  ${FC} ${FFLAGS} -c $<

.F.o:
	cd ${OBJDIR};  ${FC} ${FFLAGS} -c $<

.f.mod:
	cd ${OBJDIR};  ${FC} ${FFLAGS} -c $<


#  ---------------------------  Dependencies:  --------------------------
# libraries must be up-to-date to build executables

${EXE} :  ${MLIB} ${SLIB} ${FLIB}       # libraries must be up-to-date to build executables

# modules must be up-to-date to build ".o"-files

alocatbl.o    : $(MODXREF)
alocchrt.o    : $(MODXREF)
alocctbl.o    : $(MODXREF)
alocetbl.o    : $(MODXREF)
alocgtbl.o    : $(MODXREF)
alocmtbl.o    : $(MODXREF)
alocptbl.o    : $(MODXREF)
alocstbl.o    : $(MODXREF)
alocttbl.o    : $(MODXREF)
applumat.o    : $(MODGRID)
bldenams.o    : $(MODLISTS)
chkgrid.o     : $(MODGRID)
chkmetem.o    : $(MODGRDLIB)
chkptdef.o    : $(MODINFO)
dscsprof.o    : $(MODINFO)
dscsprof.o    : $(MODSPRO)
efsetup.o     : $(MODEMFAC)
fillatbl.o    : $(MODXREF)
fillchrt.o    : $(MODXREF)
fillchrt.o    : $(MODINFO)
fillctbl.o    : $(MODXREF)
filletbl.o    : $(MODXREF)
fillgtbl.o    : $(MODXREF)
fillmtbl.o    : $(MODXREF)
fillptbl.o    : $(MODXREF)
fillstbl.o    : $(MODXREF)
fillttbl.o    : $(MODXREF)
fltrxref.o    : $(MODSOURC) $(MODLISTS) $(MODINFO)
fmtcsrc.o     : $(MODINFO)
genptvcel.o   : $(MODGRID)  $(MODGRDLIB)
genuslst.o    : $(MODSOURC) $(MODLISTS) $(MODINFO)
getbasyr.o    : $(MODSOURC)
getctgry.o    : $(MODINFO)
getntisz.o    : $(MODLISTS)
getsinfo.o    : $(MODFILESET) $(MODINFO)
gettzone.o    : $(MODSTCY)
grd2cnty.o    : $(MODSTCY) $(MODSURG) $(MODGRID)
grdfips.o     : $(MODSURG) $(MODGRID) $(MODMET)
normtpro.o    : $(MODTMPRL)
prclinrc.o    : $(MODREPRT) $(MODINFO)
procspro.o    : $(MODSPRO)
rdar2pt.o     : $(MODXREF) $(MODLISTS) $(MODAR2PT)
rdcodnam.o    : $(MODLISTS)
rddates.o     : $(MODTMPRL)
rdeproc.o     : $(MODEMFAC) $(MODINFO)
rdgeocodes.o  : $(MODSTCY)
rdgref.o      : $(MODXREF) $(MODINFO)
rdhdays.o     : $(MODTMPRL)
rdinvchr.o    : $(MODSOURC)
rdinvmap.o    : $(MODINFO)
rdinvpol.o    : $(MODFILESET)
rdmactdsc.o   : $(MODLISTS)
rdmappol.o    : $(MODINFO) $(MODFILESET)
rdnaicsdsc.o  : $(MODLISTS)
rdorsdsc.o    : $(MODSTCY)
rdpelv.o      : $(MODELEV)
rdsccdsc.o    : $(MODLISTS)
rdsccmap.o    : $(MODMOBIL)
rdsconv.o     : $(MODLISTS) $(MODSPRO) $(MODINFO)
rdsicdsc.o    : $(MODLISTS)
rdsprofdsc.o  : $(MODSOURC)
rdsmat.o      : $(MODFILESET)
rdspdprof.o   : $(MODMOBIL)
rdspdref.o    : $(MODXREF)
rdsprof.o     : $(MODSPRO)  $(MODLISTS)
rdsrcgrps.o   : $(MODMERGE) $(MODLISTS) $(MODGRID) $(MODELEV) $(MODSURG) $(MODSOURC)
rdsref.o      : $(MODXREF)  $(MODINFO)
rdsrgdesc.o   : $(MODSURG)
rdsrg.o       : $(MODSURG) $(MODGRID)
rdstcy.o      : $(MODSTCY)
rdtref.o      : $(MODXREF) $(MODINFO) $(MODTMPRL)
rdxclude.o    : $(MODXREF) $(MODINFO)
setscctype.o  : $(MODINFO)
setsrcdy.o    : $(MODINFO)
srcgrpcnt.o   : $(MODMERGE) $(MODGRID) $(MODELEV)
tagtable.o    : $(MODXREF)  $(MODSPRO) $(MODTAG) $(MODINFO)
wrchrscc.o    : $(MODLISTS) $(MODINFO)
wrorlout.o    : $(MODSOURC) $(MODSTCY) $(MODLISTS) $(MODINFO)
wrsrcgrps.o   : $(MODMERGE) $(MODGRID) $(MODELEV) $(MODGRDLIB)
xfractbl.o    : $(MODSPRO)
xreftbl.o     : $(MODXREF) $(MODINFO)

modgrdlib.o   : $(MODGRID)

chkfileset.o  : $(MODFILESET)
chksetdesc.o  : $(MODFILESET)
cleanup.o     : $(MODFILESET)
closeset.o    : $(MODFILESET)
createset.o   : $(MODFILESET)
descset.o     : $(MODFILESET)
modfileset.o  : $(MODFILESET)
openset.o     : $(MODFILESET)
promptset.o   : $(MODFILESET)
readset.o     : $(MODFILESET)
writeset.o    : $(MODFILESET)

normbeis314.o : $(MODBEIS3) $(MODGRDLIB)
normbeis361.o : $(MODBEIS3) $(MODGRDLIB)
tmpbeis314.o  : $(MODSPRO) $(MODBEIS3)
tmpbeis361.o  : $(MODSPRO) $(MODBEIS3)

alocpkts.o    : $(MODCNTRL)
asgncntl.o    : $(MODSOURC) $(MODXREF)
asgncntl.o    : $(MODINFO)  $(MODCNTRL)
cntlmat.o     : $(MODINFO)  $(MODFILESET)
fillcdat.o    : $(MODCNTRL) $(MODXREF)  $(MODINFO)
genmultc.o    : $(MODSOURC) $(MODCNTRL) $(MODINFO)
genproj.o     : $(MODSOURC) $(MODCNTRL) $(MODINFO)
genreact.o    : $(MODSOURC) $(MODCNTRL) $(MODINFO) $(MODSPRO)
opencmat.o    : $(MODCNTRL) $(MODINFO)
opencmat.o    : $(MODFILESET)
openpmat.o    : $(MODINFO)  $(MODCNTRL) $(MODFILESET)
openrmat.o    : $(MODINFO)  $(MODFILESET)
pktloop.o     : $(MODXREF)  $(MODINFO)
procpkts.o    : $(MODSOURC) $(MODXREF) $(MODCNTRL) $(MODINFO)
rdpacket.o    : $(MODINFO)
wcntlrep.o    : $(MODSOURC) $(MODCNTRL) $(MODINFO)
wrctmp.o      : $(MODINFO)

asgnbins.o    : $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODGRID)
asgnbins.o    : $(MODELEV)  $(MODSTCY)  $(MODINFO)  $(MODREPRT) $(MODREPBN) $(MODTMPRL)
bldrepidx.o   : $(MODCNTRL) $(MODINFO)  $(MODFILESET)
genrprt.o     : $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODTMPRL) $(MODCNTRL) $(MODINFO)
openrepin.o   : $(MODREPRT) $(MODTMPRL) $(MODREPBN) $(MODCNTRL) $(MODINFO)  $(MODFILESET)
openrepout.o  : $(MODREPRT)
qarepin.o     : $(MODREPRT) $(MODSTCY)  $(MODINFO)
rdgrps.o      : $(MODREPRT) $(MODGRID)  $(MODSTCY)
rdrepin.o     : $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODCNTRL) $(MODELEV)
rdrepin.o     : $(MODLISTS) $(MODGRID)  $(MODINFO)  $(MODFILESET)
rdrprts.o     : $(MODREPRT) $(MODINFO)
rdssup.o      : $(MODREPRT) $(MODINFO)  $(MODSOURC) $(MODSPRO)
repmrggrd.o   : $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODGRID)
repunits.o    : $(MODREPRT) $(MODREPBN) $(MODTMPRL) $(MODINFO)
scanrepc.o    : $(MODREPRT) $(MODINFO)
selectsrc.o   : $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODELEV)  $(MODINFO) $(MODSPRO)
smkreport.o   : $(MODREPRT) $(MODREPBN) $(MODGRID)  $(MODINFO)
wrrephdr.o    : $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODSTCY) $(MODGRID) $(MODINFO)
wrrepout.o    : $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODSTCY) $(MODINFO)

aggwndw.o     : $(MODGRID)
bluesky2inv.o : $(MODSTCY)
gentpro.o     : $(MODINFO) $(MODSOURC) $(MODGRID) $(MODSURG) $(MODSTCY)
geofac.o      : $(MODFILESET)
inlinto2d.o   : $(MODGRID) $(MODGRDLIB)
invsplit.o    : $(MODINFO)
smk2emis.o    : $(MODINFO)
surgtool.o    : $(MODSURG) $(MODGRID)

asgnsurg.o    : $(MODSOURC) $(MODXREF) $(MODINFO) $(MODSURG)
genagmat.o    : $(MODSOURC) $(MODGRID) $(MODINFO) $(MODSURG) $(MODLISTS) $(MODXREF) $(MODGRDLIB)
genggmat.o    : $(MODSOURC) $(MODGRID) $(MODINFO) $(MODGRDLIB)
genlgmat.o    : $(MODSOURC) $(MODGRID) $(MODINFO) $(MODGCTP) $(M3UTILIO)
genmgmat.o    : $(MODSOURC) $(MODGRID) $(MODINFO) $(MODSURG) $(MODLISTS) $(MODXREF)
grdmat.o      : $(MODSOURC) $(MODGRID) $(MODINFO) $(MODXREF) $(MODSURG)  $(MODFILESET) $(MODGRDLIB)
opengmat.o    : $(MODINFO)
rdsrg4grd.o   : $(MODSURG)  $(MODGRID)
setfrac.o     : $(MODSURG)  $(MODXREF)
sizgmat.o     : $(MODSURG)  $(MODGRID) $(MODGRDLIB)

mbldmrgidx.o  : $(MODMERGE) $(MODMVSMRG)
bldprocidx.o  : $(MODMERGE) $(MODMVSMRG)
bldsrccell.o  : $(MODMVSMRG)
mgetmrgev.o   : $(MODMERGE) $(MODINFO) $(MODMVSMRG)
hourmet.o     : $(MODMBSET) $(MODMET)  $(MODSURG)
minitstcy.o   : $(MODMERGE) $(MODSTCY) $(MODSOURC) $(MODMVSMRG) $(MODLISTS)
met4moves.o   : $(MODMBSET) $(MODGRID) $(MODSURG)  $(MODSTCY)   $(MODMET)
movesmrg.o    : $(MODMERGE) $(MODGRID) $(MODSTCY)  $(MODLISTS)  $(MODMBSET)   $(MODSOURC)$(MODMVSMRG)
mmrgonams.o   : $(MODMERGE)
mmrgunits.o   : $(MODMERGE) $(MODMVSMRG)
mmrgvnams.o   : $(MODMERGE) $(MODLISTS) $(MODMVSMRG)
mopenmrgin.o  : $(MODMERGE) $(MODINFO)  $(MODSOURC)  $(MODEMFAC) $(MODGRID) $(MODMVSMRG) $(MODFILESET)
mopenmrgout.o : $(MODMERGE) $(MODGRID)  $(MODMVSMRG) $(MODFILESET)
rdcfpro.o     : $(MODMERGE) $(MODLISTS) $(MODMBSET)  $(MODMVSMRG)
rdfmref.o     : $(MODMBSET)
rdmrclist.o   : $(MODMBSET) $(MODMVSMRG)
rdmxref.o     : $(MODMBSET)
rdrpdemfacs.o : $(MODMBSET)  $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
rdrphemfacs.o : $(MODMBSET)  $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
rdrppemfacs.o : $(MODMBSET)  $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
rdrpvemfacs.o : $(MODMBSET)  $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
rdspdpro.o    : $(MODMVSMRG) $(MODLISTS)  $(MODMOBIL)
rdspdist.o    : $(MODMVSMRG) $(MODLISTS)  $(MODMOBIL)
setrefcnty.o  : $(MODMVSMRG) $(MODINFO)   $(MODSOURC) $(MODLISTS) $(MODMBSET)
mwrmrgrep.o   : $(MODMERGE)  $(MODSTCY)   $(MODLISTS) $(MODMVSMRG)
wrmrgrep.o    : $(MODGRID)
wrtemprof.o   : $(MODMET)  $(MODMBSET)

asgngrps.o    : $(MODSOURC) $(MODELEV) $(MODINFO) $(MODGRDLIB)
elevpoint.o   : $(MODSOURC) $(MODELEV) $(MODINFO) $(MODGRID) $(MODGRDLIB)
laypoint.o    : $(MODSOURC) $(MODELEV) $(MODINFO) $(MODGRID) $(MODGRDLIB) $(MODFILESET)
mxgrpemis.o   : $(MODSOURC) $(MODELEV) $(MODINFO) $(MODGRID) $(MODGRDLIB) $(MODFILESET)
openeout.o    : $(MODELEV)  $(MODINFO)
openlayout.o  : $(MODELEV)  $(MODINFO) $(MODGRID)
rpelvcfg.o    : $(MODREPRT) $(MODELEV) $(MODINFO)
wpingstk.o    : $(MODELEV)

asgnar2pt.o   : $(MODSOURC) $(MODELEV) $(MODINFO)
asgnnhapx.o   : $(MODSOURC) $(MODELEV) $(MODINFO)
fixstk.o      : $(MODSOURC)
formlist.o    : $(MODLISTS) $(MODINFO)
genmedsout.o  : $(MODINFO)  $(MODDAYHR)
genpdout.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODSTCY)
gethdr.o      : $(MODLISTS) $(MODSTCY)  $(MODINFO)
getpdinfo.o   : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
grwinven.o    : $(MODSOURC) $(MODINFO)  $(MODFILESET)
initinfo.o    : $(MODINFO)
opengrwout.o  : $(MODINFO)
openinvin.o   : $(MODINFO)  $(MODLISTS) $(MODXREF) $(MODDAYHR) $(MODMOBIL)
openinvout.o  : $(MODINFO)  $(MODFILESET)
openpdout.o   : $(MODINFO)
procar2pt.o   : $(MODSOURC) $(MODXREF)  $(MODINFO) $(MODAR2PT)
procinven.o   : $(MODSOURC) $(MODLISTS) $(MODINFO)
procinvsrcs.o : $(MODSOURC) $(MODINFO)
rdcempd.o     : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
rdcemsum.o    : $(MODDAYHR)
rddataff10ar.o: $(MODINFO) $(MODDAYHR)
rddataff10mb.o: $(MODINFO) $(MODLISTS)
rddataff10pt.o: $(MODINFO) $(MODDAYHR)
rddatamedspt.o: $(MODINFO) $(MODSOURC)
rddatantiar.o : $(MODINFO)
rddatantifr.o : $(MODINFO)
rddatantimb.o : $(MODINFO)
rddatantinp.o : $(MODINFO)
rddatantipt.o : $(MODINFO)
rdemspd.o     : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
rdff10pd.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODSTCY)
rdgrdapi.o    : $(MODSOURC) $(MODINFO)
rdgrdncf.o    : $(MODSOURC) $(MODINFO)
rdinvdata.o   : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODMOBIL)
rdinvsrcs.o   : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODMOBIL)
rdlooppd.o    : $(MODINFO)  $(MODLISTS)
rdmedsinfo.o  : $(MODSOURC)
rdmedspd.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODSTCY)
rdorlfr.o     : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
rdsrcff10ar.o : $(MODLISTS) $(MODDAYHR)
rdsrcff10mb.o : $(MODLISTS)
rdsrcff10pt.o : $(MODLISTS) $(MODDAYHR)
rdsrcmedspt.o : $(MODSOURC)
rdsrcntiar.o  : $(MODLISTS)
rdsrcntifr.o  : $(MODINFO) $(MODLISTS)
rdsrcntimb.o  : $(MODLISTS)
rdsrcntinp.o  : $(MODLISTS)
rdsrcntipt.o  : $(MODLISTS)
setnonhap.o   : $(MODSOURC) $(MODLISTS) $(MODXREF) $(MODINFO)
smkinven.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
srcmem.o      : $(MODSOURC)
wrepinven.o   : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODAR2PT)
wrinvchr.o    : $(MODSOURC) $(MODLISTS) $(MODINFO)
wrinvemis.o   : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODFILESET)
wrinvpol.o    : $(MODFILESET)
wrpdemis.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
wrptref.o     : $(MODINFO)

allocmrg.o    : $(MODMERGE) $(MODELEV) $(MODCNTRL) $(MODSTCY) $(MODGRID)
bldmrgidx.o   : $(MODMERGE)
getmrgev.o    : $(MODMERGE) $(MODINFO)
initstcy.o    : $(MODMERGE) $(MODSTCY)
mrgelev.o     : $(MODMERGE) $(MODCNTRL) $(MODELEV)
mrggrid.o     : $(MODGRID)  $(MODFILESET)
mrgmult.o     : $(MODMERGE) $(MODELEV)  $(MODSTCY)
mrgonams.o    : $(MODMERGE)
mrgunits.o    : $(MODMERGE)
mrgvnams.o    : $(MODMERGE) $(MODLISTS) 
openmrgin.o   : $(MODMERGE) $(MODINFO)  $(MODELEV) $(MODGRID) $(MODFILESET)
openmrgout.o  : $(MODMERGE) $(MODELEV)  $(MODGRID) $(MODFILESET)
rdmrginv.o    : $(MODMERGE) $(MODLISTS) $(MODSOURC)
smkmerge.o    : $(MODMERGE) $(MODCNTRL) $(MODELEV) $(MODLISTS) $(MODSTCY) $(MODSURG) $(MODGRID)
wmrgelev.o    : $(MODMERGE) $(MODSOURC) $(MODELEV) $(MODGRID)
wmrgemis.o    : $(MODMERGE) $(MODELEV)
wrmrgrep.o    : $(MODMERGE) $(MODSTCY)  $(MODGRID)

asgnspro.o    : $(MODLISTS) $(MODXREF)  $(MODSPRO) $(MODINFO)
asgntag.o     : $(MODSOURC) $(MODTAG)   $(MODINFO)
chknonhap.o   : $(MODLISTS) $(MODSPRO)
opensmat.o    : $(MODINFO)  $(MODTAG)   $(MODFILESET)
rdcombo.o     : $(MODSPRO)  $(MODLISTS)
rdtag.o       : $(MODXREF)  $(MODSPRO)  $(MODINFO) $(MODLISTS) $(MODTAG)
spcmat.o      : $(MODSPRO)  $(MODEMFAC) $(MODINFO) $(MODLISTS) $(MODTAG) $(MODFILESET)

asgntpro.o    : $(MODSOURC) $(MODTMPRL) $(MODINFO)
genhemis.o    : $(MODSOURC) $(MODXREF)  $(MODINFO) $(MODTMPRL) $(MODDAYHR)
mktmat.o      : $(MODSOURC) $(MODTMPRL) $(MODDAYHR)
opentmp.o     : $(MODINFO)  $(MODTMPRL) $(MODFILESET)
opentmpin.o   : $(MODINFO)
proctpro.o    : $(MODSOURC) $(MODTMPRL) $(MODINFO)
setdaylt.o    : $(MODSOURC) $(MODSTCY)  $(MODINFO)
temporal.o    : $(MODSOURC) $(MODTMPRL) $(MODINFO) $(MODLISTS) $(MODEMFAC) $(MODDAYHR) $(MODMBSET)
tmnamunt.o    : $(MODINFO)
updtmat.o     : $(MODSOURC) $(MODTMPRL) $(MODDAYHR)
wrtsup.o      : $(MODSOURC) $(MODTMPRL) $(MODINFO)


######################################################################
#                            builds:

${OBJDIR}:
	echo "Making " ${OBJDIR} ; mkdir -p ${OBJDIR}

##############################  Libraries

${MLIB}:  ${MODOBJ}
	cd ${OBJDIR}; $(AR) ${ARFLAGS} ${MLIB} ${MODOBJ}

${SLIB}: ${LIBOBJ}
	cd ${OBJDIR}; $(AR) ${ARFLAGS} ${SLIB} ${LIBOBJ}

${FLIB}: ${FILOBJ}
	cd ${OBJDIR}; echo ${FILOBJ}; $(AR) ${ARFLAGS} ${FLIB} ${FILOBJ}

##############################  SMOKE programs

tmpbeis3:  czangle.o getpar.o getparb.o hrbeis_361.o hrbeis3.o prebmet.o rdbpro.o \
 tmpbeis3.o tmpbeis314.o tmpbeis361.o hrno.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

normbeis3: normbeis3.o rdb3fac.o rdb4fac.o rdb4fac_csv.o normbeis314.o normbeis361.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

cntlmat: $(CTLOBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

smkreport: $(RPTOBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

grdmat: $(GRDOBJ) 
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

movesmrg: movesmrg.o mbldmrgidx.o bldprocidx.o bldsrccell.o \
 mgetmrgev.o minitstcy.o mmrgonams.o mmrgunits.o  mmrgvnams.o   \
 mopenmrgin.o mopenmrgout.o rdfmref.o  rdmrclist.o rdmxref.o    \
 rdrpdemfacs.o rdrppemfacs.o rdrpvemfacs.o rdrphemfacs.o rdspdpro.o     \
 rdcfpro.o setrefcnty.o mwrmrgrep.o rdspdist.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) $^ ${LIBS} -o $@

met4moves: met4moves.o rdmxref.o rdfmref.o hourmet.o wrtemprof.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) $^ ${LIBS} -o $@

laypoint: laypoint.o openlayout.o preplm.o plmris.o postplm.o plsprd.o \
          fire_preplm.o fire_plmris.o fire_postplm.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

elevpoint: elevpoint.o asgngrps.o openeout.o plumesti.o wpingstk.o  \
           mxgrpemis.o rpelvcfg.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

smkinven: adjustinv.o asgnar2pt.o asgnnhapx.o chklstfl.o fixstk.o formlist.o \
 genmedsout.o genpdout.o gethdr.o getpdinfo.o initinfo.o      \
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

smkmerge: smkmerge.o allocmrg.o bldmrgidx.o getmrgev.o mrgmult.o mrgonams.o \
 mrgvnams.o openmrgin.o openmrgout.o wmrgemis.o wrmrgrep.o initstcy.o \
 rdmrginv.o mrgelev.o mrgunits.o mrgbio.o wmrgelev.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

mrggrid: mrggrid.o setoutdate.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

mrgelev: mrgasciielev.o setoutdate.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

mrgpt: mrgpt.o setoutdate.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)
 
spcmat: $(SPCOBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

temporal: $(TMPOBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

##############################  emutil programs:

aggwndw: aggwndw.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

beld3to2: beld3to2.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

bluesky2inv: bluesky2inv.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

cemscan: cemscan.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

extractida: extractida.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

gcntl4carb: gcntl4carb.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

gentpro: gentpro.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

geofac: geofac.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

inlineto2d: inlineto2d.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

invsplit: invsplit.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

layalloc: layalloc.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

metcombine: metcombine.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

metscan: metscan.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

pktreduc: pktreduc.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

saregroup: saregroup.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

smk2emis: smk2emis.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

surgtool: surgtool.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

uam2ncf: uam2ncf.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)


