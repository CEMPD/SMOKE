# Version "@(#)$Id$ $Source$ $Date$ 
#
all: opts
#
opts: emisfac grdmmat prediur rawmobil spcmmat tmpmobil adjmobil
#
LIBS=$(LDIRS) -lemlib -lioapi -lnetcdf 
#
LMOB=$(LDIRS) -lemlib -lioapi -lnetcdf -lmo5b
#
#
#
EMIS=emisfac.F efcombo.F moinbld.F mobil5b.F

EMISI=MBDIMS3.EXT CHDGEN3.EXT

emisfac: $(EMIS) $(EMISI)
	echo EMISFAC:
	f77 $(FFLAGS)  -I$(SMKROOT)/lib/mo5b -o $(SCENDIR)/emisfac $(EMIS) $(LMOB)
#
#
#
GRMS=grdmmat.F lnk2grd.F rdmvadj.F adjmv.F
#
GRMI=MBDIMS3.EXT GRDIMS3.EXT
#
grdmmat: $(GRMS) $(GRMI)
	echo GRDMMAT:
	f77 $(FFLAGS) -o $(SCENDIR)/grdmmat $(GRMS) $(LIBS)
#
#
#
PRED=prediur.F rdmplist.F

PREDI=MBDIMS3.EXT GRDIMS3.EXT

prediur: $(PRED) $(PREDI)
	echo PREDIUR:
	f77 $(FFLAGS) -o $(SCENDIR)/prediur $(PRED) $(LIBS)
#
#
#
RAW=rawmobil.F 
#
RAWI=MBDIMS3.EXT TMDIMS3.EXT
#
rawmobil: $(RAW) $(RAWI)
	echo RAWMOBIL:
	f77 $(FFLAGS)  -o $(SCENDIR)/rawmobil $(RAW) $(LIBS)
#
#
#
SPCS=spcmmat.F csvname.F rdsprof.F rdmconv.F
#
SPCI=MBDIMS3.EXT CHDGEN3.EXT
#
spcmmat: $(SPC) $(SPCI)
	echo SPCMMAT:
	f77 $(FFLAGS) -o $(SCENDIR)/spcmmat $(SPCS) $(LIBS)
#
#
#
TMPV=tmpmobil.F mkvmat.F rdmplist.F csvname.F rdvmix.F
#
TMPI=MBDIMS3.EXT CHDGEN3.EXT TMDIMS3.EXT GRDIMS3.EXT
#
tmpmobil: $(TMPV) $(TMPI)
	echo TMPMOBIL:
	f77 $(FFLAGS) -o $(SCENDIR)/tmpmobil $(TMPV) $(LIBS)
#
#
#
ADJM=adjmobil.f ../lib/emlib/initem.F getfline.f padzero.f checkmem.f parsline.f
#
ADJI=MBDIMS3.EXT TMDIMS3.EXT
#
adjmobil: $(ADJM) $(ADJI)
	echo ADJMOBIL:
	f90 $(FFLAGS) -o $(SCENDIR)/adjmobil $(ADJM) \
         -L/hera/edss/subsys/smokev1/HP_UX11f90 -lioapi.debug -lnetcdf.debug
#
#
#
clean:
	rm *.o
