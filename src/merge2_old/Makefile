# Version "@(#)$Id$ $Source$ $Date$ 
#
all: opts
#
opts:   csgmerge csglmerge csglmerge2
#
LDIRS=-L$(em) \
      -L$(io) \
      -L$(ncf)
#
LIBS=$(LDIRS) -lemstuff -lm3io -lnetcdf
#
FFLAGS = -C
#
#
CSGS=csgmerge.F  amrgvec.F pmrgvec.F msmatvec.F

CSGI=ARDIMS3.EXT PTDIMS3.EXT GRDIMS3.EXT CHDIMS3.EXT MBDIMS3.EXT BIODIMS3.EXT
#
csgmerge: $(CSGS) $(CSGI)
	echo CSGLMERGE:
	f77 $(FFLAGS) -o $(SCENDIR)/csgmerge $(CSGS) $(LIBS)
#
#
#
CSGLS=csglmerge.F  amrgvec.F  plmvec.F   msmatvec.F

CSGLI=ARDIMS3.EXT PTDIMS3.EXT LADIMS3.EXT GRDIMS3.EXT CHDIMS3.EXT \
      MBDIMS3.EXT BIODIMS3.EXT
#
csglmerge: $(CSGLS) $(CSGLI)
	f77 $(FFLAGS) -o $(SCENDIR)/csglmerge $(CSGLS) $(LIBS)
#
#
#
CSGLS2=csglmerge.2.F amrgvec.F plmvec.F msmatvec.F csvname.F
#
csglmerge2: $(CSGLS2) $(CSGLI)
	f77 $(FFLAGS) -o $(SCENDIR)/csglmerge2 $(CSGLS2) $(LIBS)
#
