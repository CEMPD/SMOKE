# Version "@(#)$Id$ $Source$ $Date$ 
#
all: opts
#
opts:   csgmerge csglmerge csgldaymrg
#
#
LIBS=$(LDIRS) -lemlib -lioapi -lnetcdf
#
#
#
CSGS=csgmerge.F  amrgvec.F pmrgvec.F

CSGI=ARDIMS3.EXT PTDIMS3.EXT GRDIMS3.EXT CHDIMS3.EXT MBDIMS3.EXT BIODIMS3.EXT
#
csgmerge: $(CSGS) $(CSGI)
	echo CSGLMERGE:
	f77 $(FFLAGS) -o $(SCENDIR)/csgmerge $(CSGS) $(LIBS)
#
#
#
CSGLS=csglmerge.F  amrgvec.F  plmvec.F

CSGLI=ARDIMS3.EXT PTDIMS3.EXT VGRD3.EXT GRDIMS3.EXT CHDIMS3.EXT \
      MBDIMS3.EXT BIODIMS3.EXT
#
csglmerge: $(CSGLS) $(CSGLI)
	f77 $(FFLAGS) -o $(SCENDIR)/csglmerge $(CSGLS) $(LIBS)
#
#
#
CSGLS2=csgldaymrg.F amrgvec.F plmvec.F csvname.F
#
CSGLI2=ARDIMS3.EXT PTDIMS3.EXT VGRD3.EXT GRDIMS3.EXT CHDIMS3.EXT \
       MBDIMS3.EXT BIODIMS3.EXT TMDIMS3.EXT CHDGEN3.EXT
#
csgldaymrg: $(CSGLS2) $(CSGLI2)
	f77 $(FFLAGS) -o $(SCENDIR)/csgldaymrg $(CSGLS2) $(LIBS)
#
