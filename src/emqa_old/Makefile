# Version "@(#)$Id$ $Source$ $Date$ 
#
# Programs listed in alphabetical order
#
all: getrecs pntsum repplht reptpro
#
#
LIBS=$(LDIRS) -lemstuff -lm3io -lnetcdf
#
LIBS2=$(LDIRS) -lemstuff -lm3io -lnetcdf -leps2
#
#
#
OBJS = getrecs.F

INCS = PARMS3.EXT FDESC3.EXT GRDIMS3.EXT IODECL3.EXT PTDIMS3.EXT
 
getrecs: $(OBJS) $(INCS)
	echo GETRECS:
	f77 $(FFLAGS) -o $(SCENDIR)/getrecs $(OBJS) $(LIBS)
#
#
#
PTS=pntsum.F
#
PTSI=PTDIMS3.EXT CHDIMS3.EXT 
#
pntsum: $(PTS) $(PTSI)
	echo PNTSUM:
	f77 $(FFLAGS)   -o $(SCENDIR)/$@ $(PTS) $(LIBS)
#
#
#
RPH=repplht.F
#
RPHI=PTDIMS3.EXT CHDIMS3.EXT TMDIMS3.EXT VGRD3.EXT
#
repplht: $(RPH) $(RPHI)
	echo REPPLHT:
	f77 $(FFLAGS)   -o $(SCENDIR)/$@ $(RPH) $(LIBS)
#
#
#
RTP=reptpro.F
#
RTPI=PTDIMS3.EXT CHDIMS3.EXT TMDIMS3.EXT 
#
reptpro: $(RTP) $(RTPI)
	echo REPTPRO:
	f77 $(FFLAGS)   -o $(SCENDIR)/$@ $(RTP) $(LIBS)
#
#
#
clean:
	rm *.o
