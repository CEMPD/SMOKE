# Version "@(#)$Id$ $Source$ $Date$ 
#
all: getrecs smk2emis smk2pts mkmplist
#
# compile-link options for subroutines
#
LDIRS=-L$(em) -L$(io) -L$(ncf) -L$(eps)
#
LIBS=$(LDIRS) -lemstuff -lm3io -lnetcdf -leps2
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
EMIS = smk2emis.F trans.F chkgrd.F chktmp.F  chkvar.F

EINC = PARMS3.EXT FDESC3.EXT GRDIMS3.EXT IODECL3.EXT PTDIMS3.EXT \
       uamcntl.inc uamregn.inc param.inc mrguam.inc

smk2emis: $(EMIS) $(EINC)
	echo SMK2EMIS:
	f77 $(FFLAGS) -o $(SCENDIR)/smk2emis $(EMIS) $(LIBS)
#
#
#
PTS = smk2pts.F chktmp.F chkvar.F rdepsd.F
#
smk2pts: $(OBJS) $(EINC)
	echo SMK2PTS:
	f77 $(FFLAGS) -o $(SCENDIR)/smk2pts $(PTS) $(LIBS)
#
#
#
MKMP = mkmplist.F

MKMI = PARMS3.EXT FDESC3.EXT GRDIMS3.EXT IODECL3.EXT PTDIMS3.EXT
#
mkmplist: $(MKMP) $(MKMI)
	echo MKPMLIST:
	f77 $(FFLAGS) -o $(SCENDIR)/mkmplist $(MKMP) $(LIBS)
#
#
clean:
	rm *.o
