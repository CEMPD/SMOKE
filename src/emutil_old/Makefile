# Version "@(#)$Id$ $Source$ $Date$ 
#
# Note: Programs smk2emis and smk2pts require library eps2 which is not
#       provided with the standard SMOKE distribution.
#
all: faketemp mappoint mkmplist srggrid smk2emis smk2pts 
#
# compile-link options for subroutines
#
LIBS=$(LDIRS) -lemlib -lioapi -lnetcdf
#
LIBS2=$(LDIRS) $(LIBS) -leps2
#
#
#
FTMP=faketemp.F
#
FTMPI=PARMS3.EXT FDESC3.EXT IODECL3.EXT GRDIMS3.EXT
#
faketemp: $(FTMP) $(FTMPI)
	echo FAKETEMP:
	f77 $(FFLAGS)   -o $(SCENDIR)/faketemp $(FTMP) $(LIBS)
#
#
MAP= mappoint.F
#
MAPI= PTDIMS3.EXT CHDIMS3.EXT GRDIMS3.EXT
#
mappoint: $(MAP) $(MAPI)
	echo MAPPOINT:
	f77 $(FFLAGS)   -o $(SCENDIR)/mappoint $(MAP) $(LIBS)
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
#
EMIS = smk2emis.F trans.F chkgrd.F chktmp.F  chkvar.F

SHAREDINC = $(SMKROOT)/inc/Shared

EINC = PARMS3.EXT FDESC3.EXT GRDIMS3.EXT IODECL3.EXT PTDIMS3.EXT \
       $(SHAREDINC)/uamcntl.inc $(SHAREDINC)/uamregn.inc \
	   $(SHAREDINC)/param.inc $(SHAREDINC)/mrguam.inc
	   
smk2emis: $(EMIS) $(EINC)
	echo SMK2EMIS:
	f77 $(FFLAGS) -I$(SHAREDINC) -o $(SCENDIR)/smk2emis $(EMIS) $(LIBS2)
#
#
#
PTS = smk2pts.F chktmp.F chkvar.F rdepsd.F
#
smk2pts: $(OBJS) $(EINC)
	echo SMK2PTS:
	f77 $(FFLAGS) -I$(SHAREDINC) -o $(SCENDIR)/smk2pts $(PTS) $(LIBS2)
#
#
#
SRGP = srggrid.F

SRGI = ARDIMS3.EXT GRDIMS3.EXT
#
srggrid: $(SRGP) $(SRGI)
	echo SRGGRID:
	f77 $(FFLAGS) -o $(SCENDIR)/srggrid $(SRGP) $(LIBS)
#
#
#
clean:
	rm *.o
