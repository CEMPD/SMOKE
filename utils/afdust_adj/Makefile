.SUFFIXES:
.SUFFIXES: .f .F .o .x


## Assumes you have already sourced ASSIGN file used for making SMOKE
## this will supply FC, FFLAGS, LDIRS

#IOBASE = /work/EMIS/smoke/test/repo/ioapi-3.2
#IODIR   = ${IOBASE}/ioapi
#IOBIN   = /work/EMIS/smoke/test/ioapi
#IOINC   = ${IODIR}/fixed_src
#IOLIB   = ${IOBIN}

#NCFLIB = /work/EMIS/smoke/test/repo/netcdf/lib/


#FC=ifort
#FFLAGS= -O -132 -check -debug -g 
#LFLAGS=  -g -debug
#LIBDIR= -L/usr/lib  -L${IOAPIDIR} -L${NETCDFDIR} 
#LLIBS=  -lioapi -lnetcdf -lunwind 
#IFLAGS = -I${IOINC} -I${INCDIR} -I${IOBIN}  #  Intel, Absoft, Pathscale, UNIX Fortran
#EFLAG = -extend-source 132 -auto -zero -xHost -traceback #  Intel Fortran

FFLAGS  = ${IFLAGS} ${EFLAG} ${DEFINEFLAGS} ${PARFLAGS} ${FOPTFLAGS} ${ARCHFLAGS}
LDFLAGS = ${IFLAGS} ${DEFINEFLAGS} ${ARCHFLAGS}


SRC = .
BIN = .

SRCNAME = apply_precip_adj

OBJS = $(SRCNAME).o

all: $(OBJS) $(SRCNAME).x

.F.o:
	$(FC) $(FFLAGS) -c $(SRC)/$<

.f.o:
	$(FC) $(FFLAGS) -c $(SRC)/$<

.o.x:
	$(FC) $(LFLAGS) -o $(BIN)/$@ $< $(LDIRS) $(LLIBS)


clean:
	- rm -f *.o
	- rm $(BIN)/*.x
