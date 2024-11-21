.SUFFIXES:
.SUFFIXES: .f .F .o .x


## Assumes you have already sourced ASSIGN file used for making SMOKE
## this will supply FC, FFLAGS, LDIRS

IOBASE = /work/EMIS/smoke/test/repo/ioapi-3.2
IODIR   = ${IOBASE}/ioapi
IOBIN   = /work/EMIS/smoke/test/ioapi
IOINC   = ${IODIR}/fixed_src
IOLIB   = ${IOBIN}

NCFLIB = /work/EMIS/smoke/test/repo/netcdf/lib/


FC=ifort
LFLAGS=  -g -debug
LIBDIR= -L/usr/lib  -L${IOBIN} -L${NCFLIB} 
LLIBS=  -lioapi -lnetcdff -lnetcdf #-lunwind 
IFLAGS = -I${IOINC} -I${IOBIN} #  Intel, Absoft, Pathscale, UNIX Fortran
EFLAG = -extend-source 132 -auto -zero -xHost -traceback -Bstatic #  Intel Fortran
FFLAGS  = ${IFLAGS} ${EFLAG} ${DEFINEFLAGS} ${PARFLAGS} ${FOPTFLAGS} ${ARCHFLAGS}

LDFLAGS = ${IFLAGS} ${DEFINEFLAGS} ${ARCHFLAGS}


SRC = .
BIN = .

SRCNAME = apply_precip_adj_wrf

OBJS = $(SRCNAME).o

all: $(OBJS) $(SRCNAME).x


.f.o:
	$(FC) $(FFLAGS) -c $(SRC)/$< 

.o.x:
	$(FC) $(FFLAGS) -o $(BIN)/$@ $< $(LIBDIR) $(LLIBS)


clean:
	- rm -f *.o
	- rm $(BIN)/*.x
