#!/bin/tcsh -f

## Assuming SMOKE and I/O API were compiled using instruction in https://github.com/CEMPD/SMOKE/wiki/B.-Instructions-for-SMOKE-Installation
setenv BIN Linux2_x86_64ifx
source $SMKDEV/ioapi-3.2.github/IOAPI.config.csh $BIN  

## Modify the following settings as needed, espcially if SMOKE and I/O API were not compiles using the above instruction
## Also comment out the corresponding settings in Makefile and Makefile.wrf
setenv NETCDFINC "-I${NETCDFC}/include -I${NETCDFF}/include"
setenv IOBASE ${IOAPI_HOME}
setenv IODIR  ${IOBASE}/ioapi
#setenv IOBIN  ${IOAPI_HOME}/$BIN
setenv IOINC  ${IODIR}/fixed_src
setenv IOLIB  "-L${IOBIN} -lioapi ${NCFLIBS}"
 
## Commment out correspoding settings the following 
setenv FC     'ifx'
setenv LFLAGS '-auto -warn -extend-source 132 -zero -check bounds -traceback -qopenmp'
setenv LDIRS  '-L/usr/lib'
setenv LLIBS  '${IOLIB} -lunwind'
setenv IFLAGS '-I${IOINC} ${NETCDFINC} -I${IOBIN} -module ${IOBIN}'
setenv EFLAG  '-auto -warn -extend-source 132 -zero -check bounds -traceback -qopenmp'

#> compile apply_precip_adj
make all -f Makefile

#> compile apply_precip_adj_wrf
setenv FFLAGS "$LFLAGS $IFLAGS"
make all -f Makefile.wrf


