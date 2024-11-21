#!/bin/tcsh -f

setenv SMKBIN /proj/ie/proj/SMOKE/htran/SMOKE_INSTALL/Linux2_x86_64ifort # Path to compiled SMOKE executable
setenv SMKSRC /proj/ie/proj/SMOKE/htran/SMOKE_INSTALL/src                # Path to SMOKE src folder

setenv IOAPI_HOME /proj/ie/proj/SMOKE/htran/ioapi-3.2                    # Path to I/O API
setenv BIN Linux2_x86_64ifort                                            # If using intel fortran compiler (ifort) to compile I/O API

setenv NETCDFC ${IOAPI_HOME}/dependencies/netcdfc.4.9.2_intel            # Path to netcdf-C that was used to compile I/O API
setenv NETCDFF ${IOAPI_HOME}/dependencies/netcdff.4.6.1_intel            # Path to netcdf-fortran that was used to compile I/O API
setenv HDF5    ${IOAPI_HOME}/dependencies/hdf5.1.14.3_intel              # Optional: Path to HDF-5 that was used to compile I/O API with netcdf-4 support
setenv SZIP    ${IOAPI_HOME}/dependencies/szip.2.1.1_intel               # Optional: Path to szip library that was used to compile I/O API with netcdf-4 support 
setenv ZLIB    ${IOAPI_HOME}/dependencies/zlib.1.3_intel                 # Optional: Path to zlib library that was used to compile I/O API with netcdf-4 support

setenv srccode gen_afdust_tfrac_v2.f
setenv exe     gen_afdust_tfrac_v2

setenv depndnt "${SMKBIN}/getfline.o ${SMKBIN}/getcfdsc.o"
setenv FC /proj/ie/proj/SMOKE/htran/external_software/intel_oneapi/install/2024.1/bin/ifort

rm -rf $exe.o $exe

$FC -auto -warn notruncated_source -Bstatic -static-intel -qopenmp -qopenmp-link=static -I${SMKBIN} -I${IOAPI_HOME}/ioapi/fixed_src -I${SMKSRC}/inc -I${IOAPI_HOME}/${BIN} -extend-source 132 -zero -static-intel -debug -check bounds -traceback -O2 -unroll -stack-temps -safe-cray-ptr -convert big_endian -assume byterecl -DAUTO_ARRAYS=1 -DF90=1 -DFLDMN=1 -DFSTR_L=int -DIOAPI_NO_STDOUT=1 -DAVOID_FLUSH=1 -DBIT32=1  -c $srccode

$FC -auto -warn notruncated_source -Bstatic -static-intel -qopenmp -qopenmp-link=static -I${SMKBIN} -I${IOAPI_HOME}/ioapi/fixed_src -I${SMKSRC}/inc -I${IOAPI_HOME}/${BIN} -extend-source 132 -zero -static-intel -debug -check bounds -traceback -O2 -unroll -stack-temps -safe-cray-ptr -convert big_endian -assume byterecl -DAUTO_ARRAYS=1 -DF90=1 -DFLDMN=1 -DFSTR_L=int -DIOAPI_NO_STDOUT=1 -DAVOID_FLUSH=1 -DBIT32=1 -o $exe ${exe}.o $depndnt -L/usr/lib  -L${IOAPI_HOME}/${BIN} -lioapi -L${NETCDFF}/lib -lnetcdff -L${NETCDFC}/lib -L${HDF5}/lib -L${SZIP}/lib -L${ZLIB}/lib -lnetcdf -lhdf5_hl -lhdf5 -lm -lz -lsz -L${SMKBIN} -lfileset -lsmoke -lemmod
