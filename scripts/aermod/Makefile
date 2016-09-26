IODIR=../../../ioapi/Linux2_x86_64ifort

convert_phour: convert_phour.f90
	ifort -o convert_phour convert_phour.f90 -I$(IODIR) -L$(IODIR) -lioapi -lnetcdf -qopenmp
