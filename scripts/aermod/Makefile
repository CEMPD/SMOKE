IODIR=../../../ioapi/Linux2_x86_64ifort

all: convert_phour cmv_output_hourly

convert_phour: convert_phour.f90
	ifort -o convert_phour convert_phour.f90 -I$(IODIR) -L$(IODIR) -lioapi -lnetcdff -lnetcdf -qopenmp

cmv_output_hourly: cmv_output_hourly.f90
	ifort -o cmv_output_hourly cmv_output_hourly.f90 -I$(IODIR) -L$(IODIR) -lioapi -lnetcdff -lnetcdf -qopenmp
