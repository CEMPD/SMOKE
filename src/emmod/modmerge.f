
        MODULE MODMERGE

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public variables and allocatable arrays 
!     used only in the merge module.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 3/99 by M. Houyoux
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
!
! COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
! All Rights Reserved
! 
! Carolina Environmental Program
! University of North Carolina at Chapel Hill
! 137 E. Franklin St., CB# 6116
! Chapel Hill, NC 27599-6116
! 
! smoke@unc.edu
!
! Pathname: $Source$
! Last updated: $Date$ 
!
!****************************************************************************

        IMPLICIT NONE

        INCLUDE 'EMPRVT3.EXT'

!.........  CONTROL VARIABLES...
!.........  Flags to indicate source categories being processed
        LOGICAL, PUBLIC :: AFLAG
        LOGICAL, PUBLIC :: BFLAG
        LOGICAL, PUBLIC :: MFLAG
        LOGICAL, PUBLIC :: PFLAG
        LOGICAL, PUBLIC :: XFLAG   ! multi-sources

!.........  Flags to indicate multiplicative controls being used
        LOGICAL, PUBLIC :: AUFLAG
        LOGICAL, PUBLIC :: MUFLAG
        LOGICAL, PUBLIC :: PUFLAG
        LOGICAL, PUBLIC :: TUFLAG  ! total

!.........  Flags to indicate reactivity controls being used
        LOGICAL, PUBLIC :: ARFLAG
        LOGICAL, PUBLIC :: MRFLAG
        LOGICAL, PUBLIC :: PRFLAG
        LOGICAL, PUBLIC :: TRFLAG  ! total

!.........  Flags to indicate inventory emissions are projected or not
        LOGICAL, PUBLIC :: PFACFLAG = .FALSE.     ! "pfac" variable flag
        LOGICAL, PUBLIC :: APRJFLAG = .FALSE.
        LOGICAL, PUBLIC :: MPRJFLAG = .FALSE.
        LOGICAL, PUBLIC :: PPRJFLAG = .FALSE.

!.........  Flags to indicate if by-day hourly emissions files are being used
        LOGICAL, PUBLIC :: AFLAG_BD = .FALSE.
        LOGICAL, PUBLIC :: MFLAG_BD = .FALSE.
        LOGICAL, PUBLIC :: PFLAG_BD = .FALSE.

!.........  Flags for controlling on-off features
        LOGICAL, PUBLIC :: TFLAG   ! use temporalized (hourly) emissions
        LOGICAL, PUBLIC :: SFLAG   ! use speciation
        LOGICAL, PUBLIC :: LFLAG   ! use layer fractions
        LOGICAL, PUBLIC :: HFLAG   ! create hourly outputs
        LOGICAL, PUBLIC :: PINGFLAG! create plume-in-grid hourly file
        LOGICAL, PUBLIC :: INLINEFLAG ! create in-line ptsource file
        LOGICAL, PUBLIC :: ELEVFLAG! create ASCII elevated sources file
        LOGICAL, PUBLIC :: EXPLFLAG! use PHOUR file to get explicit plume rise
        LOGICAL, PUBLIC :: SMATCHK ! use SMAT int output file for model species calc (Movesmrg only)
        LOGICAL, PUBLIC :: SUBSECFLAG ! use sub-sectorgrouping
        LOGICAL, PUBLIC :: SRCGRPFLAG ! use source grouping
        LOGICAL, PUBLIC :: VARFLAG ! use variable grid definition
        LOGICAL, PUBLIC :: LMETCHK ! compare met information in headers
        LOGICAL, PUBLIC :: LMKTPON ! true: use market penetration
        LOGICAL, PUBLIC :: LGRDOUT ! output gridded I/O API NetCDF file(s)
        LOGICAL, PUBLIC :: LREPSTA ! output state-total emissions in report(s)
        LOGICAL, PUBLIC :: LREPCNY ! output county-total emissions in report(s)
        LOGICAL, PUBLIC :: LREPSCC ! output SCC-total emissions in report(s) (Movesmrg only)
        LOGICAL, PUBLIC :: LREPSRC ! output source-level emissions in report(s) (Movesmrg only)
        LOGICAL, PUBLIC :: LREPANY ! true: LREPSTA = T or LREPCNY = T
        LOGICAL, PUBLIC :: LREPINV ! include inventory emissions in report(s)
        LOGICAL, PUBLIC :: LREPSPC ! include speciated emissions in report(s)
        LOGICAL, PUBLIC :: LREPCTL ! include controlled emissions in report(s)
        LOGICAL, PUBLIC :: LAVEDAY ! use average day emissions from the inven

!.........  FILE NAMES AND UNIT NUMBERS...

!.........  Input file unit numbers
        INTEGER     , PUBLIC :: ASDEV  ! area ASCII inventory input
        INTEGER     , PUBLIC :: CDEV   ! state/county names file
        INTEGER     , PUBLIC :: EDEV   ! elevated/PinG ASCII input file
        INTEGER     , PUBLIC :: GDEV   ! gridding surrogates file
        INTEGER     , PUBLIC :: MSDEV  ! mobile ASCII inventory input
        INTEGER     , PUBLIC :: PDEV  = 0  ! inventory table
        INTEGER     , PUBLIC :: PSDEV  ! point ASCII inventory input
        INTEGER     , PUBLIC :: CFDEV  ! control factor input file
        INTEGER     , PUBLIC :: SGDEV  ! source grouping input file

!.........  Output file unit numbers
        INTEGER     , PUBLIC :: ARDEV  ! area ASCII report output
        INTEGER     , PUBLIC :: BRDEV  ! bioegnic ASCII report output
        INTEGER     , PUBLIC :: EVDEV  ! elevated source ASCII output
        INTEGER     , PUBLIC :: MRDEV  ! mobile ASCII report output
        INTEGER     , PUBLIC :: PRDEV  ! point ASCII report output
        INTEGER     , PUBLIC :: TRDEV  ! total ASCII report output

!.........  Logical file names
        CHARACTER(16), PUBLIC :: AENAME ! area inventory input
        CHARACTER(16), PUBLIC :: ATNAME( 7 ) ! area temporal input
        CHARACTER(16), PUBLIC :: AGNAME ! area gridding matrix input
        CHARACTER(16), PUBLIC :: ASNAME ! area speciation matrix input 
        CHARACTER(16), PUBLIC :: ARNAME ! area reactivity control matrix input
        CHARACTER(16), PUBLIC :: AUNAME ! area multiplicative cntl matrix input
        CHARACTER(16), PUBLIC :: AONAME ! area output gridded file
        CHARACTER(16), PUBLIC :: AREPNAME ! area report file name

        CHARACTER(16), PUBLIC :: BONAME ! biogenic output for units conversion
        CHARACTER(16), PUBLIC :: BTNAME ! biogenic gridded, temporal input
        CHARACTER(16), PUBLIC :: BREPNAME ! biogenic report file name

        CHARACTER(16), PUBLIC :: MENAME ! mobile inventory input
        CHARACTER(16), PUBLIC :: MTNAME( 7 ) ! mobile temporal input
        CHARACTER(16), PUBLIC :: MGNAME ! mobile gridding matrix input
        CHARACTER(16), PUBLIC :: MSNAME ! mobile speciation matrix input
        CHARACTER(16), PUBLIC :: MRNAME ! mobile reactivity control matrix input
        CHARACTER(16), PUBLIC :: MUNAME ! mobile multiplicative cntl matrix input
        CHARACTER(16), PUBLIC :: MONAME ! mobile output gridded file
        CHARACTER(16), PUBLIC :: MREPNAME ! mobile report file name
        CHARACTER(16), PUBLIC :: MTMPNAME ! Movesmrg mtmp intmediate output file name

        CHARACTER(16), PUBLIC :: PENAME ! point inventory input
        CHARACTER(16), PUBLIC :: PTNAME( 7 ) ! point temporal input
        CHARACTER(16), PUBLIC :: PGNAME ! point gridding matrix input
        CHARACTER(16), PUBLIC :: PSNAME ! point speciation matrix input 
        CHARACTER(16), PUBLIC :: PLNAME ! point layer fractions file 
        CHARACTER(16), PUBLIC :: PRNAME ! point reactivity control matrix input
        CHARACTER(16), PUBLIC :: PUNAME ! point multiplicative cntl matrix input
        CHARACTER(16), PUBLIC :: PONAME ! point output gridded file
        CHARACTER(16), PUBLIC :: PVNAME ! point elevated stack groups file
        CHARACTER(16), PUBLIC :: PHNAME ! point hourly explicit plume info file
        CHARACTER(16), PUBLIC :: PINGNAME ! plume-in-grid emissions file name
        CHARACTER(16), PUBLIC :: INLINENAME ! in-line pt source emissions file name
        CHARACTER(16), PUBLIC :: PELVNAME ! elevated ASCII emissions file name
        CHARACTER(16), PUBLIC :: PREPNAME ! point report file name

        CHARACTER(16), PUBLIC :: TONAME ! total output gridded file
        CHARACTER(16), PUBLIC :: TREPNAME ! total report file name

!.........  DIMENSIONS AND SOURCE CATEGORY PROPERTIES...
!.........  Number of sources
        INTEGER, PUBLIC :: NASRC = 0   ! number of area sources
        INTEGER, PUBLIC :: NMSRC = 0   ! number of mobile sources
        INTEGER, PUBLIC :: NPSRC = 0   ! number of point sources
        INTEGER, PUBLIC :: JSTACK = 0  ! pt src pos'n of stack in src chars

!.........  Number of emissions layers
        INTEGER, PUBLIC :: EMLAYS = 1

!.........  Number of inventory pollutants
        INTEGER, PUBLIC :: ANIPOL = 0  ! area
        INTEGER, PUBLIC :: BNIPOL = 0  ! biogenic
        INTEGER, PUBLIC :: MNIPOL = 0  ! mobile
        INTEGER, PUBLIC :: PNIPOL = 0  ! point
        INTEGER, PUBLIC :: NIPOL  = 0  ! global

!.........  Number of inventory activities
        INTEGER, PUBLIC :: MNIACT = 0  ! mobile

!.........  Number of inventory pollutants and activities
        INTEGER, PUBLIC :: MNIPPA = 0  ! mobile
        INTEGER, PUBLIC :: NIPPA  = 0  ! global

!.........  Number of pol-to-species combinations
        INTEGER, PUBLIC :: ANSMATV = 0  ! area
        INTEGER, PUBLIC :: BNSMATV = 0  ! biogenics (actually, no. of species)
        INTEGER, PUBLIC :: MNSMATV = 0  ! mobile
        INTEGER, PUBLIC :: PNSMATV = 0  ! point
        INTEGER, PUBLIC :: NSMATV  = 0  ! global

!.........  Number of multiplicative array variables
        INTEGER, PUBLIC :: ANUMATV = 0  ! area
        INTEGER, PUBLIC :: MNUMATV = 0  ! mobile
        INTEGER, PUBLIC :: PNUMATV = 0  ! point

!.........  Number of reactivity array variables
        INTEGER, PUBLIC :: ANRMATV = 0  ! area
        INTEGER, PUBLIC :: MNRMATV = 0  ! mobile
        INTEGER, PUBLIC :: PNRMATV = 0  ! point

!.........  Number of species 
        INTEGER, PUBLIC :: ANMSPC = 0   ! area
        INTEGER, PUBLIC :: BNMSPC = 0   ! biogenics
        INTEGER, PUBLIC :: MNMSPC = 0   ! mobile
        INTEGER, PUBLIC :: PNMSPC = 0   ! point
        INTEGER, PUBLIC :: NMSPC  = 0   ! all

!.........  Number of gridding matrix entries
        INTEGER, PUBLIC :: ANGMAT = 0   ! area
        INTEGER, PUBLIC :: MNGMAT = 0   ! mobile

!.........  Reactivtity matrix sizes
        INTEGER, PUBLIC :: ANSREAC = 0   ! area  , no. src w/ reac controls
        INTEGER, PUBLIC :: MNSREAC = 0   ! mobile, no. src w/ reac controls
        INTEGER, PUBLIC :: PNSREAC = 0   ! point , no. src w/ reac controls
        INTEGER, PUBLIC :: ARNMSPC = 0   ! area  , no. reac species
        INTEGER, PUBLIC :: MRNMSPC = 0   ! mobile, no. reac species
        INTEGER, PUBLIC :: PRNMSPC = 0   ! point , no. reac species
 
!.........  Pollutant names from inventory/temporal file
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: AEINAM( : ) ! area
        CHARACTER(IOVLEN3)             , PUBLIC :: BEINAM( 2 ) ! biogenic
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: MEINAM( : ) ! mobile
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: PEINAM( : ) ! point
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: EINAM ( : ) ! all

!.........  Associated pollutant variable names from inventory file, O=other
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: AONAMES( : ) ! area
        CHARACTER(IOVLEN3),              PUBLIC :: BONAMES( 2 ) ! biogenic
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: MONAMES( : ) ! mobile
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: PONAMES( : ) ! point
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: TONAMES( : ) ! all

!.........  Associated pollutant variable units from inventory file
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: AOUNITS( : ) ! area
        CHARACTER(IOULEN3),              PUBLIC :: BOUNITS( 2 ) ! biogenic
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: MOUNITS( : ) ! mobile
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: POUNITS( : ) ! point
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: TOUNITS( : ) ! all

!.........  Pollutant and activity names from inventory/temporal file
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: MEANAM( : ) ! mobile
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: EANAM ( : ) ! all

!.........  Pollutant-to-species names from speciation file
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: ASVDESC( : ) ! area
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: BSVDESC( : ) ! biogenic
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: MSVDESC( : ) ! mobile
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: PSVDESC( : ) ! point
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: TSVDESC( : ) ! all

!.........  Pollutant-to-species units from speciation file
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: ASVUNIT( : ) ! area
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: BSVUNIT( : ) ! biogenic
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: MSVUNIT( : ) ! mobile
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: PSVUNIT( : ) ! point
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: TSVUNIT( : ) ! all

!.........  Pollutant-to-species names from reactivity file
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: ARVDESC( : ) ! area
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: MRVDESC( : ) ! mobile
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: PRVDESC( : ) ! point

!.........  Pollutant names from multiplicative control matrix
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: AUVNAMS( : ) ! area
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: MUVNAMS( : ) ! mobile
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: PUVNAMS( : ) ! point

!.........  Species names
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: AEMNAM  ( : ) ! area
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: BEMNAM  ( : ) ! mobile
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: MEMNAM  ( : ) ! biogenic
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: PEMNAM  ( : ) ! point
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC ::  EMNAM  ( : ) ! all

!.........  Index between all model species names and all inventory pollutants
        INTEGER, ALLOCATABLE, PUBLIC :: EMIDX( : )

!.........  Number of unique units needed
!           If using speciation, NUNITS = NMPSC; otherwise, NUNITS = NIPPA
        INTEGER, PUBLIC :: NUNITS = 0

!.........  Map file pollutants list and physical file names
        INTEGER,                         PUBLIC :: ANMAP = 0
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: AMAPNAM( : )
        CHARACTER(PHYLEN3), ALLOCATABLE, PUBLIC :: AMAPFIL( : )
        INTEGER,                         PUBLIC :: MNMAP = 0
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: MMAPNAM( : )
        CHARACTER(PHYLEN3), ALLOCATABLE, PUBLIC :: MMAPFIL( : )
        INTEGER,                         PUBLIC :: PNMAP = 0
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: PMAPNAM( : )
        CHARACTER(PHYLEN3), ALLOCATABLE, PUBLIC :: PMAPFIL( : )

!.........  NAMES AND INDICES FOR GROUP STRUCTURE

!.........  Grouped variable arrays
        INTEGER, ALLOCATABLE, PUBLIC :: VGRPCNT( : )  ! count per group
        INTEGER, ALLOCATABLE, PUBLIC :: IDVGP( : )    ! group ID number
        CHARACTER(IODLEN3), ALLOCATABLE, PUBLIC :: GVNAMES( :,: ) ! var name
        INTEGER, ALLOCATABLE, PUBLIC :: SIINDEX( :,: ) ! index to EANAM
        INTEGER, ALLOCATABLE, PUBLIC :: SPINDEX( :,: ) ! index to EMNAM
        LOGICAL, ALLOCATABLE, PUBLIC :: GVLOUT ( :,: ) ! output points

!.........  Group indices for inventory emissions
        INTEGER, ALLOCATABLE, PUBLIC :: A_EXIST( :,: )   ! area
        INTEGER, ALLOCATABLE, PUBLIC :: M_EXIST( :,: )   ! mobile
        INTEGER, ALLOCATABLE, PUBLIC :: P_EXIST( :,: )   ! point

!.........  Group indices for multiplicative control matrices
        INTEGER, ALLOCATABLE, PUBLIC :: AU_EXIST( :,: )   ! area
        INTEGER, ALLOCATABLE, PUBLIC :: MU_EXIST( :,: )   ! mobile
        INTEGER, ALLOCATABLE, PUBLIC :: PU_EXIST( :,: )   ! point

!.........  Group indices for reactivity control matrices
        INTEGER, ALLOCATABLE, PUBLIC :: AR_EXIST( :,: )   ! area
        INTEGER, ALLOCATABLE, PUBLIC :: MR_EXIST( :,: )   ! mobile
        INTEGER, ALLOCATABLE, PUBLIC :: PR_EXIST( :,: )   ! point

!.........  Group indices for speciation matrices
        INTEGER, ALLOCATABLE, PUBLIC :: AS_EXIST( :,: )   ! area
        INTEGER, ALLOCATABLE, PUBLIC :: BS_EXIST( :,: )   ! biogenic 
        INTEGER, ALLOCATABLE, PUBLIC :: MS_EXIST( :,: )   ! mobile
        INTEGER, ALLOCATABLE, PUBLIC :: PS_EXIST( :,: )   ! point

!.........  EPISODE SETTINGS...

!.........  Episode information
        INTEGER, PUBLIC :: SDATE  = 0     ! Julian start date of episode
        INTEGER, PUBLIC :: STIME  = 0     ! start time of episode
        INTEGER, PUBLIC :: EDATE  = 0     ! Julian end date
        INTEGER, PUBLIC :: ETIME  = 0     ! end time of episode
        INTEGER, PUBLIC :: NSTEPS = 1     ! number of time loop iterations
        INTEGER, PUBLIC :: TSTEP  = 0     ! hourly time steps
        INTEGER, PUBLIC :: TZONE  = -1    ! time zone
        INTEGER, PUBLIC :: BYEAR  = 0     ! base inventory year
        INTEGER, PUBLIC :: PYEAR  = 0     ! projected inventory year
        INTEGER, PUBLIC :: PVSDATE= 0     ! Julian start date in STACK_GROUPS
        INTEGER, PUBLIC :: PVSTIME= 0     ! start time in STACK_GROUPS

!.........  Dates for by-day hourly emissions input
        INTEGER, PUBLIC :: ASDATE( 7 )    ! Julian start dates of each atmp file
        INTEGER, PUBLIC :: MSDATE( 7 )    ! Julian start dates of each mtmp file
        INTEGER, PUBLIC :: PSDATE( 7 )    ! Julian start dates of each ptmp file

!.........  Units conversions information
        REAL             , PUBLIC :: BIOGFAC     ! conv fac for gridded bio
        REAL             , PUBLIC :: BIOTFAC     ! conv fac for bio totals
        REAL, ALLOCATABLE, PUBLIC :: GRDFAC( : ) ! for spc/pol/act grid outputs
        REAL, ALLOCATABLE, PUBLIC :: TOTFAC( : ) ! for spc/pol/act st/co outputs

        CHARACTER(IOULEN3), PUBLIC :: BIOUNIT = ' '  ! biogenic input units
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: GRDUNIT( : ) ! gridded output units
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: SPCUNIT( : ) ! speciation units
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: TOTUNIT( : ) ! st/co total units

!.........  INPUT ARRAYS ...
!.........  Country/State/County codes
        INTEGER,            ALLOCATABLE, PUBLIC :: MCFIP( : ) ! mobile code idex for Movesmrg
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: AIFIP( : ) ! area codes, dim: nasrc
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: MIFIP( : ) ! mobile codes, dim: nmsrc
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: PIFIP( : ) ! point codes, dim: npsrc

!.........  Emissions (either inventory or hourly)
        REAL   , ALLOCATABLE, PUBLIC :: AEMSRC( :,: ) ! area: nasrc, mxpolpgp
        REAL   , ALLOCATABLE, PUBLIC :: MEMSRC( :,: ) ! mobile: nmsrc, mxpolpgp
        REAL   , ALLOCATABLE, PUBLIC :: PEMSRC( :,: ) ! point: npsrc, mxpolpgp

!.........  Emissions (for inventory only)
        REAL   , ALLOCATABLE, PUBLIC :: AEISRC( :,: ) ! area: nasrc, mxpolpgp
        REAL   , ALLOCATABLE, PUBLIC :: MEISRC( :,: ) ! mobile: nmsrc, mxpolpgp
        REAL   , ALLOCATABLE, PUBLIC :: PEISRC( :,: ) ! point: npsrc, mxpolpgp

!.........  Biogenic gridded emissions
        REAL   , ALLOCATABLE, PUBLIC :: BEMGRD( : )   ! dim: ngrid

!.........  Layer fractions
        REAL   , ALLOCATABLE, PUBLIC :: LFRAC( :,: )  ! dim: npsrc, emlays

!.........  Elevated adjustments buffer for Mrgmult
        REAL   , ALLOCATABLE, PUBLIC :: ELEVADJ( : )! max(nasrc,nmsrc,npsrc)

!.........  Gridding matrices, dim depends on src category
        REAL   , ALLOCATABLE, PUBLIC :: AGMATX( : ) ! contiguous area gridding
        REAL   , ALLOCATABLE, PUBLIC :: MGMATX( : ) ! contiguous mobile gridding
        REAL   , ALLOCATABLE, PUBLIC :: PGMATX( : ) ! contiguous point gridding

!.........  Speciation matrices, dim n*src, 
        REAL   , ALLOCATABLE, PUBLIC :: ASMATX( :,: ) ! area
        REAL   , ALLOCATABLE, PUBLIC :: MSMATX( :,: ) ! mobile
        REAL   , ALLOCATABLE, PUBLIC :: PSMATX( :,: ) ! point

!.........  Reactivity matrix temporary source-based arrays, 1=emis, 2=mkt pen
        REAL   , ALLOCATABLE, PUBLIC :: ARINFO( :,: ) ! area
        REAL   , ALLOCATABLE, PUBLIC :: MRINFO( :,: ) ! mobile
        REAL   , ALLOCATABLE, PUBLIC :: PRINFO( :,: ) ! point

!.........  OUTPUT ARRAYS ...

!.........  Gridded emissions (or for mobile, VMT)
        REAL   , ALLOCATABLE, PUBLIC :: AEMGRD( : )   ! area  , dim: ngrid
        REAL   , ALLOCATABLE, PUBLIC :: MEMGRD( : )   ! mobile, dim: ngrid
        REAL   , ALLOCATABLE, PUBLIC :: PEMGRD( :,: ) ! point , dim: ngrid
        REAL   , ALLOCATABLE, PUBLIC :: TEMGRD( :,: ) ! all,3d, dim: ngrid,nlay

!.........  In the following lists, ndim can either by nipol (number of 
!           inventory pollutants) or nmspc (number of model species)

!.........  State total inventory of speciated emissions (or VMT for mobile)
        REAL   , ALLOCATABLE, PUBLIC :: AEBSTA( :,: ) ! area  , dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: BEBSTA( :,: ) ! biogen, dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MEBSTA( :,: ) ! mobile, dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PEBSTA( :,: ) ! point , dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: TEBSTA( :,: ) ! all   , dim nsta, ndim

!.........  State total multiplicative-controlled emissions
        REAL   , ALLOCATABLE, PUBLIC :: AEUSTA( :,: ) ! area  , dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MEUSTA( :,: ) ! mobile, dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PEUSTA( :,: ) ! point , dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: TEUSTA( :,: ) ! all   , dim nsta, ndim

!.........  State total reactivity-controlled emissions
        REAL   , ALLOCATABLE, PUBLIC :: AERSTA( :,: ) ! area  , dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MERSTA( :,: ) ! mobile, dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PERSTA( :,: ) ! point , dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: TERSTA( :,: ) ! all   , dim nsta, ndim

!.........  State total all-controlled emissions
        REAL   , ALLOCATABLE, PUBLIC :: AECSTA( :,: ) ! area  , dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MECSTA( :,: ) ! mobile, dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PECSTA( :,: ) ! point , dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: TECSTA( :,: ) ! all   , dim nsta, ndim

!.........  County total inventory of speciated emissions (or VMT for mobile)
        REAL   , ALLOCATABLE, PUBLIC :: AEBCNY( :,: ) ! area  , dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: BEBCNY( :,: ) ! biogen, dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MEBCNY( :,: ) ! mobile, dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PEBCNY( :,: ) ! point , dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: TEBCNY( :,: ) ! all   , dim ncny, ndim

!.........  County total multiplicative-controlled emissions
        REAL   , ALLOCATABLE, PUBLIC :: AEUCNY( :,: ) ! area  , dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MEUCNY( :,: ) ! mobile, dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PEUCNY( :,: ) ! point , dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: TEUCNY( :,: ) ! all   , dim ncny, ndim

!.........  County total reactivity-controlled emissions
        REAL   , ALLOCATABLE, PUBLIC :: AERCNY( :,: ) ! area  , dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MERCNY( :,: ) ! mobile, dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PERCNY( :,: ) ! point , dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: TERCNY( :,: ) ! all   , dim ncny, ndim

!.........  County total all-controlled emissions
        REAL   , ALLOCATABLE, PUBLIC :: AECCNY( :,: ) ! area  , dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MECCNY( :,: ) ! mobile, dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PECCNY( :,: ) ! point , dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: TECCNY( :,: ) ! all   , dim ncny, ndim

!.........  Source and SCC total speciated emissions
        REAL   , ALLOCATABLE, PUBLIC :: MEBSCC( :,: ) ! mobile, dim nscc, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MEBSUM( :,: ) ! mobile, dim nscc, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MEBSRC( :,:,: ) ! mobile, dim nsrc, ndim by hour
        REAL   , ALLOCATABLE, PUBLIC :: MEBSTC( :,:,: ) ! mobile, dim nsta, nscc, ndim

        REAL   , ALLOCATABLE, PUBLIC :: PEBSCC( :,: ) ! point, dim nscc, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PEBSUM( :,: ) ! point, dim nscc, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PEBSRC( :,:,: ) ! point, dim nsrc, ndim by hour
        REAL   , ALLOCATABLE, PUBLIC :: PEBSTC( :,:,: ) !  point, dim nsta, nscc, ndim

!.........  Source apportionment storage
        INTEGER,              PUBLIC :: NSRCGRP       ! total number of source group criteria
        INTEGER,              PUBLIC :: NGRPS         ! total number of groups
        INTEGER, ALLOCATABLE, PUBLIC :: IUGRPNUM( : ) ! list of unique source group numbers, dim: ngrps
        INTEGER, ALLOCATABLE, PUBLIC :: IUGRPIDX( : ) ! idx for each unique source group, dim: max group num
        INTEGER, ALLOCATABLE, PUBLIC :: IGRPNUM( : )  ! list of source group numbers, dim: nsrcgrp
        INTEGER, ALLOCATABLE, PUBLIC :: ISRCGRP( : )  ! source group idx for each source
        REAL,    ALLOCATABLE, PUBLIC :: EMGGRD( :,: ) ! emissions by grid cell and source group
        REAL,    ALLOCATABLE, PUBLIC :: EMGGRDSPC( :,:,: )    ! EMGGRD + species
        REAL,    ALLOCATABLE, PUBLIC :: EMGGRDSPCT( :,:,:,: ) ! EMGGRDSPC + time step
        INTEGER,              PUBLIC :: NSGOUTPUT     ! number of output records
        INTEGER, ALLOCATABLE, PUBLIC :: GRPCNT( :,: ) ! num srcs matching grid cell and group
        CHARACTER(16),        PUBLIC :: SRCGRPNAME    ! source group output file name (stack groups)
        CHARACTER(16),        PUBLIC :: SGINLNNAME    ! in-line src grp emissions file name
        CHARACTER(16), ALLOCATABLE, PUBLIC :: SUBOUTNAME( : )    ! sub-sector emissions output 

        END MODULE MODMERGE
