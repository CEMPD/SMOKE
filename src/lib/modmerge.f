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
! COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
! All Rights Reserved
!
! See file COPYRIGHT for conditions of use.
!
! Environmental Programs Group
! MCNC--North Carolina Supercomputing Center
! P.O. Box 12889
! Research Triangle Park, NC  27709-2889
!
! env_progs@mcnc.org
!
! Pathname: $Source$
! Last updated: $Date$ 
!
!****************************************************************************

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

!.........  Flags to indicate additive controls being used
        LOGICAL, PUBLIC :: AAFLAG
        LOGICAL, PUBLIC :: MAFLAG
        LOGICAL, PUBLIC :: PAFLAG
        LOGICAL, PUBLIC :: TAFLAG  ! total

!.........  Flags to indicate reactivity controls being used
        LOGICAL, PUBLIC :: ARFLAG
        LOGICAL, PUBLIC :: MRFLAG
        LOGICAL, PUBLIC :: PRFLAG
        LOGICAL, PUBLIC :: TRFLAG  ! total

!.........  Flags to indicate inventory emissions are projected or not
        LOGICAL, PUBLIC :: APRJFLAG = .FALSE.
        LOGICAL, PUBLIC :: MPRJFLAG = .FALSE.
        LOGICAL, PUBLIC :: PPRJFLAG = .FALSE.

!.........  Flags for controlling on-off features
        LOGICAL, PUBLIC :: TFLAG   ! use temporalized (hourly) emissions
        LOGICAL, PUBLIC :: SFLAG   ! use speciation
        LOGICAL, PUBLIC :: LFLAG   ! use layer fractions
        LOGICAL, PUBLIC :: VFLAG   ! process VMT (not emission)
        LOGICAL, PUBLIC :: LMETCHK ! compare met information in headers
        LOGICAL, PUBLIC :: LMKTPON ! true: use market penetration
        LOGICAL, PUBLIC :: LGRDOUT ! output gridded I/O API NetCDF file(s)
        LOGICAL, PUBLIC :: LREPSTA ! output state-total emissions in report(s)
        LOGICAL, PUBLIC :: LREPCNY ! output county-total emissions in report(s)
        LOGICAL, PUBLIC :: LREPINV ! include inventory emissions in report(s)
        LOGICAL, PUBLIC :: LREPSPC ! include speciated emissions in report(s)
        LOGICAL, PUBLIC :: LREPCTL ! include controlled emissions in report(s)

!.........  FILE NAMES AND UNIT NUMBERS...

!.........  Unit numbers
        INTEGER     , PUBLIC :: ASDEV  ! area ASCII inventory input
        INTEGER     , PUBLIC :: MSDEV  ! mobile ASCII inventory input
        INTEGER     , PUBLIC :: PSDEV  ! point ASCII inventory input
        INTEGER     , PUBLIC :: PDEV   ! inventory pollutants list
        INTEGER     , PUBLIC :: ARDEV  ! area ASCII report output
        INTEGER     , PUBLIC :: MRDEV  ! mobile ASCII report output
        INTEGER     , PUBLIC :: PRDEV  ! point ASCII report output
        INTEGER     , PUBLIC :: TRDEV  ! total ASCII report output

!.........  Logical file names
        CHARACTER*16, PUBLIC :: AENAME ! area inventory input
        CHARACTER*16, PUBLIC :: ATNAME ! area temporal input
        CHARACTER*16, PUBLIC :: AGNAME ! area gridding matrix input
        CHARACTER*16, PUBLIC :: ASNAME ! area speciation matrix input 
        CHARACTER*16, PUBLIC :: AANAME ! area additive control matrix input
        CHARACTER*16, PUBLIC :: ARNAME ! area reactivity control matrix input
        CHARACTER*16, PUBLIC :: AUNAME ! area multiplicative cntl matrix input
        CHARACTER*16, PUBLIC :: AONAME ! area output gridded file
        CHARACTER*16, PUBLIC :: AREPNAME ! area report file name

        CHARACTER*16, PUBLIC :: BTNAME ! biogenic gridded, temporal input
        CHARACTER*16, PUBLIC :: BREPNAME ! biogenic report file name

        CHARACTER*16, PUBLIC :: MENAME ! mobile inventory input
        CHARACTER*16, PUBLIC :: MTNAME ! mobile temporal input
        CHARACTER*16, PUBLIC :: MGNAME ! mobile gridding matrix input
        CHARACTER*16, PUBLIC :: MSNAME ! mobile speciation matrix input
        CHARACTER*16, PUBLIC :: MANAME ! mobile additive control matrix input
        CHARACTER*16, PUBLIC :: MRNAME ! mobile reactivity control matrix input
        CHARACTER*16, PUBLIC :: MUNAME ! mobile multiplicative cntl matrix input
        CHARACTER*16, PUBLIC :: MONAME ! mobile output gridded file
        CHARACTER*16, PUBLIC :: MREPNAME ! mobile report file name

        CHARACTER*16, PUBLIC :: PENAME ! point inventory input
        CHARACTER*16, PUBLIC :: PTNAME ! point temporal input
        CHARACTER*16, PUBLIC :: PGNAME ! point gridding matrix input
        CHARACTER*16, PUBLIC :: PSNAME ! point speciation matrix input 
        CHARACTER*16, PUBLIC :: PLNAME ! point layer fractions file 
        CHARACTER*16, PUBLIC :: PANAME ! point additive control matrix input
        CHARACTER*16, PUBLIC :: PRNAME ! point reactivity control matrix input
        CHARACTER*16, PUBLIC :: PUNAME ! point multiplicative cntl matrix input
        CHARACTER*16, PUBLIC :: PONAME ! point output gridded file
        CHARACTER*16, PUBLIC :: PREPNAME ! point report file name

        CHARACTER*16, PUBLIC :: TONAME ! total output gridded file

!.........  DIMENSIONS AND SOURCE CATEGORY PROPERTIES...
!.........  Number of sources
        INTEGER, PUBLIC :: NASRC = 0   ! number of area sources
        INTEGER, PUBLIC :: NMSRC = 0   ! number of mobile sources
        INTEGER, PUBLIC :: NPSRC = 0   ! number of point sources

!.........  Number of inventory pollutants
        INTEGER, PUBLIC :: ANIPOL = 0  ! area
        INTEGER, PUBLIC :: MNIPOL = 0  ! mobile
        INTEGER, PUBLIC :: PNIPOL = 0  ! point
        INTEGER, PUBLIC :: NIPOL  = 0  ! global

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

!.........  Number of additive array variables
        INTEGER, PUBLIC :: ANAMATV = 0  ! area
        INTEGER, PUBLIC :: MNAMATV = 0  ! mobile
        INTEGER, PUBLIC :: PNAMATV = 0  ! point

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

!.........  State/county information
        INTEGER, PUBLIC :: NSTA = 0     ! number of country/state combinations
        INTEGER, PUBLIC :: NCNY = 0     ! number of counties

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
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: AEINAM( : ) ! area
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: MEINAM( : ) ! mobile
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: PEINAM( : ) ! point
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: EINAM ( : ) ! all

!.........  Pollutant-to-species names from speciation file
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE, PUBLIC :: ASVDESC( : ) ! area
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE, PUBLIC :: BSVDESC( : ) ! biogenic
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE, PUBLIC :: MSVDESC( : ) ! mobile
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE, PUBLIC :: PSVDESC( : ) ! point
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE, PUBLIC :: TSVDESC( : ) ! all

!.........  Pollutant-to-species names from reactivity file
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE, PUBLIC :: ARVDESC( : ) ! area
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE, PUBLIC :: MRVDESC( : ) ! mobile
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE, PUBLIC :: PRVDESC( : ) ! point

!.........  Pollutant names from multiplicative control matrix
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: AUVNAMS( : ) ! area
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: MUVNAMS( : ) ! mobile
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: PUVNAMS( : ) ! point

!.........  Pollutant names from additive control matrix
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: AAVNAMS( : ) ! area
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: MAVNAMS( : ) ! mobile
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: PAVNAMS( : ) ! point

!.........  Species names
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: AEMNAM  ( : ) ! area
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: BEMNAM  ( : ) ! mobile
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: MEMNAM  ( : ) ! biogenic
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: PEMNAM  ( : ) ! point
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC ::  EMNAM  ( : ) ! all

!.........  NAMES AND INDICES FOR GROUP STRUCTURE

!.........  Grouped pollutant arrays
        INTEGER, ALLOCATABLE, PUBLIC :: PLCNT( : )  ! count per group
        INTEGER, ALLOCATABLE, PUBLIC :: IDPGP( : )  ! group ID number
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: PLNAMES( :,: ) ! pol name

!.........  Grouped pol-to-species
        INTEGER, ALLOCATABLE, PUBLIC :: NSMPPG( :,: )    ! No. spcs per pol/grp 
        INTEGER, ALLOCATABLE, PUBLIC :: SMINDEX( :,:,: ) ! index to TSVDESC
        INTEGER, ALLOCATABLE, PUBLIC :: SPINDEX( :,:,: ) ! index to EMNAM

!.........  Group indices for inventory emissions
        INTEGER, ALLOCATABLE, PUBLIC :: A_EXIST( :,: )   ! area
        INTEGER, ALLOCATABLE, PUBLIC :: M_EXIST( :,: )   ! mobile
        INTEGER, ALLOCATABLE, PUBLIC :: P_EXIST( :,: )   ! point

!.........  Group indices for multiplicative control matrices
        INTEGER, ALLOCATABLE, PUBLIC :: AU_EXIST( :,: )   ! area
        INTEGER, ALLOCATABLE, PUBLIC :: MU_EXIST( :,: )   ! mobile
        INTEGER, ALLOCATABLE, PUBLIC :: PU_EXIST( :,: )   ! point

!.........  Group indices for additive control matrices
        INTEGER, ALLOCATABLE, PUBLIC :: AA_EXIST( :,: )   ! area
        INTEGER, ALLOCATABLE, PUBLIC :: MA_EXIST( :,: )   ! mobile
        INTEGER, ALLOCATABLE, PUBLIC :: PA_EXIST( :,: )   ! point

!.........  Group indices for reactivity control matrices
        INTEGER, ALLOCATABLE, PUBLIC :: AR_EXIST( :,:,: )   ! area
        INTEGER, ALLOCATABLE, PUBLIC :: MR_EXIST( :,:,: )   ! mobile
        INTEGER, ALLOCATABLE, PUBLIC :: PR_EXIST( :,:,: )   ! point

!.........  Group indices for speciation matrices
        INTEGER, ALLOCATABLE, PUBLIC :: AS_EXIST( :,:,: )   ! area
        INTEGER, ALLOCATABLE, PUBLIC :: BS_EXIST( :,:,: )   ! biogenic 
        INTEGER, ALLOCATABLE, PUBLIC :: MS_EXIST( :,:,: )   ! mobile
        INTEGER, ALLOCATABLE, PUBLIC :: PS_EXIST( :,:,: )   ! point

!.........  GRID AND EPISODE SETTINGS...
!.........  Grid information
        CHARACTER(LEN=IOVLEN3), PUBLIC :: GRDNM = ' '  ! grid name
        CHARACTER(LEN=IOVLEN3), PUBLIC :: COORD = ' '  ! coord system name
        REAL   , PUBLIC :: GDTYP = -1     ! i/o api grid type code
        REAL   , PUBLIC :: P_ALP = 0.D0   ! projection alpha
        REAL   , PUBLIC :: P_BET = 0.D0   ! projection beta
        REAL   , PUBLIC :: P_GAM = 0.D0   ! projection gamma
        REAL   , PUBLIC :: XCENT = 0.D0   ! x-center of projection
        REAL   , PUBLIC :: YCENT = 0.D0   ! y-center of projection
        REAL   , PUBLIC :: XORIG = 0.D0   ! x-origin of grid
        REAL   , PUBLIC :: YORIG = 0.D0   ! y-origin of grid
        REAL   , PUBLIC :: XCELL = 0.D0   ! x-dim of cells
        REAL   , PUBLIC :: YCELL = 0.D0   ! y-dim of cells
        INTEGER, PUBLIC :: NCOLS = 0      ! number of columns in grid
        INTEGER, PUBLIC :: NROWS = 0      ! number of rows in grid
        INTEGER, PUBLIC :: NGRID = 0      ! number of cells in grid

!.........  Vertical structure information
        INTEGER, PUBLIC :: EMLAYS = 1      ! number of emissions layers
        INTEGER, PUBLIC :: VGTYP  = -1     ! type of vertical coordinates
        REAL   , PUBLIC :: VGTOP  = 0.0    ! model-top, for sigma coord types
        REAL   , ALLOCATABLE, PUBLIC :: VGLVS( : ) ! vertical coordinate values

!.........  Episode information
        INTEGER, PUBLIC :: SDATE  = 0     ! Julian start date of episode
        INTEGER, PUBLIC :: STIME  = 0     ! start time of episode
        INTEGER, PUBLIC :: NSTEPS = 1     ! number of time loop iterations
        INTEGER, PUBLIC :: TSTEP  = 10000 ! hourly time steps
        INTEGER, PUBLIC :: TZONE  = -1    ! time zone
        INTEGER, PUBLIC :: BYEAR  = 0     ! base inventory year
        INTEGER, PUBLIC :: PYEAR  = 0     ! projected inventory year

!.........  Units
        CHARACTER(LEN=IOVLEN3), PUBLIC :: INVUNIT = ' ' ! emis units w/time
        CHARACTER(LEN=IOVLEN3), PUBLIC :: SPCUNIT = ' ' ! speciation units

!.........  INPUT ARRAYS ...
!.........  Country/State/County codes
        INTEGER, ALLOCATABLE, PUBLIC :: AIFIP( : )    ! area codes, dim: nasrc
        INTEGER, ALLOCATABLE, PUBLIC :: MIFIP( : )    ! mobile codes, dim: nmsrc
        INTEGER, ALLOCATABLE, PUBLIC :: PIFIP( : )    ! point codes, dim: npsrc

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
        REAL   , ALLOCATABLE, PUBLIC :: PEMGRD( : )   ! point , dim: ngrid
        REAL   , ALLOCATABLE, PUBLIC :: TEMGRD( :,: ) ! all,3d, dim: ngrid,nlay

!.........  State by source and county by source pointers for summing emissions
!           by state and county. They point to state and county sum arrays
        INTEGER, ALLOCATABLE, PUBLIC :: AIXSTA( : )  ! area, dim NASRC
        INTEGER, ALLOCATABLE, PUBLIC :: BIXSTA( : )  ! biogen, dim NGRID
        INTEGER, ALLOCATABLE, PUBLIC :: MIXSTA( : )  ! mobile, dim NMSRC
        INTEGER, ALLOCATABLE, PUBLIC :: PIXSTA( : )  ! point, dim NPSRC

        INTEGER, ALLOCATABLE, PUBLIC :: AIXCNY( : )  ! area, dim NASRC
        INTEGER, ALLOCATABLE, PUBLIC :: BIXCNY( : )  ! biogen, dim NGRID
        INTEGER, ALLOCATABLE, PUBLIC :: MIXCNY( : )  ! mobile, dim NMSRC
        INTEGER, ALLOCATABLE, PUBLIC :: PIXCNY( : )  ! point, dim NPSRC

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

!.........  State total additive-controlled emissions
        REAL   , ALLOCATABLE, PUBLIC :: AEASTA( :,: ) ! area  , dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MEASTA( :,: ) ! mobile, dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PEASTA( :,: ) ! point , dim nsta, ndim
        REAL   , ALLOCATABLE, PUBLIC :: TEASTA( :,: ) ! all   , dim nsta, ndim

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

!.........  County total additive-controlled emissions
        REAL   , ALLOCATABLE, PUBLIC :: AEACNY( :,: ) ! area  , dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: MEACNY( :,: ) ! mobile, dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: PEACNY( :,: ) ! point , dim ncny, ndim
        REAL   , ALLOCATABLE, PUBLIC :: TEACNY( :,: ) ! all   , dim ncny, ndim

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

        END MODULE MODMERGE
