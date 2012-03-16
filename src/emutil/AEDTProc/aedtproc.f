      PROGRAM AEDTPROC

C***********************************************************************
C  program body starts at line  
C
C  DESCRIPTION:
C       This program processes hourly high resolution AEDT data into 
C       CMAQ-ready hourly emissions for aircraft assessment studies.
C
C  PRECONDITIONS REQUIRED:
C       Preprocessed AEDT airport sources emissions
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 11/2006 by B.H. Baek
C                  
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
C.........  MODULES
        USE MODGRID, ONLY: NCOLS, NROWS, XORIG, YORIG, XOFF, YOFF,
     &                     GDTYP, XCELL, YCELL, XCENT, YCENT, NLAYS,
     &                     P_ALP, P_BET, P_GAM, GRDNM, NGRID, COORD

        USE MODSTCY, ONLY: NCOUNTY, CNTYCOD, CNTYTZNM

        USE MODLISTS, ONLY: ITMSPC, ITCASA, ITVTSA, ITEXPL, ITNAMA,
     &                      NINVTBL, MXIDAT, INVDNAM


        IMPLICIT NONE

C.........  INCLUDES
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'   ! emissions constant parameters
        
C.........  EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL       ENVINT
        INTEGER       FIND1
        INTEGER       FINDC
        INTEGER       GETFLINE
        INTEGER       INDEX1
        INTEGER       JULIAN
        INTEGER       PROMPTFFILE
        INTEGER       STR2INT
        LOGICAL       BLKORCMT
        LOGICAL       CHKREAL
        LOGICAL       INGRID
        LOGICAL       DSCM3GRD

        REAL          STR2REAL, ENVREAL

        CHARACTER(2)  CRLF
        CHARACTER(16) PROMPTMFILE

        EXTERNAL      CRLF,  CHKREAL, FIND1, FINDC, INDEX1, GETFLINE, 
     &                JULIAN, INGRID, STR2INT, STR2REAL, BLKORCMT, 
     &                PROMPTFFILE, PROMPTMFILE, DSCM3GRD, ENVINT, ENVREAL

C.........  LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: 
     &   CVSW = '$Name$'  ! CVS release tag
        INTEGER, PARAMETER :: MXSEG = 33    ! number of segments in line
        INTEGER, PARAMETER :: MXOUT = 11    ! number of output values
        INTEGER, PARAMETER :: NVARS = 43    ! number of unique modeling species

        REAL,    PARAMETER :: FSC = 600.0         ! fuel surfur concentration
        REAL,    PARAMETER :: E   = 2.0           ! surfur mole fractoin

C.........  LOCAL VARIABLES

C.........  Static arrays
        CHARACTER(125)  INP_SEG( 4 )        ! parsed input line
        CHARACTER(32)   SEGMENT( MXSEG )    ! parsed input line

C.........  Allocatable arrays
        REAL,          ALLOCATABLE :: LFRAC ( : )       ! model layer fractions
        REAL,          ALLOCATABLE :: TERRAIN( : )      ! terrain height (m)
        REAL,          ALLOCATABLE :: ZZF   ( :,: )     ! layer's full height (m)
        REAL,          ALLOCATABLE :: SFCHGT( : )       ! surface height (m)
        REAL,          ALLOCATABLE :: VGLVLS( : )       ! gridded mask values to be output
        REAL,          ALLOCATABLE :: CVHAPT( :,: ) 
        REAL,          ALLOCATABLE :: CVHAPP( :,: )
        REAL,          ALLOCATABLE :: CVSPCT( :,: ) 
        REAL,          ALLOCATABLE :: CVSPCP( :,: )
        REAL,          ALLOCATABLE :: TMP3D ( :,:,:,: ) ! tmp emissions

        REAL,          ALLOCATABLE :: APRT_ELEV( : )
        CHARACTER(16), ALLOCATABLE :: APRT_CODE( : )

        CHARACTER(16), ALLOCATABLE :: NPHAPP( : )
        CHARACTER(16), ALLOCATABLE :: NPHAPT( : )
        CHARACTER(16), ALLOCATABLE :: NSSPCP( : )
        CHARACTER(16), ALLOCATABLE :: NSSPCT( : )
        CHARACTER(15), ALLOCATABLE :: ENGINE_TYP( :,: )  ! AEDT Engine types
        CHARACTER(15), ALLOCATABLE :: FLIGHT_ENG( :,: )  ! AEDT FlightIDs and engine types

C.........   Local arrays dimensioned by subroutine arguments
C.........   Note that the NGRID dimension could conceivably be too small if
C            a link winds through the whole domain, but this is a case that
C            is not worth going to extra trouble for since it is not realistic
        INTEGER,       ALLOCATABLE :: ACEL( : )    ! number of cell intersections per src
        REAL,          ALLOCATABLE :: AFAC( : )    ! fraction of link in cell

C.........  File units and logical names
        INTEGER      :: APDEV= 0            ! a list of of airport elevations input file
        INTEGER      :: FDEV = 0            ! a list of of AEDT flight information/emissions files
        INTEGER      :: FLDEV= 0            ! Actual AEDT flight information/emissions files
        INTEGER      :: SDEV = 0            ! a list of of AEDT segment information files
        INTEGER      :: SGDEV= 0            ! Actual AEDT segment information files
        INTEGER      :: TDEV = 0            ! turbine engine HAP list
        INTEGER      :: PDEV = 0            ! piston engine HAP list
        INTEGER      :: TSDEV= 0            ! turbine engine HAP list-new 1098 profile
        INTEGER      :: PSDEV= 0            ! piston engine HAP list-1099 profile
        INTEGER      :: LDEV = 0            ! unit number for log file
        INTEGER      :: RDEV = 0            ! report file

        CHARACTER*16    GNAME               ! grid-point layered met file
        CHARACTER*16    XNAME               ! cross-point layered met file
        CHARACTER*16    ONAME               ! CMAQ-ready gridded emissions file
        
C.........  Other local variables
        INTEGER         C, F, I, IC, J, K, L 
        INTEGER         N, NA, NC, ND, NF, NL, NP, NS, T, V, NV    ! counters and indices
        INTEGER         IOS                 ! i/o status
        INTEGER         IREC                ! line counter
        

        INTEGER      :: N_FLG = 0           ! number of unique flight IDs
        INTEGER      :: N_SEG = 0           ! number of flight segments

        INTEGER         SEGID               ! segment id 
        INTEGER         MODID               ! flight mode id (LTO from departure and to arrival)
        INTEGER         NCEL                ! tmp number of cells
        INTEGER         ORG_CELLID          ! origin cell id
        INTEGER         END_CELLID          ! ending cell id

        INTEGER         NFILES
        INTEGER         NSTEPS

        INTEGER         ROW, COL 

        INTEGER         VGTYP

        INTEGER         HOUR                ! current hour 
        INTEGER         DAY                 ! current day
        INTEGER         MON                 ! current month
        INTEGER         YEAR                ! current modeling year
        INTEGER      :: DYEAR = -1           ! Year diff = Overyear-Modelyear (2005-2006)

        INTEGER         JDATE, TDATE        ! Processing date
        INTEGER         JTIME               ! Processing time

        INTEGER         BDATE, SDATE, EDATE ! Modeling date
        INTEGER         STIME, ETIME        ! Modeling time
        INTEGER         LTOP, LBOT, DDP     ! layer# for top/bottom cells

        INTEGER         STRID, ENDID, INCID ! start/end grid cells and increment

        INTEGER         NAPRT_ELEV          ! A number of airports

        INTEGER         NFAC                ! A number of conversion factors
        INTEGER         NFAC_TURBINE        ! A number of conversion factors for Turbine HAPs
        INTEGER         NFAC_PISTON         ! A number of conversion factors for Piston HAPs

        INTEGER         NSPC                ! A number of model species
        INTEGER         NSPC_TURBINE        ! A number of model species for Turbine HAPs
        INTEGER         NSPC_PISTON         ! A number of model species for Piston HAPs

        INTEGER         NFLFILES            ! Max number of FLIGHT_FILELIST input files
        INTEGER         NSGFILES            ! Max number of SEGMENT_FILELIST input files

        REAL            CO, NOX, NO, NO2, HONO, SO2, POC, PEC, PSO4, TOG

        REAL         :: LTOALT = 0.0        ! LTO operations altitude (ft)
        REAL         :: CUTOFF = 0.0        ! cutoff altitude (ft) for modeling
        REAL         :: AVGELEV= 0.0        ! average of depart/arrival airport elevations

        REAL         :: FAC    = 0.0        ! tmp factor value 
        REAL         :: ALEN   = 0.0        ! link length
        REAL         :: LATVAL = 0.0        ! absolute latitude
        REAL         :: LONVAL = 0.0        ! absolute longitude
        REAL         :: PRVLAT = 0.0        ! previous absolute latitude
        REAL         :: PRVLON = 0.0        ! previous absolute longitude
        REAL         :: PRVHGT = 0.0        ! previous altitude
        REAL         :: HEIGHT = 0.0        ! altitude
        REAL         :: DELTAZ = 0.0        ! delta altitude
        REAL         :: SEG_TIME = 0.0      ! segment flight duratoin
        REAL         :: SEG_DIST = 0.0      ! segment flight distance
        REAL         :: FUELBURN = 0.0      ! fuel burn
        
        REAL            Z, Zo, Zh, ZBOT, ZTOP, PDIFF
        REAL            ZFRAC, PFRAC, VGTOP, LTOT
        REAL            TMPVAL

        LOGICAL      :: EFLAG = .FALSE.     ! true: ERROR

        CHARACTER(32)   TMPCHAR             
        CHARACTER(1056) LINE                ! input line buffer
        CHARACTER(500)  FLG_BUF             ! flight input line buffer
        CHARACTER(500)  SEG_BUF             ! segment input line buffer
        CHARACTER(256)  MESG                ! message buffer
        CHARACTER(16)   POLNAM              ! tmp pollutant name
        CHARACTER(16)   DPRTID              ! tmp departure airport code
        CHARACTER(16)   ARRVID              ! tmp arriving airport code
        CHARACTER(16)   FLGID               ! aircraft flight identifier
        CHARACTER(3)    TZONE               ! time zone
        CHARACTER(16)   COORUNIT            !  coordinate system projection units
        CHARACTER(80)   GDESC               !  grid description


        CHARACTER(16) :: PROGNAME = 'AEDTPROC' ! program name

C*************************************************************************
C   begin body of program AEDTPROC

        LDEV = INIT3()

C.........  Write copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get logical value from the environment
        MESG = 'Enter year difference between MCIP modeling year ' //
     &          'and AEDT modeling year [MCIP_YEAR-AEDT_YEAR]'
        DYEAR = ENVINT( 'YEAR_DIFF', MESG, 0, IOS )

        MESG = 'Determine the max altitude (feet) for LTO operations [default: 10,000]'
        LTOALT = ENVREAL( 'LTO_ALTITUDE', MESG, 10000.0, IOS )

        MESG = 'Determine the cutoff altitude (feet) for modeling [default: 50,000]'
        CUTOFF = ENVREAL( 'CUTOFF_ALTITUDE', MESG, 50000.0, IOS )
        
C........  Open unit numbers of input files
        MESG = 'Enter logical name for aircraft flight information file list'
        FDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'FLIGHT_FILELIST', PROGNAME )
        
        MESG = 'Enter logical name for aircraft segment input file list'
        SDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SEGMENT_FILELIST', PROGNAME )

C.........  stores HAP factors for Turbine engines
        MESG = 'Enter logical name for a file of airport elevation input file'
        APDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'APRT_ELEVATION',PROGNAME )
        CALL READ_AIRPORT_ELEVATION

C.........  stores HAP factors for Turbine engines
        MESG = 'Enter logical name for a file of a list of TURBIN E HAPs'
        TDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'TURBINE_HAP',PROGNAME )
        CALL READ_TURBINE_HAP_FACTOR

C.........  stores HAP factors for Piston engines
        MESG = 'Enter logical name for a file of a list of PISTON HAPs'
        PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'PISTON_HAP', PROGNAME )
        CALL READ_PISTON_HAP_FACTOR
        
C.........  Store Turbine engine speciation profiles : new 1098 developed by UNC-FAA
        MESG = 'Enter logical name for a file of a list of TURBINE speciation profiles'
        TSDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'TURBINE_SPC',PROGNAME )
        CALL READ_TURBINE_SPECIATION

C.........  Store Piston engine speciation profiles : 1099 (2011)
        MESG = 'Enter logical name for a file of a list of PISTON speciation profiles'
        PSDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'PISTON_SPC', PROGNAME )
        CALL READ_PISTON_SPECIATION

        MESG = 'Enter logical name for report file'
        RDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'REPORT', PROGNAME )

C.........  Open Met input files
        MESG = 'Enter name for CROSS-POINT GRID file'
        GNAME = PROMPTMFILE( MESG, FSREAD3, 'GRD_CRO_2D', PROGNAME )

        MESG = 'Enter name for CROSS-POINT LAYERED MET file'
        XNAME = PROMPTMFILE( MESG, FSREAD3, 'MET_CRO_3D', PROGNAME )


C.........  Get grid name from the environment and read grid parameters
        IF ( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD, GDTYP3D,
     &       COORUNIT, P_ALP3D, P_BET3D, P_GAM3D, XCENT3D,
     &       YCENT3D, XORIG3D, YORIG3D, XCELL3D,
     &       YCELL3D, NCOLS3D, NROWS3D, NTHIK3D)) THEN
 
            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........   Check grid information
        CALL CHKGRID( GDNAM3D, 'GRIDDESC', 1, EFLAG )

        IF ( EFLAG ) THEN
             MESG = 'Problem with gridded input data.'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Read description of 3d file for defining layer structure
        IF( .NOT. DESC3( XNAME ) ) THEN
             MESG = 'Could not get description of file "' //
     &               TRIM( XNAME ) // '"'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Store local layer info
        ALLOCATE( VGLVLS( 0:MXLAYS3 ), STAT= IOS)
        CALL CHECKMEM( IOS, 'VGLVLS', PROGNAME )

        SDATE  = SDATE3D
        EDATE  = SDATE3D
        STIME  = STIME3D
        ETIME  = STIME3D
        NSTEPS = MXREC3D
        NLAYS  = NLAYS3D
        VGTYP  = VGTYP3D
        VGTOP  = VGTOP3D
        VGLVLS = VGLVS3D   ! array
        
        DO T = 1,NSTEPS - 1
            CALL NEXTIME( EDATE, ETIME, 10000 )
        END DO

C.........  Determine beginning processing date
        YEAR = INT( SDATE/1000 )
        BDATE = SDATE - ( YEAR * 1000 )
        IF( BDATE == 1 ) THEN
            BDATE = (YEAR-1) * 1000 + JULIAN( (YEAR-1), 12, 31 )
        ELSE
            BDATE = SDATE - 1
        END IF
       print*,BDATE,SDATE,EDATE,YEAR,'B,S,Edate'
C.........  Allocate arrays
        ALLOCATE( SFCHGT( NGRID ), STAT= IOS)
        CALL CHECKMEM( IOS, 'SFCHGT', PROGNAME )
        ALLOCATE( TERRAIN( NGRID ), STAT= IOS)
        CALL CHECKMEM( IOS, 'TERRAIN', PROGNAME )
        ALLOCATE( ZZF( NGRID,NLAYS ), STAT= IOS)
        CALL CHECKMEM( IOS, 'ZZF', PROGNAME )
        ALLOCATE( TMP3D( NGRID,NLAYS,NVARS,NSTEPS ), STAT= IOS)
        CALL CHECKMEM( IOS, 'TMP3D', PROGNAME )
        ALLOCATE( LFRAC( NLAYS ), STAT= IOS)
        CALL CHECKMEM( IOS, 'LFRAC', PROGNAME )
        ALLOCATE( ACEL( NGRID ), STAT= IOS)
        CALL CHECKMEM( IOS, 'ACEL', PROGNAME )
        ALLOCATE( AFAC( NGRID ), STAT= IOS)
        CALL CHECKMEM( IOS, 'AFAC', PROGNAME )
        ACEL   = 0
        AFAC   = 0.0
        SFCHGT = 0.0
        TERRAIN= 0.0
        ZZF    = 0.0
        LFRAC  = 0.0
        TMP3D  = 0.0

C.............  Initialize I/O API output file headers
c        CALL HDRMISS3

        FDESC3D = ' '

        FDESC3D( 1 ) = 'AEDT aircraft emissions file'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ AEDT from FAA'
        FDESC3D( 4 ) = '/Output file from AEDTPROC ' //
     &      'Program that creates 3D CMAQ-ready emissions using AEDT '//
     &      '3D Chorded emissions'
c        FDESC3D( 21 ) = '/INVEN FROM/ ' // 'AEDT'
c        FDESC3D( 22 ) = '/INVEN VERSION/ ' // 'FAA-AEDT'

        SDATE3D = SDATE
        STIME3D = STIME
        NVARS3D = NVARS
        VTYPE3D = M3REAL
        UNITS3D = 'moles/s'

        UNITS3D( 6  )  = 'g/s'    ! POC in g/s
        UNITS3D( 7  )  = 'g/s'    ! PEC in g/s
        UNITS3D( 8  )  = 'g/s'    ! PSO4 in g/s
        UNITS3D( 43 )  = 'g/s'    ! TOG in original g/s (no MW is available for a proper conversion to mole
        
        VNAME3D( 1  ) = 'CO              '
        VDESC3D( 1  ) = 'Model species CO'
        VNAME3D( 2  ) = 'SO2             '
        VDESC3D( 2  ) = 'Model species SO2'
        VNAME3D( 3  ) = 'NO              '
        VDESC3D( 3  ) = 'Model species NO'
        VNAME3D( 4  ) = 'NO2             '
        VDESC3D( 4  ) = 'Model species NO2'
        VNAME3D( 5  ) = 'HONO            '
        VDESC3D( 5  ) = 'Model species HONO'
        VNAME3D( 6  ) = 'POC             '
        VDESC3D( 6  ) = 'Model species POC'
        VNAME3D( 7  ) = 'PEC             '
        VDESC3D( 7  ) = 'Model species PEC'
        VNAME3D( 8  ) = 'PSO4            '
        VDESC3D( 8  ) = 'Model species PSO4'
        VNAME3D( 9  ) = 'ALD2            '
        VDESC3D( 9  ) = 'Model species ALD2'
        VNAME3D( 10 ) = 'ALDX            '
        VDESC3D( 10 ) = 'Model species ALDX'
        VNAME3D( 11 ) = 'CH4             '
        VDESC3D( 11 ) = 'Model species CH4'
        VNAME3D( 12 ) = 'ETH             '
        VDESC3D( 12 ) = 'Model species ETH'
        VNAME3D( 13 ) = 'ETHA            '
        VDESC3D( 13 ) = 'Model species ETHA'
        VNAME3D( 14 ) = 'ETOH            '
        VDESC3D( 14 ) = 'Model species ETOH'
        VNAME3D( 15 ) = 'FORM            '
        VDESC3D( 15 ) = 'Model species FORM'
        VNAME3D( 16 ) = 'IOLE            '
        VDESC3D( 16 ) = 'Model species IOLE'
        VNAME3D( 17 ) = 'MEOH            '
        VDESC3D( 17 ) = 'Model species MEOH'
        VNAME3D( 18 ) = 'OLE             '
        VDESC3D( 18 ) = 'Model species OLE'
        VNAME3D( 19 ) = 'PAR             '
        VDESC3D( 19 ) = 'Model species PAR'
        VNAME3D( 20 ) = 'TOL             '
        VDESC3D( 20 ) = 'Model species TOL'
        VNAME3D( 21 ) = 'UNR             '
        VDESC3D( 21 ) = 'Model species UNR'
        VNAME3D( 22 ) = 'XYL             '
        VDESC3D( 22 ) = 'Model species XYL'
        VNAME3D( 23 ) = 'FORM_PRIMARY    '
        VDESC3D( 23 ) = 'Model species FORM_PRIMARY'
        VNAME3D( 24 ) = 'ALD2_PRIMARY    '
        VDESC3D( 24 ) = 'Model species ALD2_PRIMARY'
        VNAME3D( 25 ) = 'BENZENE         '
        VDESC3D( 25 ) = 'Model species BENZENE'
        VNAME3D( 26 ) = 'TOLUENE         '
        VDESC3D( 26 ) = 'Model species TOLUENE'
        VNAME3D( 27 ) = 'ACROLEIN        '
        VDESC3D( 27 ) = 'Model species ACROLEIN'
        VNAME3D( 28 ) = 'ACROLEIN_PRIMARY'
        VDESC3D( 28 ) = 'Model species ACROLEIN_PRIMARY'
        VNAME3D( 29 ) = 'BUTADIENE13     '
        VDESC3D( 29 ) = 'Model species BUTADIENE13'
        VNAME3D( 30 ) = 'PXYL            '
        VDESC3D( 30 ) = 'Model species PXYL'
        VNAME3D( 31 ) = 'OXYL            '
        VDESC3D( 31 ) = 'Model species OXYL'
        VNAME3D( 32 ) = 'MXYL            '
        VDESC3D( 32 ) = 'Model species MXYL'
        VNAME3D( 33 ) = 'NAPHTHALENE       '
        VDESC3D( 33 ) = 'Model species NAPHTHALENE'
        VNAME3D( 34 ) = 'PROPIONAL       '
        VDESC3D( 34 ) = 'Model species PROPIONAL'
        VNAME3D( 35 ) = 'ETHYLBENZ       '
        VDESC3D( 35 ) = 'Model species ETHYBENZ'
        VNAME3D( 36 ) = 'STYRENE         '
        VDESC3D( 36 ) = 'Model species STRYENE'
        VNAME3D( 37 ) = 'PHENOL          '
        VDESC3D( 37 ) = 'Model species PHENOL'
        VNAME3D( 38 ) = 'METHANOL        '
        VDESC3D( 38 ) = 'Model species METHANOL'
        VNAME3D( 39 ) = 'CUMENE          '
        VDESC3D( 39 ) = 'Model species CUMEME'
        VNAME3D( 40 ) = 'MTHYLNAP2       '
        VDESC3D( 40 ) = 'Model species MTHYLNAP2'
        VNAME3D( 41 ) = 'BENZALD         '
        VDESC3D( 41) = 'Model species BENZALD'
        VNAME3D( 42 ) = 'ETHYLENE        '
        VDESC3D( 42 ) = 'Model species ETHYLENE'
        VNAME3D( 43 ) = 'TOG_ORG         '
        VDESC3D( 43 ) = 'Model species TOG Original'

        MESG = 'Enter logical name for output file'
        ONAME = PROMPTMFILE( MESG, FSUNKN3, 'OUTPUT', PROGNAME )

C.............  Open new file
        IF( .NOT. OPEN3( 'OUTPUT', FSUNKN3, PROGNAME ) ) THEN
            MESG = 'Could not create new output file OUTPUT'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Reading a input list file
        NFLFILES = GETFLINE( FDEV, 'FLIGHT_FILELIST input file' )
        NSGFILES = GETFLINE( SDEV, 'SEGMENT_FILELIST input file' )

        IF( NFLFILES /= NSGFILES ) THEN
            MESG = 'ERROR: No of files from both FLIGHT_FILELIST '//
     &            'and SEGMENT_FILELIST must be same.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Process files in input list
        I = 0
        J = 0
        NF = 0

        DO N = 1, NFLFILES

C.............  Read AEDT Flight info file list
            READ( FDEV, 93000, IOSTAT=IOS ) FLG_BUF
            J = J + 1

C.............  Check for i/o errors
            IF( IOS /= 0 ) THEN
                WRITE( MESG,94010 ) 'I/O error', IOS,
     &              'reading input file at line', J
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C............. parse line to compute julian date          
            CALL PARSLINE( FLG_BUF, 4, INP_SEG )

            MON = STR2INT( TRIM( INP_SEG( 1 ) ) )
            DAY = STR2INT( TRIM( INP_SEG( 2 ) ) )
            YEAR= STR2INT( TRIM( INP_SEG( 3 ) ) )

            YEAR = YEAR + DYEAR  ! Override modeling year by user

            JDATE = YEAR * 1000 + JULIAN( YEAR, MON, DAY )

C.............  Read AEDT segment input file list
            READ( SDEV, 93000, IOSTAT=IOS ) SEG_BUF
            I = I + 1

C.............  Check for i/o errors
            IF( IOS /= 0 ) THEN
                WRITE( MESG,94010 ) 'I/O error', IOS,
     &              'reading input file at line', I
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip non-processing dates in flielist file
            IF( JDATE < BDATE ) CYCLE
            IF( JDATE > EDATE ) CYCLE
            NF = NF + 1

C.............  Open AEDT flight info input file      
            OPEN( FLDEV, FILE = INP_SEG( 4 ), STATUS='OLD', IOSTAT=IOS )

            IF( IOS /= 0 ) THEN
                MESG = 'Could not open file:' // 
     &                 CRLF() // BLANK5 // TRIM( FLG_BUF )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            MESG = 'Successfully opened AEDT Flight information input file:' //
     &             CRLF() // BLANK5 // TRIM( FLG_BUF )
            CALL M3MSG2( MESG )

C.............  Store AEDT flight information (flightID and Engine types)
            CALL READ_FLIGHT_ENGINE( FLDEV )
            
            CLOSE( FLDEV )

C............. parse line to compute julian date          
            CALL PARSLINE( SEG_BUF, 4, INP_SEG )

            MON = STR2INT( TRIM( INP_SEG( 1 ) ) )
            DAY = STR2INT( TRIM( INP_SEG( 2 ) ) )
            YEAR= STR2INT( TRIM( INP_SEG( 3 ) ) )

            YEAR = YEAR + DYEAR  ! Override modeling year by user

            TDATE = YEAR * 1000 + JULIAN( YEAR, MON, DAY )
            
            IF( TDATE /= JDATE ) THEN
                MESG = 'ERROR: Must process same date for both Flight'//
     &              ' information and emission input files'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Open actual Segment input file      
            OPEN( SGDEV, FILE = INP_SEG( 4 ), STATUS='OLD', IOSTAT=IOS )

            IF( IOS /= 0 ) THEN
                MESG = 'Could not open file:' // 
     &                 CRLF() // BLANK5 // TRIM( SEG_BUF )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            MESG = 'Successfully opened AEDT Flight Emission input file:'//
     &             CRLF() // BLANK5 // TRIM( SEG_BUF )
            CALL M3MSG2( MESG )

c            N_SEG = GETFLINE( SGDEV, 'Flight Emissions input file' )

C.............  Initialize values before reading file
            IREC  = 0
            PRVLAT = 0.0
            PRVLON = 0.0
            PRVHGT = 0.0

C.............  Read through input file          
            DO
            
                READ( SGDEV, 93000, END=299 ) LINE

                IREC = IREC + 1

                IF( MOD( IREC,1000000 ) == 0 ) THEN
                    WRITE( MESG, 94010 ) 'Processing line at', IREC
                    CALL M3MSG2( MESG )
                END IF

                IF( IOS > 0 ) EFLAG = .TRUE.

C.................  Skip blank/comment/non-HOUREMIS lines
                IF( IREC == 1 ) CYCLE

C.................  Read data from LINE
                CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.................  Skip when it is out of range of episode dates
                FLGID = TRIM( SEGMENT( 1 ) )       ! Flight ID
                SEGID = STR2INT( SEGMENT( 2 ) )    ! segment ID for flight
                MODID = STR2INT( SEGMENT( 12 ) )    ! segment ID for flight
                
                MON   = STR2INT( SEGMENT( 6 ) )    ! integer current month
                DAY   = STR2INT( SEGMENT( 7 ) )    ! integer current hour
                HOUR  = STR2INT( SEGMENT( 9 ) )    ! integer current hour
                YEAR  = STR2INT( SEGMENT( 8 ) )    ! integer current year

                YEAR = YEAR + DYEAR  ! Override modeling year by user
                
                T = HOUR + 1                       ! starting from 1 to 23
                JDATE = YEAR * 1000 + JULIAN( YEAR, MON, DAY )  ! segment julian date
                JTIME = HOUR * 10000

                SEG_TIME =  3600.0                 ! one hours duration (seconds)

C.................  store lat/lon coordinates (starting point)
                LATVAL = STR2REAL( SEGMENT( 3 ) )
                LONVAL = STR2REAL( SEGMENT( 4 ) )
                HEIGHT = FT2M * STR2REAL( SEGMENT( 5 ) )   ! altitude in unit of feet

C.................  Convert mobile source coordinates from lat-lon to output grid
                CALL CONVRTXY( 1, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,
     &                         XCENT, YCENT, LONVAL, LATVAL )

C.................  define start/end coordinates for each link
                IF( SEGID == 0 ) THEN
                    PRVLAT = LATVAL
                    PRVLON = LONVAL
                    PRVHGT = HEIGHT
                END IF

                Zo     = PRVHGT                ! origin height
                Zh     = HEIGHT                ! end height

C.................  If source is in the domain, get cell number and store
                ORG_CELLID = 0
                END_CELLID = 0
                IF( INGRID( PRVLON,PRVLAT,NCOLS,NROWS,COL,ROW ) ) THEN
                    ORG_CELLID = COL + ( ROW - 1 ) * NCOLS
                END IF
                
                IF( INGRID( LONVAL,LATVAL,NCOLS,NROWS,COL,ROW ) ) THEN
                    END_CELLID = COL + ( ROW - 1 ) * NCOLS
                END IF

C.................  Skip non-processing dates
                IF( ( JDATE <  SDATE ) .OR. ( JDATE > EDATE ) .OR.
     &              ( JDATE == EDATE .AND. JTIME > ETIME )        ) THEN
                    PRVLAT = LATVAL  ! store previous lat coordinate
                    PRVLON = LONVAL  ! store previous lon coordinate
                    PRVHGT = HEIGHT  ! store previous height
                    CYCLE
                END IF

C.................  If link source, determine the number of cells for this source
                NCEL = 0
                ACEL = 0
                AFAC = 0.0
                CALL LNK2GRD( NGRID, PRVLON, PRVLAT, LONVAL, LATVAL,
     &                        NCEL, ACEL, AFAC, ALEN, EFLAG)

C.................   previous values needed to be updated before the source gets skipped
                PRVLAT = LATVAL  ! store previous lat coordinate
                PRVLON = LONVAL  ! store previous lon coordinate
                PRVHGT = HEIGHT  ! store previous height

C.................  Make sure that there was enough storage
                IF( EFLAG ) THEN
                    WRITE( MESG,94010 )
     &                  'INTERNAL ERROR: Overflow for source at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
C.................  Skip if there is no intersected grid cell
                IF( NCEL == 0 ) CYCLE

C.................  Vertical allocation within gridded x-y cell(s)
C                   Need to compute each gridded x-y cell actual bottom and top heights
C                   and then vertically allocate gridded x-y link into x-z layers 
C                   using ZF layer full heights

C.................  Read layers top height (meter)
                IF( .NOT. READ3( XNAME, 'ZF', -1,
     &                           JDATE, JTIME, ZZF ) ) THEN
                    MESG = 'ERROR : Could not read ZF from file '
     &                     // TRIM( XNAME )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Read terrain height (m) to convert MSL to AGL
                IF( .NOT. READ3( GNAME, 'HT', -1,
     &                           JDATE, JTIME, TERRAIN ) ) THEN
                    MESG = 'ERROR : Could not read TERRAIN from file '
     &                     // TRIM( GNAME )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Convert MSL altitude to AGL (MSL-Terrain Height)
C                   Subtract airport elevation when altitude is below 3000ft (LTO mode)
C                   Subtract terrain height when altitude is higher than 3000ft
C                   Over 3000ff, calculate sigma level using surface pressure and 
C                   estimated pressure based on calculated altitude
C                   1000ft * 0.3048 = 3048 meter

                F = INDEX1( FLGID, N_FLG, FLIGHT_ENG( :,1 ) )  ! Flight ID
                DPRTID = FLIGHT_ENG( F,2 )                     ! Departure Airport ID
                ARRVID = FLIGHT_ENG( F,3 )                     ! Departure Airport ID
                    
                ND = INDEX1( DPRTID, NAPRT_ELEV, APRT_CODE )
                NA = INDEX1( ARRVID, NAPRT_ELEV, APRT_CODE )

                IF( Zo <= LTOALT * FT2M ) THEN

                    IF( MODID < 4 ) THEN
                        Zo = Zo - APRT_ELEV( ND )   ! departure airport elev
                    ELSE IF( MODID > 7 ) THEN
                        Zo = Zo - APRT_ELEV( NA )   ! arrival airport elev
                    ELSE
                        AVGELEV=( APRT_ELEV(ND) + APRT_ELEV(NA) ) / 2.0  !avg elev
                        Zo = Zo - AVGELEV
                    ENDIF

                ELSE
                    Zo = Zo - TERRAIN( ORG_CELLID )
                END IF

                IF( Zh <= LTOALT * FT2M ) THEN
                    IF( MODID < 4 ) THEN
                        Zh = Zh - APRT_ELEV( ND )   ! departure airport elev
                    ELSE IF( MODID > 7 ) THEN
                        Zh = Zh - APRT_ELEV( NA )   ! arrival airport elev
                    ELSE
                        AVGELEV=( APRT_ELEV(ND) + APRT_ELEV(NA) ) / 2.0  !avg elev
                        Zh = Zh - AVGELEV
                    ENDIF
                ELSE
                    Zh = Zh - TERRAIN( END_CELLID )
                END IF

                IF( Zo < 0.0 ) Zo = 0.0
                IF( Zh < 0.0 ) Zh = 0.0

                ZBOT = Zo

                DELTAZ = Zh - Zo   ! delta z (<0:langind, >0:climbing)

C.................  Sort the order of grid cell processing. Starting from org_cellid....
                IF( ACEL( 1 ) /= ORG_CELLID ) THEN
                    STRID = NCEL
                    ENDID = 1
                    INCID = -1
                ELSE
                    STRID = 1
                    ENDID = NCEL
                    INCID = 1
                END IF


                DO NC = STRID, ENDID, INCID

                    LFRAC = 0.0         ! Gridded x-y link vertical fractions 

                    C     = ACEL( NC )  ! cell id
                    ZFRAC = AFAC( NC )  ! cell-fraction value

                    IF( ZFRAC < 0.0 ) THEN
                        MESG = 'ERROR: Can not process negative x-y '//
     &                      'link fraction'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    ENDIF

                    Z = ZFRAC * DELTAZ  ! actual vertical length within single x-y cell

C.....................  Update bottom and top layers
                    IF( DELTAZ < 0 ) THEN   ! aircraft landing mode
                        ZTOP = ZBOT 
                        ZBOT = ZTOP + Z
                    ELSE                    ! aircraft climbing mode
                        ZBOT = ZBOT
                        ZTOP = ZBOT + Z
                    END IF

C.....................  retrieve col/row from cellid using C=(ROW-1)*NCOLS+COL
                    ROW = C / NCOLS
                    IF( MOD( C, NCOLS ) .GT. 0. ) ROW = ROW + 1
                    COL = C - ( ROW-1 ) * NCOLS

C.....................  Looping through layers to determine associated layer for each link
                    DO L = 1, NLAYS - 1

                        IF ( ZBOT <= ZZF( C,L ) ) THEN
                            LBOT = L
                            GO TO  111   ! end loop and skip reset of LBOT
                        END IF

                    END DO

                    LBOT = NLAYS           !  fallback

111                 CONTINUE                !  loop exit:  bottom found at LBOT
 
                    IF ( ZTOP <= ZZF( C,LBOT ) ) THEN  !  plume in this layer
 
                        PFRAC = 1.0
                        LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
                        LTOP = LBOT
c                write(*,'(9i8,f15.7,a)')segid,c,org_cellid,end_cellid,col,row,LBOT,LTOP,LBOT,pfrac,' onelayer'

                    ELSE IF( LBOT == NLAYS ) THEN    ! plume above top layer
 
                        PFRAC = 1.0
                        LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
                        LTOP = NLAYS
c                write(*,'(9i8,f15.7,a)')segid,c,org_cellid,end_cellid,col,row,LBOT,LTOP,LTOP,pfrac,' toplayer'
                
                    ELSE                               ! plume crosses layers
 
                        DO L = LBOT + 1, NLAYS
                            IF ( ZTOP <= ZZF( C,L ) ) THEN
                                LTOP = L
                                GO TO 222  ! end loop and skip reset of LTOP
                            END IF
                        END DO
                        LTOP = NLAYS

222                     CONTINUE
 
C.........................  Calculate between layer 
                        PDIFF = ZTOP - ZBOT
                    
C.........................  Calculate a fraction for the bottom layer
                        PFRAC = ( ( ZZF( C,LBOT ) - ZBOT )
     &                              / PDIFF )
                        LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
c                write(*,'(9i8,f15.7,a)')segid,c,org_cellid,end_cellid,col,row,LBOT,LTOP,LBOT,pfrac
                    
C.........................  Calculate a fraction for the top layer
                        PFRAC = ( ( ZTOP - ZZF( C,LTOP-1 ) )
     &                             / PDIFF )
                        LFRAC( LTOP ) = LFRAC( LTOP ) + PFRAC

                        DDP = LTOP - LBOT

                        IF( DDP >= 2 ) THEN

                            DO L = LBOT+1, LTOP-1 !  layers in plume
                            
                                PFRAC = ( ( ZZF(C,L) -ZZF(C,L-1) )
     &                                      / PDIFF )
c                write(*,'(9i8,f15.7,a)')segid,c,org_cellid,end_cellid,col,row,LBOT,LTOP,L,pfrac
                                LFRAC( L ) = LFRAC( L ) + PFRAC
                            END DO

                        ENDIF
c                write(*,'(9i8,f15.7,a)')segid,c,org_cellid,end_cellid,col,row,LBOT,LTOP,LTOP,lfrac(ltop)
                    
                    END IF
c       write(*,'(a,7I8,3F15.8)')trim(FLGID),SEGID,c,col,row,nc,lbot,ltop,zbot,ztop,zfrac
c       print*,zbot,ztop,height+deltaz,height,deltaz

                    IF( DELTAZ > 0 ) THEN
                        ZBOT = ZTOP   ! landing mode: ZBOT needs to be updated with ZTOP
                    ELSE
                        ZTOP = ZBOT
                    END IF

C..................... Before applying layer fractions make sure that they add to 1.0
                    LTOT = 0.0
                    DO NL = 1, NLAYS
                        LTOT = LTOT + LFRAC( NL )
                    ENDDO

                    IF ( LTOT < 0.999 ) THEN
                        MESG = 'ERROR: Total of layer fractions '//
     &                         'are less than 1.0.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )                  
                    ELSE IF ( LTOT > 1.001 ) THEN
                        MESG = 'ERROR: Total of layer fractions '//
     &                         'are greater than 1.0.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )                  
                    ENDIF

C.....................  Loop over allocated layers by link
                    DO L = LBOT, LTOP
c         write(*,'(4i5,f10.7)')segid,col,row,L,LFRAC(l)
C.........................  SPECIATION PROCESSING
C.........................  Store emission factors by pollutant
                        CO   = STR2REAL( SEGMENT( 22 ) ) / SEG_TIME
                        NOX  = STR2REAL( SEGMENT( 24 ) ) / SEG_TIME
                        PEC  = STR2REAL( SEGMENT( 25 ) ) / SEG_TIME
                        POC  = STR2REAL( SEGMENT( 27 ) ) / SEG_TIME
                        TOG  = 1.16 * ( STR2REAL( SEGMENT( 23 ) ) / SEG_TIME )

C.........................  Compute SO2 and SO4 using FSC and E
                        FUELBURN = STR2REAL( SEGMENT( 21 ) ) / SEG_TIME
                        SO2  = (FSC/1000.) * ((100.-E)/100.) * FUELBURN * (64/32)
                        PSO4 = (FSC/1000.) * (E/100.)        * FUELBURN * (96/32)

C.........................  Estimated layer-emissions multiply by layer fractions * cell fraction
                        CO   = CO   * ZFRAC * LFRAC( L )
                        NOX  = NOX  * ZFRAC * LFRAC( L )
                        PEC  = PEC  * ZFRAC * LFRAC( L )
                        PSO4 = PSO4 * ZFRAC * LFRAC( L )
                        POC  = POC  * ZFRAC * LFRAC( L )
                        SO2  = SO2  * ZFRAC * LFRAC( L )
                        TOG  = TOG  * ZFRAC * LFRAC( L )

C.........................  conversion from mass to mole
                        CO   = CO  / 28.0
                        NOX  = NOX / 46.0   ! NOx is in NO2 equilvalent
                        SO2  = SO2 / 64.0

C.........................  SPECIATION Section
C.........................  Apply 3D gridded link factors to pollutants for speciation
                        TMP3D( C,L,1,T ) = TMP3D( C,L,1,T ) + CO               ! CO in mole/s
                        TMP3D( C,L,2,T ) = TMP3D( C,L,2,T ) + SO2              ! SO2 in mole/s

                        IF( HEIGHT <= LTOALT * FT2M ) THEN
                            TMP3D( C,L,3,T ) = TMP3D( C,L,3,T ) + 0.76 * NOX   ! NO in mole/s
                            TMP3D( C,L,4,T ) = TMP3D( C,L,4,T ) + 0.23 * NOX   ! NO2 in mole/s
                            TMP3D( C,L,5,T ) = TMP3D( C,L,5,T ) + 0.01 * NOX   ! HONO in mole/s
                            TMP3D( C,L,6,T ) = TMP3D( C,L,6,T ) + 1.0  * POC   ! POC in unit of g/s
                            TMP3D( C,L,7,T ) = TMP3D( C,L,7,T ) + 1.0  * PEC   ! PEC in unit of g/s

                        ELSE
                            TMP3D( C,L,3,T ) = TMP3D( C,L,3,T ) + 0.90 * NOX   ! NO
                            TMP3D( C,L,4,T ) = TMP3D( C,L,4,T ) + 0.09 * NOX   ! NO2
                            TMP3D( C,L,5,T ) = TMP3D( C,L,5,T ) + 0.01 * NOX   ! HONO
                            TMP3D( C,L,6,T ) = TMP3D( C,L,6,T ) + 0.03 * FUELBURN   ! POC in unit of g/s
                            TMP3D( C,L,7,T ) = TMP3D( C,L,7,T ) + 0.03 * FUELBURN   ! PEC in unit of g/s
                        END IF 

                        TMP3D( C,L,8, T ) = TMP3D( C,L, 8,T ) + PSO4           ! POS4 in unit of g/s
                        TMP3D( C,L,43,T ) = TMP3D( C,L,43,T ) + TOG            ! TOG in original g/s

C.........................  Compute model species in mole/s from TOG
C.........................  Determine whether engine is piston or turbine
                        F = INDEX1( FLGID, N_FLG, FLIGHT_ENG( :,1 ) )

C.........................  Only two piston engine out of entire engine types available
                        IF( FLIGHT_ENG( F,4 ) == 'IO360' .OR.
     &                      FLIGHT_ENG( F,4 ) == 'R1820'       ) THEN

C.............................  Estimate Model species using speciation profile factors from TOG
                            DO NS = 1, NSPC_PISTON

                                POLNAM =  NSSPCP( NS )

                                NV = INDEX1( POLNAM, NVARS, VNAME3D )

                                IF( NV < 1 ) THEN
                                    MESG = 'ERROR: ' // POLNAM //
     &                                  ' is not a model species:: ' //
     &                                  'PLEASE check the name of species'
                                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                                END IF

                                TMPVAL = CVSPCP( NS,1 ) / CVSPCP( NS,2 ) * TOG
                                TMP3D( C,L,NV,T ) = TMP3D( C,L,NV,T ) + TMPVAL

                            END DO    ! end of loop

C.............................  Estimate HAPs using HAP's profile factors from TOG
                            DO NP = 1, NFAC_PISTON

                                POLNAM = NPHAPP( NP )

                                NV = INDEX1( POLNAM, NVARS, VNAME3D )

                                IF( NV < 1 ) THEN
                                    MESG = 'ERROR: ' // POLNAM //
     &                                  ' is not a model species:: ' //
     &                                  'PLEASE check the name of species'
                                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                                END IF
                                
                                TMPVAL = CVHAPP( NP,1 ) / CVHAPP( NP,2 ) * TOG
                                TMP3D( C,L,NV,T ) = TMP3D( C,L,NV,T ) + TMPVAL

                            END DO    ! end of loop

C..........................  processing turbine engine types for TOG speciation
                        ELSE

C.............................  Estimate Model species using speciation profile factors from TOG
                            DO NS = 1, NSPC_TURBINE

                                POLNAM =  NSSPCT( NS )

                                NV = INDEX1( POLNAM, NVARS, VNAME3D )

                                IF( NV < 1 ) THEN
                                    MESG = 'ERROR: ' // POLNAM //
     &                                  ' is not a model species:: ' //
     &                                  'PLEASE check the name of species'
                                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                                END IF

                                TMPVAL = CVSPCT( NS,1 ) / CVSPCT( NS,2 ) * TOG
                                TMP3D( C,L,NV,T ) = TMP3D( C,L,NV,T ) + TMPVAL

                            END DO    ! end of loop


C.............................  Estimate HAPs using HAP's profile factors from TOG
                            DO NP = 1, NFAC_TURBINE
                            
                                POLNAM = NPHAPT( NP )

                                NV = INDEX1( POLNAM, NVARS, VNAME3D )

                                IF( NV < 1 ) THEN
                                    MESG = 'ERROR: ' // POLNAM //
     &                                  ' is not a model species:: ' //
     &                                  'PLEASE check the name of species'
                                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                                END IF

                                TMPVAL = CVHAPT( NP,1 ) / CVHAPT( NP,2 ) * TOG
                                TMP3D( C,L,NV,T ) = TMP3D( C,L,NV,T ) + TMPVAL

                            END DO    ! end of loop

                        END IF

                    END DO     ! end of layer loop per link

                    DO V = 1, NVARS
                    DO L = 1, NLAYS
                        WRITE( TMPCHAR,'(F32.16)') TMP3D(C,L,V,T)
                        IF( .NOT. CHKREAL( TMPCHAR ) ) THEN
                            TMP3D( C,L,V,T ) = 0.0
                        END IF
                    END DO
                    END DO

                END DO      ! end of link loop

            ENDDO       ! end of flight loop

299         CONTINUE ! Exit from read loop

            CLOSE( FLDEV )
            IF( EFLAG ) THEN 
                WRITE( MESG,94010 ) 'I/O error', IOS,
     &               'reading input file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END DO     ! end of flight segment input files.

        IF( NF == 0 ) THEN
            MESG = 'ERROR: No flight emission files are processed :: '
     &          // 'CHECK FLIGHT_FILELIST and SEGMENT_FILELIST input files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get output file name using environment variable
C.........  Define top layer for output file
        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS
        
            DO V = 1, NVARS

                POLNAM = VNAME3D( V )
            
                IF ( .NOT. WRITE3( ONAME, POLNAM, JDATE, JTIME, 
     &                            TMP3D( :,:,V,T ) ) ) THEN
                    WRITE( MESG, 93000 ) 'Could not write to "'
     &                    // TRIM( ONAME ) // '".'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            END DO

            CALL NEXTIME( JDATE, JTIME, 10000 )

        END DO

        IF ( .NOT. CLOSE3( ONAME ) ) THEN
            CALL M3ERR( PROGNAME, 0, 0,
     &          TRIM( ONAME ) // '".', .TRUE. )
        END IF      !  if close3() failed 

C.........   End of program:
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C.........  Formatted file I/O formats...... 93xxx
93000   FORMAT( A )
93010   FORMAT( I2.2,I3.3,',"',A,'","',I5,'","',I10,'",,"","',A10,
     &        '",,,"',I10,'",,,,,,,,"L",',F9.4,',',F9.4,',',I2,
     &        ',"',A,'",,,,,,,,,,,,,,,,,,' )
93020   FORMAT(I2.2, I3.3, A15, 2I12,12X, A5, A8, A3, 24E7.1,E8.2,1X,A10,
     &       1X, A16)

C.......  Internal buffering formats...... 94xxx
94010   FORMAT( 10 ( A, :, I8, :, 2X  ) )
94011   FORMAT( 10 ( A, :, F10.3, :, 2X  ) )
94012   FORMAT( 50F15.3 )

C******************  INTERNAL SUBPROGRAMS  *****************************

       CONTAINS

C      This subroutine stores a list of unique engine types

       SUBROUTINE READ_FLIGHT_ENGINE( ADEV )

C..........  Subroutine arguments
         INTEGER, INTENT ( IN ) :: ADEV

C..........   Local variables
         INTEGER         I, N                  ! indices and counters

         INTEGER      :: NLINES = 0            ! number of lines in input file

         CHARACTER(256)  LINE                  ! Read buffer for a line
         CHARACTER(300)  MESG                  ! Message buffer
         CHARACTER(60)   SEGMENT( MXSEG )      ! line parsing array

C......................................................................

C...........  Get the number of lines
         NLINES = GETFLINE( ADEV, 'AEDT Flight information input file' )

C..........  Count a list of unique AEDT Engine types 
         N_FLG = 0
         IREC  = 0
         DO I = 1, NLINES

             READ ( ADEV, 93000, IOSTAT=IOS ) LINE
             IREC = IREC + 1

             IF ( IOS .GT. 0 ) THEN
                 WRITE( MESG, 94010)
     &                 'I/O error', IOS, 'reading ICAO_FIPS '//
     &                 'description file at line', IREC
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             END IF

C.............  Skip blank and comment lines
             IF( LINE( 1:19 ) == 'FLIGHT_ID' ) CYCLE

             N_FLG = N_FLG + 1

         END DO    ! end of loop

         IF( N_FLG == 0 ) THEN 
             MESG = 'ERROR: No valid entries for AEDT Engine types'
             CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
         END IF

         REWIND( ADEV )

C...........  Determine number of lines in filelist; this will be the maximum
C             number of airport sources
         IF( ALLOCATED( FLIGHT_ENG ) ) DEALLOCATE( FLIGHT_ENG )
         ALLOCATE( FLIGHT_ENG( N_FLG,4 ), STAT=IOS )
         CALL CHECKMEM( IOS, 'FLIGHT_ENG', PROGNAME )
         FLIGHT_ENG = ' '

C.........  Read the AEDT Engine type file and find 

         N_FLG = 0
         IREC  = 0
        
         DO I = 1, NLINES

             READ ( ADEV, 93000, IOSTAT=IOS ) LINE
             IREC = IREC + 1

             IF ( IOS .GT. 0 ) THEN
                 WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading ICAO_FIPS '//
     &                'description file at line', IREC
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             END IF

C.............  Left adjust line
             LINE = TRIM( LINE )

C.............  Skip blank and comment lines
             IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
             CALL PARSLINE( LINE, MXSEG, SEGMENT )

             N_FLG = N_FLG + 1
             FLIGHT_ENG( N_FLG,1 ) =  ADJUSTL( SEGMENT( 1 ) )  ! flight IDs
             FLIGHT_ENG( N_FLG,2 ) =  ADJUSTL( SEGMENT( 2 ) )  ! Airport IDs
             FLIGHT_ENG( N_FLG,3 ) =  ADJUSTL( SEGMENT( 4 ) )  ! Airport IDs
             FLIGHT_ENG( N_FLG,4 ) =  ADJUSTL( SEGMENT( 8 ) )  ! Engine types

         END DO    ! end of loop

         RETURN
        
C...................  FORMAT  STATEMENTS   ............................

C.......  Formatted file I/O formats...... 93xxx
93000    FORMAT( A )

C.......  Internal buffering formats...... 94xxx
94010    FORMAT( 10 ( A, :, I8, :, 2X  ) )

       END SUBROUTINE READ_FLIGHT_ENGINE

C******************  INTERNAL SUBPROGRAMS  *****************************

C  Read and store airport elevations to compute AGL of below 3000ft LTO sources

       SUBROUTINE READ_AIRPORT_ELEVATION

C..........   Local variables
         INTEGER         I, N                  ! indices and counters

         INTEGER      :: NLINES = 0            ! number of lines in input file

         CHARACTER(256)  LINE                  ! Read buffer for a line
         CHARACTER(300)  MESG                  ! Message buffer
         CHARACTER(16 ):: SEGMENT( 3 ) = ' '   ! line parsing array

C......................................................................

C.........  Get the number of lines
         NLINES = GETFLINE( APDEV, 'Airport elevation input file' )

C...........  Determine number of lines in filelist; this will be the maximum
C             number of airport sources
         ALLOCATE( APRT_CODE( NLINES ), STAT=IOS )
         CALL CHECKMEM( IOS, 'APRT_CODE', PROGNAME )
         ALLOCATE( APRT_ELEV( NLINES ), STAT=IOS )
         CALL CHECKMEM( IOS, 'APRT_ELEV', PROGNAME )

         APRT_CODE = ' '
         APRT_ELEV = 0.0
        
         N = 0
         IREC  = 0

         DO I = 1, NLINES

             READ ( APDEV, 93000, IOSTAT=IOS ) LINE
             IREC = IREC + 1

             IF ( IOS .GT. 0 ) THEN
                 WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading APRT_ELEVATION '//
     &                'description file at line', IREC
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             END IF

C.............  Left adjust line
             LINE = TRIM( LINE )

C.............  Skip blank and comment lines
             IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
             CALL PARSLINE( LINE, 3, SEGMENT )

             N = N + 1
             APRT_CODE( N ) = SEGMENT( 1 )
             APRT_ELEV( N ) = STR2REAL( SEGMENT( 3 ) )  * FT2M

         END DO    ! end of loop

         IF( N == 0 ) THEN 
             MESG = 'ERROR: No entries of APRT_ELEVATION'
             CALL M3MSG2( MESG )
         END IF

         NAPRT_ELEV = N

         CLOSE( APDEV )

         RETURN

C...................  FORMAT  STATEMENTS   ............................

C.......  Formatted file I/O formats...... 93xxx
93000    FORMAT( A )

C.......  Internal buffering formats...... 94xxx
94010    FORMAT( 10 ( A, :, I8, :, 2X  ) )

       END SUBROUTINE READ_AIRPORT_ELEVATION

C******************  INTERNAL SUBPROGRAMS  *****************************
       
C      This subroutine stores a list of aircraft AEDT HAPs to convert 
C      to species mass emission based upon TOG emissions

       SUBROUTINE READ_TURBINE_HAP_FACTOR

C..........   Local variables
         INTEGER         I, N                  ! indices and counters

         INTEGER      :: NLINES = 0            ! number of lines in input file

         CHARACTER(256)  LINE                  ! Read buffer for a line
         CHARACTER(300)  MESG                  ! Message buffer
         CHARACTER(16 ):: SEGMENT( 3 ) = ' '   ! line parsing array

C......................................................................

C.........  Get the number of lines
         NLINES = GETFLINE( TDEV, 'TURBINE_HAP input file' )

C...........  Determine number of lines in filelist; this will be the maximum
C             number of airport sources
         ALLOCATE( CVHAPT( NLINES,2 ), STAT=IOS )
         CALL CHECKMEM( IOS, 'CVHAPT', PROGNAME )
         ALLOCATE( NPHAPT( NLINES ), STAT=IOS )
         CALL CHECKMEM( IOS, 'NPHAPT', PROGNAME )

         NPHAPT = ' '
         CVHAPT = 0.0
        
         N = 0
         IREC  = 0

         DO I = 1, NLINES

             READ ( TDEV, 93000, IOSTAT=IOS ) LINE
             IREC = IREC + 1

             IF ( IOS .GT. 0 ) THEN
                 WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading HAP_FACTOR '//
     &                'description file at line', IREC
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             END IF

C.............  Left adjust line
             LINE = TRIM( LINE )

C.............  Skip blank and comment lines
             IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
             CALL PARSLINE( LINE, 3, SEGMENT )

             N = N + 1
             NPHAPT( N ) = SEGMENT( 1 )
             CVHAPT( N,1 ) = STR2REAL( SEGMENT( 2 ) )
             CVHAPT( N,2 ) = STR2REAL( SEGMENT( 3 ) )

         END DO    ! end of loop

         IF( N == 0 ) THEN 
             MESG = 'ERROR: No entries of TURBINE HAPs'
             CALL M3MSG2( MESG )
         END IF

         NFAC_TURBINE = N

         CLOSE( TDEV )

         RETURN

C...................  FORMAT  STATEMENTS   ............................

C.......  Formatted file I/O formats...... 93xxx
93000    FORMAT( A )

C.......  Internal buffering formats...... 94xxx
94010    FORMAT( 10 ( A, :, I8, :, 2X  ) )

       END SUBROUTINE READ_TURBINE_HAP_FACTOR

C******************  INTERNAL SUBPROGRAMS  *****************************

C      This subroutine stores a list of aircraft AEDT HAPs to convert 
C      to species mass emission based upon TOG emissions

       SUBROUTINE READ_PISTON_HAP_FACTOR

C...........   Local variables
         INTEGER         I, N                  ! indices and counters

         INTEGER      :: NLINES = 0            ! number of lines in input file

         CHARACTER(256)  LINE                  ! Read buffer for a line
         CHARACTER(300)  MESG                  ! Message buffer
         CHARACTER(16 ):: SEGMENT( 3 ) = ' '   ! line parsing array

C......................................................................

C.........  Get the number of lines
         NLINES = GETFLINE( PDEV, 'PISTON_HAP input file' )

C...........  Determine number of lines in filelist; this will be the maximum
C             number of airport sources
         ALLOCATE( CVHAPP( NLINES,2 ), STAT=IOS )
         CALL CHECKMEM( IOS, 'CVHAPP', PROGNAME )
         ALLOCATE( NPHAPP( NLINES ), STAT=IOS )
         CALL CHECKMEM( IOS, 'NPHAPP', PROGNAME )

         NPHAPP = ' '
         CVHAPP = 0.0
        
         N = 0
         IREC  = 0

         DO I = 1, NLINES

             READ ( PDEV, 93000, IOSTAT=IOS ) LINE
             IREC = IREC + 1

             IF ( IOS .GT. 0 ) THEN
                 WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading HAP_FACTOR '//
     &                'description file at line', IREC
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             END IF

C.............  Left adjust line
             LINE = TRIM( LINE )

C.............  Skip blank and comment lines
             IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
             CALL PARSLINE( LINE, 3, SEGMENT )

             N = N + 1
             NPHAPP( N ) = SEGMENT( 1 )
             CVHAPP( N,1 ) = STR2REAL( SEGMENT( 2 ) )
             CVHAPP( N,2 ) = STR2REAL( SEGMENT( 3 ) )

         END DO    ! end of loop

         IF( N == 0 ) THEN 
             MESG = 'ERROR: No entries of PISTON HAPs'
             CALL M3MSG2( MESG )
         END IF

         NFAC_PISTON = N

         CLOSE( PDEV )

         RETURN

C...................  FORMAT  STATEMENTS   ............................

C.......  Formatted file I/O formats...... 93xxx
93000    FORMAT( A )

C.......  Internal buffering formats...... 94xxx
94010    FORMAT( 10 ( A, :, I8, :, 2X  ) )

       END SUBROUTINE READ_PISTON_HAP_FACTOR

C******************  INTERNAL SUBPROGRAMS  *****************************
       
C      This subroutine stores a list of aircraft AEDT HAPs to convert 
C      to species mass emission based upon TOG emissions

       SUBROUTINE READ_TURBINE_SPECIATION

C..........   Local variables
         INTEGER         I, N                  ! indices and counters

         INTEGER      :: NLINES = 0            ! number of lines in input file

         CHARACTER(256)  LINE                  ! Read buffer for a line
         CHARACTER(300)  MESG                  ! Message buffer
         CHARACTER(16 ):: SEGMENT( 3 ) = ' '   ! line parsing array

C......................................................................

C.........  Get the number of lines
         NLINES = GETFLINE( TSDEV, 'TURBINE_SPC input file' )

C...........  Determine number of lines in filelist; this will be the maximum
C             number of airport sources
         ALLOCATE( CVSPCT( NLINES,2 ), STAT=IOS )
         CALL CHECKMEM( IOS, 'CVSPCT', PROGNAME )
         ALLOCATE(  NSSPCT( NLINES ), STAT=IOS )
         CALL CHECKMEM( IOS, ' NSSPCT', PROGNAME )

          NSSPCT = ' '
         CVSPCT = 0.0
        
         N = 0
         IREC  = 0

         DO I = 1, NLINES

             READ ( TSDEV, 93000, IOSTAT=IOS ) LINE
             IREC = IREC + 1

             IF ( IOS .GT. 0 ) THEN
                 WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading TURBINE_SPC '//
     &                'description file at line', IREC
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             END IF

C.............  Left adjust line
             LINE = TRIM( LINE )

C.............  Skip blank and comment lines
             IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
             CALL PARSLINE( LINE, 3, SEGMENT )

             N = N + 1
              NSSPCT( N ) = SEGMENT( 1 )
             CVSPCT( N,1 ) = STR2REAL( SEGMENT( 2 ) )
             CVSPCT( N,2 ) = STR2REAL( SEGMENT( 3 ) )

         END DO    ! end of loop

         IF( N == 0 ) THEN 
             MESG = 'ERROR: No entries of TURBINE Speciations'
             CALL M3MSG2( MESG )
         END IF

         NSPC_TURBINE = N
         
         CLOSE( TSDEV )

         RETURN

C...................  FORMAT  STATEMENTS   ............................

C.......  Formatted file I/O formats...... 93xxx
93000    FORMAT( A )

C.......  Internal buffering formats...V... 94xxx
94010    FORMAT( 10 ( A, :, I8, :, 2X  ) )

       END SUBROUTINE READ_TURBINE_SPECIATION

C******************  INTERNAL SUBPROGRAMS  *****************************

C      This subroutine stores a list of aircraft AEDT HAPs to convert 
C      to species mass emission based upon TOG emissions

       SUBROUTINE READ_PISTON_SPECIATION

C...........   Local variables
         INTEGER         I, N                  ! indices and counters

         INTEGER      :: NLINES = 0            ! number of lines in input file

         CHARACTER(256)  LINE                  ! Read buffer for a line
         CHARACTER(300)  MESG                  ! Message buffer
         CHARACTER(16 ):: SEGMENT( 3 ) = ' '   ! line parsing array

C......................................................................

C.........  Get the number of lines
         NLINES = GETFLINE( PSDEV, 'PISTON_SPC input file' )

C...........  Determine number of lines in filelist; this will be the maximum
C             number of airport sources
         ALLOCATE( CVSPCP( NLINES,2 ), STAT=IOS )
         CALL CHECKMEM( IOS, 'CVSPCP', PROGNAME )
         ALLOCATE(  NSSPCP( NLINES ), STAT=IOS )
         CALL CHECKMEM( IOS, ' NSSPCP', PROGNAME )

          NSSPCP = ' '
         CVSPCP = 0.0
        
         N = 0
         IREC  = 0

         DO I = 1, NLINES

             READ ( PSDEV, 93000, IOSTAT=IOS ) LINE
             IREC = IREC + 1

             IF ( IOS .GT. 0 ) THEN
                 WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading PISTON_SPC '//
     &                'description file at line', IREC
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             END IF

C.............  Left adjust line
             LINE = TRIM( LINE )

C.............  Skip blank and comment lines
             IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
             CALL PARSLINE( LINE, 3, SEGMENT )

             N = N + 1
              NSSPCP( N ) = SEGMENT( 1 )
             CVSPCP( N,1 ) = STR2REAL( SEGMENT( 2 ) )
             CVSPCP( N,2 ) = STR2REAL( SEGMENT( 3 ) )

         END DO    ! end of loop

         IF( N == 0 ) THEN 
             MESG = 'ERROR: No entries of PISTON Speciations'
             CALL M3MSG2( MESG )
         END IF

         NSPC_PISTON = N
         
         CLOSE( PSDEV )

         RETURN

C...................  FORMAT  STATEMENTS   ............................

C.......  Formatted file I/O formats...... 93xxx
93000    FORMAT( A )

C.......  Internal buffering formats...... 94xxx
94010    FORMAT( 10 ( A, :, I8, :, 2X  ) )

       END SUBROUTINE READ_PISTON_SPECIATION


C......................................................................

       END PROGRAM AEDTPROC
