
        PROGRAM TMPBEIS3

C***********************************************************************
C  program body starts at line  187
C
C  DESCRIPTION:
C       Computes hourly time stepped gridded biogenic emissions using 
C       normalized gridded emissions from NORMBEISv0.9 and postprocessed MM5
C       meteorology.
C
C  PRECONDITIONS REQUIRED:
C       Postprocessed MM5 meteorology that contains temperature, 
C       solar radiation, and pressure data. 
C       Normalized gridded emissions B3GRD from NORMBEIS3v0.9 
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       HRBIO, PREBMET
C
C  REVISION  HISTORY:
C      3/01: Prototype by Jeff Vukovich
C            Tested only on 36km Lambert domain 
C            Summer/winter switch file option not tested
C                  
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C MCNC-Environmental Programs Group
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
C...........   Modules for public variables
C...........   This module contains the speciation profile tables

        USE MODSPRO, ONLY: MXSPEC, SPCNAMES

C...........   This module contains BEIS3 specific arrays
        USE MODBEIS3, ONLY: AVGISOP, AVGMONO, AVGOVOC, AVGNO, AVGLAI

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     !
        INCLUDE 'B3DIMS3.EXT'     ! biogenic-related constants

C...........   PARAMETERS and their descriptions:
C...........   LOCAL PARAMETERS

        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag
        
C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         ENVINT 
        LOGICAL         ENVYN
        CHARACTER*50    GETCFDSC
        INTEGER         GETFLINE
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*16    PROMPTMFILE
        INTEGER         PROMPTFFILE
        INTEGER         TRIMLEN
        CHARACTER*16    VERCHAR

        EXTERNAL        ENVINT, ENVYN, GETFLINE, HHMMSS, INDEX1, PROMPTMFILE,  
     &                  PROMPTFFILE, TRIMLEN, VERCHAR

C.........  Latitude and longitude for zenith angle calculation

        REAL, ALLOCATABLE :: LAT  ( :, : )    !  grid lat (deg) -90 <= LAT <= 90
        REAL, ALLOCATABLE :: LON  ( :, : )    !  grid lon (deg) -180 <= LON <= 180 

C.........  Gridded meteorology data
 
        REAL, ALLOCATABLE :: TASFC ( :, : )     !  level-1 air  temperature (K)
        REAL, ALLOCATABLE :: TSOLAR ( :, :)     !  Photosynthetic Active Radiation (PAR)
        REAL, ALLOCATABLE :: COSZEN( :, : )     !  cosine of zenith angle
        REAL, ALLOCATABLE :: PRES( :, : )       !  surface pressure
        INTEGER, ALLOCATABLE :: SWITCH( :, : )  !  Seasonal switch

C.......   Gridded normalized emissions to use 

        REAL, ALLOCATABLE ::  SISOP( :, : )        !   for isoprene
        REAL, ALLOCATABLE ::  SMONO( :, : )        !   for monoterpenes
        REAL, ALLOCATABLE ::  SOVOC( :, : )        !   for other VOCs
        REAL, ALLOCATABLE ::  SNO( :, : )          !   for nitric oxide
        REAL, ALLOCATABLE ::  SLAI( :, : )         !   for leaf area

C.......   Mole and mass factors
 
        REAL, ALLOCATABLE ::  MLFAC( :, : )           !  mole factors 
        REAL, ALLOCATABLE ::  MSFAC( :, : )           !  mass factors (tons/hour)

C.......   BEIS3 internal, output species

        REAL, ALLOCATABLE :: EMPOL( :, :, : )         ! emissions of biogenic categories
        REAL, ALLOCATABLE :: EMISL( :, :, : )         ! emissions in moles/hour
        REAL, ALLOCATABLE :: EMISS( :, :, : )         ! emissions in tons/hour

        CHARACTER*5      CTZONE     ! string of time zone
        CHARACTER*4      BTMP       ! temporary variable string 
        CHARACTER*16     RADNAM     ! string for shortwave radiation reaching ground
        CHARACTER*16     TMPRNAM    ! string for temperature 
        CHARACTER*16     VTMP       ! temporary variable string
        CHARACTER*50  :: METSCEN    !  temporary string for met scenario name
        CHARACTER*50  :: CLOUDSHM   !  temporary string for cloud scheme name
        CHARACTER*50  :: LUSE       !  temporary string for land use description


        CHARACTER*16,ALLOCATABLE ::  EMSPC( : )   ! names of emitting species 
        CHARACTER(LEN=SPNLEN3)       SPPRO        ! speciation profile to use

        CHARACTER*80    PARMENU( 2 )            ! Methods to calc. PAR
        DATA     PARMENU
     &           / 'Use MM5 generated RGRND',  
     &             'Assume Clear Skies'     /

C...........   Logical names and unit numbers

        INTEGER         LDEV    !  unit number for log device
        INTEGER         RDEV    !  unit number for speciation profiles file
            
        CHARACTER*16    ENAME   !  logical name for emissions output (moles)
        CHARACTER*16    SNAME   !  logical name for emissions output (mass)
        CHARACTER*16    NNAME   !  logical name for normalized-emissions input
        CHARACTER*16    GNAME   !  logical name for GRID_CRO_2D
        CHARACTER*16    BNAME   !  logical name for frost switch input
        CHARACTER*16    M3NAME  !  logical name for MET_FILE1
        CHARACTER*16    M2NAME  !  logical name for MET_FILE2

C...........   Other variables and their descriptions:

        INTEGER         B, M    !  counters for biogenic, model species
        INTEGER         I, J, K, L, C, R  !  loop counters and subscripts
        INTEGER         HR      !  current simulation hour
        INTEGER         IOS     !  temporay IO status
        INTEGER         JDATE   !  current simulation date (YYYYDDD)
        INTEGER         JTIME   !  current simulation time (HHMMSS)
        INTEGER         LDATE   !  previous simulation date
        INTEGER         MDATE   !  met file 1 start date
        INTEGER         MSPCS   ! no. of emitting species
        INTEGER         MTIME   !  met file 1 start time
        INTEGER         MXSTEPS !  maximum number of time steps
        INTEGER         NCOLS   ! no. of grid columns
        INTEGER         NGRID   ! no. of grid cells
        INTEGER         NLINES  ! no. of lines in GSPRO speciation profiles file  
        INTEGER         NROWS   ! no. of grid rows
        INTEGER         NSTEPS  !  duration of met file
        INTEGER         BSTEPS  ! no. of hourly time steps for output
        INTEGER         PARTYPE !  method number to calculate PAR
        INTEGER         RDATE   !  met file 2 start date 
        INTEGER         RTIME   !  met file 2 start time
        INTEGER         TSTEP   !  time step set by environment
        INTEGER         TZONE   !  output-file time zone ; not used in program

        LOGICAL         EFLAG   !  error flag
        LOGICAL ::      SAMEFILE = .TRUE.   ! radiation/cld and tmpr data in same file 
        LOGICAL ::      SWITCH_FILE = .TRUE.  ! use frost switch file
        LOGICAL ::      ASSUME_SUMMER = .TRUE. ! use summer normalized emissions
        LOGICAL ::      GETATN

        CHARACTER*300   MESG    !  message buffer for M3EXIT()

        CHARACTER*16 :: PROGNAME = 'TMPBEIS3'   !  program name

C***********************************************************************
C   begin body of program TMPBEIS3

        LDEV = INIT3()
 
C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Evaluate the environment variables...

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.........  Write time zone to character string
        WRITE( CTZONE,94000 ) TZONE

        MESG = 'Speciation profile to use for biogenics'
        CALL ENVSTR( 'BIOG_SPRO', MESG, 'B3V09', SPPRO, IOS )

C........  Are the rad/cld data in the same file as tmpr data

        MESG = 'Radiation/cloud in same file as temperature data?'
        SAMEFILE = ENVYN ( 'BIOMET_SAME', MESG, .FALSE., IOS )

        RDEV = PROMPTFFILE(
     &           'Enter logical name for SPECIATION PROFILES file',
     &           .TRUE., .TRUE., 'GSPRO', PROGNAME )


C.........  Scan speciation profiles file to get all of the pollutant-species
C           combinations that are valid for the pollutants in the inventory.
C.........  The species names are sorted in ABC order for each pollutant, and
C           and the pollutants are in the same order as BIOSPC.
C.........  Also retrieve the maximum number of species per pollutant and
C           maximum number of profile entries per pollutant.

        CALL DSCSPROF( RDEV, BSPCS, BIOSPC )

        NLINES = GETFLINE( RDEV, 'Speciation profile file' )

        ALLOCATE( EMSPC( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMSPC', PROGNAME )

        MSPCS = 0

C.......... Find emitting species names

        DO I = 1, BSPCS
          DO J = 1, MXSPEC
            IF ( SPCNAMES( J, I ) .EQ. ' ' ) CYCLE

            K = INDEX1 ( SPCNAMES( J, I ), MSPCS, EMSPC ) 

            IF ( K .EQ. 0 ) THEN
               MSPCS = MSPCS + 1
               EMSPC( MSPCS ) = SPCNAMES( J, I )
            ENDIF
          ENDDO
        ENDDO

        ALLOCATE( MLFAC ( MSPCS, BSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MLFAC', PROGNAME )

        ALLOCATE( MSFAC ( MSPCS, BSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MSFAC', PROGNAME )
 
 
C.............  Read speciation profiles file

        MESG = BLANK5 // 'Reading biogenic speciation profile...'
        CALL M3MSG2( MESG )

        CALL RDBPRO( RDEV, SPPRO, BSPCS, BIOSPC, MSPCS, EMSPC, 
     &               MLFAC, MSFAC ) 

C.......   Get the method to calculate PAR

        PARTYPE = ENVINT( 'BG_CLOUD_TYPE', 
     &                    'How PAR will be calculated', 1, IOS )

        MESG =  'PAR calculation will ' // PARMENU( PARTYPE )  
        WRITE( LDEV, 92000 ) MESG

C.......   Check to see if frost date switch file to be used

        MESG = 'Using a frost date switch file?'
        SWITCH_FILE = ENVYN ( 'BIOSW_YN', MESG, .TRUE., IOS )

C.......   Get normalized emissions file, BGRD

        NNAME = PROMPTMFILE( 
     &          'Enter name for NORMALIZED EMISSIONS input file',
     &          FSREAD3, 'B3GRD', PROGNAME )

C......    Read description of normalized emissions file

        IF ( .NOT. DESC3( NNAME ) ) THEN

            MESG = 'Could not get description of file "' //
     &             NNAME( 1:TRIMLEN( NNAME ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

C.......   Initialize grid definition 

        CALL CHKGRID( NNAME, 'GRID' , 0, EFLAG )

        NCOLS = NCOLS3D 
        NROWS = NROWS3D
        NGRID = NCOLS3D * NROWS3D
        LUSE  = GETCFDSC( FDESC3D, '/LANDUSE/', .FALSE. )

        IF ( SWITCH_FILE ) THEN

C.......    Get bioseason switch file, BIOSEASON

           BNAME = PROMPTMFILE( 
     &          'Enter name for season switch input file',
     &          FSREAD3, 'BIOSEASON', PROGNAME )
           
C......    Read description of switch file

           IF ( .NOT. DESC3( BNAME ) ) THEN

              MESG = 'Could not get description of file "' //
     &             NNAME( 1:TRIMLEN( BNAME ) ) // '"'
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

           ENDIF

C.......   Check grid definition 

           CALL CHKGRID( BNAME, 'GRID' , 0 , EFLAG )

           IF ( EFLAG ) THEN
              MESG = 'Problems opening input files. See ERROR(S) above.'
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
           END IF

           ALLOCATE( SWITCH( NCOLS, NROWS), STAT=IOS )
           CALL CHECKMEM( IOS, 'SWITCH', PROGNAME )
           SWITCH = 0   ! array

        ELSE

           MESG = 'Use summer normalized emissions?'
           ASSUME_SUMMER = ENVYN ( 'SUMMER_YN', MESG, .TRUE., IOS )

        ENDIF

C.......   Get temperature file

        M3NAME = PROMPTMFILE( 
     &          'Enter name for gridded temperature input file',
     &          FSREAD3, 'MET_FILE1', PROGNAME )

C......    Get name of temperature variable to use

        MESG = 'Variable name for temperature'
        CALL ENVSTR( 'TMPR_VAR', MESG, 'TA', TMPRNAM, IOS )

C......    Read description of temperature file 
        
        IF ( .NOT. DESC3( M3NAME ) ) THEN

            MESG = 'Could not get description of file "' //
     &             M3NAME( 1:TRIMLEN( M3NAME ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

C......  Make sure requested tmpr variable available

        J = INDEX1( TMPRNAM , NVARS3D, VNAME3D )

        IF ( J .LE. 0 ) THEN

             MESG = 'Could not find ' // TMPRNAM // 'in file ' //
     &                   M3NAME( 1:TRIMLEN( M2NAME ) ) 

             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
 
        ENDIF

        CALL CHKGRID( M3NAME, 'GRID' , 0 , EFLAG ) 

C........  If grid definition does not match BGRD file then stop

        IF ( SAMEFILE .AND. EFLAG ) THEN
           MESG = 'Problems opening input files. See ERROR(S) above.'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C........  Get met and cloud scheme descriptions from M3NAME if they exist

        METSCEN  = GETCFDSC( FDESC3D, '/MET SCENARIO/', .FALSE. )
        CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. )

        IF ( PARTYPE .NE. 2 ) THEN

C.........  Open second met file if needed

           IF ( .NOT. SAMEFILE ) THEN

             M2NAME = PROMPTMFILE(
     &           'Enter name for gridded radiation/cloud input file',
     &           FSREAD3, 'MET_FILE2', PROGNAME )

C......    Read description of radiation/cloud file

             IF ( .NOT. DESC3( M2NAME ) ) THEN

                MESG = 'Could not get description of file "' //
     &             M2NAME( 1:TRIMLEN( M2NAME ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

             ENDIF

             CALL CHKGRID( M2NAME, 'GRID' , 0 , EFLAG )

C........  If grid definition does not match BGRD file then stop

             IF ( EFLAG ) THEN
              MESG = 'Problems opening input files. See ERROR(S) above.'
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             END IF

             CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. )
        
           ELSE
             
             M2NAME = M3NAME
 
           ENDIF

C......    Get name of radiation variable to use

           MESG = 'Variable name for radiation'
           CALL ENVSTR( 'RAD_VAR', MESG, 'RGRND', RADNAM, IOS )


           J = INDEX1( RADNAM, NVARS3D, VNAME3D )
           IF ( J .LE. 0 ) THEN

              MESG = 'Could not find ' // RADNAM // 'in file ' //
     &               M2NAME( 1:TRIMLEN( M2NAME ) )

              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

           ENDIF

        ENDIF

C.......   Get default time characteristic for output file:
C.......   If we're going to prompt, then set the defaults based on met
C.......      otherwise, use environment variables to set defaults

        JDATE  = SDATE3D
        JTIME  = STIME3D
        NSTEPS = MXREC3D

        CALL GETM3EPI( TZONE, JDATE, JTIME, TSTEP, NSTEPS )

C.......   Build description for, and create/open output file
C.......   (all but variables-table in description is borrowed from M3NAME)

        SDATE3D = JDATE
        STIME3D = JTIME
        BSTEPS  = NSTEPS
        MXREC3D = BSTEPS
        NVARS3D = MSPCS
        NLAYS3D = 1
        TSTEP3D = 10000   ! only 1-hour timestep supported

        DO  M = 1, MSPCS
            VNAME3D( M ) = EMSPC( M )
            UNITS3D( M ) = 'moles/hr'
            VDESC3D( M ) = 'biogenic emissions of the indicated species'
            VTYPE3D( M ) = M3REAL
        ENDDO

        FDESC3D = ' '   ! array

        FDESC3D( 1 ) = 'Gridded biogenic emissions from SMOKE-BEIS3'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        FDESC3D( 4 ) = '/TZONE/ '   // CTZONE
        FDESC3D( 5 ) = '/LANDUSE/ ' // LUSE
        FDESC3D( 6 ) = '/MET SCENARIO/ ' // METSCEN
        FDESC3D( 7 ) = '/CLOUD SCHEME/ ' // CLOUDSHM

C......  Open first output file (moles/hour)

        ENAME = PROMPTMFILE(
     &          'Enter name for B3GTS output file - moles',
     &          FSUNKN3, 'B3GTS_L', PROGNAME )

        DO M = 1, MSPCS
            UNITS3D( M ) = 'tons/hr'
        ENDDO

C....... Open second output file (tons/hour)

        SNAME = PROMPTMFILE(
     &          'Enter name for B3GTS output file - mass',
     &          FSUNKN3, 'B3GTS_S', PROGNAME )

C.......   Build name table for variables in normalized emissions file

        ALLOCATE( AVGISOP( NCOLS, NROWS, NSEASONS), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGISOP', PROGNAME )
        ALLOCATE( AVGMONO( NCOLS, NROWS, NSEASONS), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGMONO', PROGNAME )
        ALLOCATE( AVGOVOC( NCOLS, NROWS, NSEASONS), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGOVOC', PROGNAME )
        ALLOCATE( AVGNO( NCOLS, NROWS, NSEASONS), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGNO', PROGNAME )
        ALLOCATE( AVGLAI( NCOLS, NROWS, NSEASONS), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGLAI', PROGNAME )

        AVGISOP = 0.0 ! array
        AVGMONO = 0.0 ! array
        AVGOVOC = 0.0 ! array
        AVGNO   = 0.0 ! array
        AVGLAI  = 0.0 ! array

C........  Open 2-D grid parameters file to get LAT and LON

        GNAME = PROMPTMFILE( 
     &          'Enter name for 2D GRID PARAMETERS input file',
     &          FSREAD3, 'GRID_CRO_2D', 'TMPBIO' )

        IF ( .NOT. DESC3( GNAME ) ) THEN
             CALL M3EXIT( PROGNAME, 0, 0,
     &       'Could not get description of file "'
     &       // GNAME( 1:TRIMLEN( GNAME ) ) // '"', 2 )

        ENDIF
 
        CALL CHKGRID( GNAME, 'GRID' , 0 , EFLAG ) 

C........  If grid definition does not match BGRD file then stop

        IF ( EFLAG ) THEN
             MESG = 'Problems opening input files. See ERROR(S) above.'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        ALLOCATE( LAT ( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LAT', PROGNAME )

        ALLOCATE( LON ( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LON', PROGNAME )

        ALLOCATE( COSZEN ( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'COSZEN', PROGNAME )

        IF ( .NOT. READ3( GNAME, 'LAT', 1, 0, 0, LAT ) ) THEN
             MESG = 'Could not read LAT from file "' //
     &                GNAME( 1:TRIMLEN( GNAME ) ) // '"'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. READ3( GNAME, 'LON', 1, 0, 0, LON ) ) THEN
             MESG = 'Could not read LON from file "' //
     &                GNAME( 1:TRIMLEN( GNAME ) ) // '"'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.......   Loops reading the various categories of normalized emissions:
 
        DO  M = 1, NSEASONS 
 
            I = 1
            BTMP = BIOTYPES( I ) 

            VTMP = 'AVG' // BTMP( 1: TRIMLEN( BTMP ))
     &                     // SEASON( M ) 
        
            IF ( .NOT. READ3( NNAME, VTMP, 1, 0, 0, 
     &                      AVGISOP( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  VTMP( 1 : TRIMLEN( VTMP ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            I = I + 1

            BTMP = BIOTYPES( I ) 

            VTMP = 'AVG' // BTMP( 1: TRIMLEN( BTMP ))
     &                     // SEASON( M ) 
        
            IF ( .NOT. READ3( NNAME, VTMP, 1, 0, 0, 
     &                      AVGMONO( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  VTMP( 1 : TRIMLEN( VTMP ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            I = I + 1

            BTMP = BIOTYPES( I ) 

            VTMP = 'AVG' // BTMP( 1: TRIMLEN( BTMP ))
     &                     // SEASON( M ) 
        
            IF ( .NOT. READ3( NNAME, VTMP, 1, 0, 0, 
     &                      AVGOVOC( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  VTMP( 1 : TRIMLEN( VTMP ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            I = I + 1

            BTMP = BIOTYPES( I ) 

            VTMP = 'AVG' // BTMP( 1: TRIMLEN( BTMP ))
     &                     // SEASON( M ) 
        
            IF ( .NOT. READ3( NNAME, VTMP, 1, 0, 0, 
     &                      AVGNO( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  VTMP( 1 : TRIMLEN( VTMP ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            I = I + 1

            BTMP = BIOTYPES( I ) 

            VTMP = 'AVG' // BTMP( 1: TRIMLEN( BTMP ))
     &                     // SEASON( M ) 
        
            IF ( .NOT. READ3( NNAME, VTMP, 1, 0, 0, 
     &                      AVGLAI( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  VTMP( 1 : TRIMLEN( VTMP ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        ENDDO        !  end loop to read normalized emissions 

C.......   Allocate memory for met and emissions 

        ALLOCATE( TASFC( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TASFC', PROGNAME )
        ALLOCATE( TSOLAR( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TSOLAR', PROGNAME )

        ALLOCATE( EMPOL( NCOLS, NROWS, BSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMPOL', PROGNAME )

        ALLOCATE( PRES( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRES', PROGNAME )

        ALLOCATE( EMISL( NCOLS, NROWS, MSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISL', PROGNAME )
        ALLOCATE( EMISS( NCOLS, NROWS, MSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISS', PROGNAME )

        ALLOCATE( SISOP( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SISOP', PROGNAME )
        ALLOCATE( SMONO( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SMONO', PROGNAME )
        ALLOCATE( SOVOC( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SOVOC', PROGNAME )
        ALLOCATE( SNO( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SNO', PROGNAME )
        ALLOCATE( SLAI( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SLAI', PROGNAME )

C.........  Set up gridded met file(s) dates and times for specific time zone

        IF( M3NAME .NE. ' ' ) THEN
          CALL PREBMET( M3NAME, M2NAME, SAMEFILE, TZONE, TSTEP3D, JDATE, 
     &                  JTIME, NSTEPS, MDATE, MTIME, RDATE, RTIME )
        END IF

        IF ( M2NAME .EQ. M3NAME ) THEN

           RDATE = MDATE  
           RTIME = MTIME
  
        ENDIF 

C.......... Initialize normalized emissons to be used 

        IF ( ASSUME_SUMMER ) THEN

          SISOP = AVGISOP( 1:NCOLS, 1:NROWS, NSUMMER )
          SMONO = AVGMONO( 1:NCOLS, 1:NROWS, NSUMMER )
          SOVOC = AVGOVOC( 1:NCOLS, 1:NROWS, NSUMMER )
          SNO   = AVGNO( 1:NCOLS, 1:NROWS, NSUMMER )
          SLAI  = AVGLAI( 1:NCOLS, 1:NROWS, NSUMMER )

        ELSE

          SISOP = AVGISOP( 1:NCOLS, 1:NROWS, NWINTER )
          SMONO = AVGMONO( 1:NCOLS, 1:NROWS, NWINTER )
          SOVOC = AVGOVOC( 1:NCOLS, 1:NROWS, NWINTER )
          SNO   = AVGNO( 1:NCOLS, 1:NROWS, NWINTER )
          SLAI  = AVGLAI( 1:NCOLS, 1:NROWS, NWINTER )

        ENDIF

        LDATE = 0     

C.........  Loop thru the number of time steps (hourly)
 
        DO  HR = 1, BSTEPS

           EMISL = 0   !  array
           EMISS = 0   !  array
           EMPOL = 0   !  array

           IF( JDATE .NE. LDATE ) THEN

               CALL WRDAYMSG( JDATE, MESG )               

C..........  If new date, read season switch 

               IF ( SWITCH_FILE ) THEN
                 
                  MESG = 'Reading gridded season switch data..'  
                  CALL M3MESG( MESG ) 
                  IF ( .NOT. READ3( BNAME, 'SEASON', 1, 
     &                JDATE, 0, SWITCH ) ) THEN
                     MESG = 'Could not read SEASON from file ' //
     &               BNAME( 1:TRIMLEN( BNAME ) )
                     CALL M3EXIT( PROGNAME, JDATE, 0, MESG, 2 )
                  END IF

                  MESG = 'Applying gridded season switch data..' 
                  CALL M3MESG( MESG )

                  DO I = 1, NCOLS
                    DO J = 1, NROWS

C..........   If switch equal to 0 use winter normalized emissions

                       IF ( SWITCH ( I, J ) .EQ. 0 ) THEN
                          
                            SISOP( I, J ) = AVGISOP( I, J, 2 ) 
                            SMONO( I, J ) = AVGMONO( I, J, 2 )
                            SOVOC( I, J ) = AVGOVOC( I, J, 2 )
                            SNO( I, J )   = AVGNO( I, J, 2 )
                            SLAI( I, J )  = AVGLAI( I, J, 2 )
    
                       ELSE

                            SISOP( I, J ) = AVGISOP( I, J, 1 ) 
                            SMONO( I, J ) = AVGMONO( I, J, 1 )
                            SOVOC( I, J ) = AVGOVOC( I, J, 1 )
                            SNO( I, J )   = AVGNO( I, J, 1 )
                            SLAI( I, J )  = AVGLAI( I, J, 1 )
 
                       ENDIF

                    ENDDO
                  ENDDO
 
               ENDIF  ! switch file endif

           ENDIF   ! new day endif

C.............  Write to screen because WRITE3 only writes to LDEV
           WRITE( *, 94030 ) HHMMSS( JTIME )

C............. Compute zenith angle

           CALL CZANGLE( JDATE, JTIME, NCOLS, NROWS, LAT, LON, 
     &                   COSZEN, GETATN ) 

C.............  Read temperature data

           IF ( .NOT. READ3( M3NAME, TMPRNAM, 1, 
     &          MDATE, MTIME, TASFC ) ) THEN
              MESG = 'Could not read ' // TMPRNAM // 'from file ' //
     &                M3NAME( 1:TRIMLEN( M3NAME ) )
              CALL M3EXIT( PROGNAME, MDATE, MTIME, MESG, 2 )
           END IF

C..............  Read surface pressure data 

           IF ( .NOT. READ3( M3NAME, 'PRSFC', 1, 
     &                        MDATE, MTIME, PRES ) ) THEN
                  MESG = 'Could not read PRES from file "' //
     &                   M3NAME( 1:TRIMLEN( M3NAME ) ) // '"'
                  CALL M3EXIT( PROGNAME, MDATE, MTIME, MESG, 2 )
           END IF

C...............  convert from Pa to millibars

           DO  C = 1, NCOLS
               DO  R = 1, NROWS
                   PRES( C, R ) = PRES( C, R ) * 0.010  ! Pa to mb
               ENDDO
           ENDDO

C............. If necessary read solar radiation

           IF ( PARTYPE .EQ. 1 ) THEN      ! Use RADNAM

              IF ( .NOT. READ3( M2NAME, RADNAM, ALLAYS3, RDATE,
     &             RTIME, TSOLAR(1,1) ) ) THEN
                 MESG = 'Could not read ' // RADNAM // 'from file ' //
     &                M2NAME( 1:TRIMLEN( M2NAME ) ) 

                 CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )
              END IF

C............. Calculate non-speciated emissions
C............. must pass met date and time here

  
              CALL HRBEIS3( MDATE, MTIME, NCOLS, NROWS, COSZEN,
     &                   SISOP, SMONO, SOVOC, SNO, SLAI,
     &                   TASFC, TSOLAR, PRES, EMPOL )
  

           ELSE  IF ( PARTYPE .EQ. 2 ) THEN     ! Use clouds or clear skies to calculate PAR

C
C....... Option of clear skies not available yet.
C

           END IF 

C............. Speciate emissions

           DO I = 1, NCOLS
             DO J = 1, NROWS
               DO L = 1, MSPCS
                 DO K = 1, BSPCS
                   EMISL( I, J, L  ) = EMISL( I ,J, L ) +
     &                             EMPOL( I, J, K ) * MLFAC( L, K )
                   EMISS( I, J, L  ) = EMISS( I, J, L ) +
     &                             EMPOL( I, J, K ) * MSFAC( L, K ) 
                 ENDDO
               ENDDO
             ENDDO
           END DO

C.............  Write out speciated emissions    

           IF ( .NOT. WRITE3( ENAME, 'ALL', 
     &                        JDATE, JTIME, EMISL ) ) THEN
               CALL M3EXIT( PROGNAME, JDATE, JTIME, 
     &                      'Error writing BEIS3 OUTPUT file' , 2 )
           END IF                              !  if write3 failed

           IF ( .NOT. WRITE3( SNAME, 'ALL',
     &                        JDATE, JTIME, EMISS ) ) THEN
               CALL M3EXIT( PROGNAME, JDATE, JTIME,
     &                      'Error writing BEIS3 OUTPUT file' , 2 )
           END IF                              !  if write3 failed

C.............. Next time step

           LDATE = JDATE
           CALL NEXTIME( JDATE, JTIME, 10000 )

           CALL NEXTIME( MDATE, MTIME, 10000 ) 

           IF ( M2NAME .EQ. M3NAME ) THEN

              RDATE = MDATE
              RTIME = MTIME

           ELSE

             CALL NEXTIME( RDATE, RTIME, 10000 ) 

           ENDIF
      ENDDO                !  end loop on hours HR

C.........   End of program:

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )

94030   FORMAT( 8X, 'at time ', A8 )

        END PROGRAM TMPBEIS3  

