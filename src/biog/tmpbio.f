
        PROGRAM TMPBIO

C***********************************************************************
C  program body starts at line 219
C
C  DESCRIPTION:
C       Computes time stepped gridded biogenic emissions in terms of 
C       normalized gridded emissions from RAWBIO and postprocessed MM5
C       meteorology.
C
C  PRECONDITIONS REQUIRED:
C       Postprocessed MM5 meteorology that contains temperature and/or
C       solar radiation/cloud data. 
C       Normalized gridded emissions BGRD from RAWBIO or GRDBIO
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       HRBIO, PREBMET
C
C  REVISION  HISTORY:
C       Prototype 11/99 by JMV from version 4.9 of TMPBIO SMOKE prototype
C                  
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************

C...........   Modules for public variables
C...........   This module contains the speciation profile tables
        USE MODSPRO

C...........   This module contains the global variables for the 3-d grid
        USE MODGRID

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     !
        INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants


C...........   PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag
     
C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF   
        INTEGER         ENVINT 
        LOGICAL         ENVYN
        CHARACTER*50    GETCFDSC
        INTEGER         GETDATE
        INTEGER         GETFLINE
        INTEGER         GETNUM
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*16    PROMPTMFILE
        INTEGER         PROMPTFFILE
        CHARACTER*16    VERCHAR

        EXTERNAL        CRLF, ENVINT, ENVYN, GETDATE, GETFLINE, GETNUM, 
     &                  GETYN, HHMMSS, INDEX1, PROMPTMFILE, 
     &                  PROMPTFFILE, VERCHAR

C.........  Gridded meteorology data
                
        REAL, ALLOCATABLE :: LAT  ( :, : )    !  grid lat (deg) -90 <= LAT <= 90
        REAL, ALLOCATABLE :: LON  ( :, : )    !  grid lon (deg) -180 <= LON <= 180 
        REAL, ALLOCATABLE :: TASFC ( :, : )    !  level-1 air  temperature (K)
        REAL, ALLOCATABLE :: PRSFC  ( :, : )    !  pressure (Pa)
        REAL, ALLOCATABLE :: TSOLAR ( :, :)     !  Photosynthetic Active Radiation (PAR)
        INTEGER, ALLOCATABLE :: SEASON( :, : )  !  Seasonal switch

C.......   Gridded normalized emissions to use in hrbio(s)

        REAL, ALLOCATABLE ::  PINE( :, :, : )         !   for pine
        REAL, ALLOCATABLE ::  DECD( :, :, : )         !   for deciduous forest
        REAL, ALLOCATABLE ::  CONF( :, :, : )         !   for coniferous forest
        REAL, ALLOCATABLE ::  AGRC( :, :, : )         !   for agricultural land
        REAL, ALLOCATABLE ::  LEAF( :, :, : )         !   for leaf area
        REAL, ALLOCATABLE ::  OTHR( :, :, : )         !   for other land

        REAL, ALLOCATABLE ::  AVLAI( :, : )           !  average LAI
        REAL, ALLOCATABLE ::  NORNO( :, :, : )        !  normalized NO emissions

C.......   Gridded winter normalized emissions from BGRDW file 

        REAL, ALLOCATABLE ::  PINEW( :, :, : )     !   for pine
        REAL, ALLOCATABLE ::  DECDW( :, :, : )     !   for deciduous forest
        REAL, ALLOCATABLE ::  CONFW( :, :, : )     !   for coniferous forest
        REAL, ALLOCATABLE ::  AGRCW( :, :, : )     !   for agricultural land
        REAL, ALLOCATABLE ::  LEAFW( :, :, : )     !   for leaf area
        REAL, ALLOCATABLE ::  OTHRW( :, :, : )     !   for other land

        REAL, ALLOCATABLE ::  AVLAIW( :, : )      !  average LAI
        REAL, ALLOCATABLE ::  NORNOW( :, :, :)    !  normalized NO emissions

C........  normalized emissions after seasonal switch 

        REAL, ALLOCATABLE ::  PINES( :, :, : )     !   for pine
        REAL, ALLOCATABLE ::  DECDS( :, :, : )     !   for deciduous forest
        REAL, ALLOCATABLE ::  CONFS( :, :, : )     !   for coniferous forest
        REAL, ALLOCATABLE ::  AGRCS( :, :, : )     !   for agricultural land
        REAL, ALLOCATABLE ::  LEAFS( :, :, : )     !   for leaf area
        REAL, ALLOCATABLE ::  OTHRS( :, :, : )     !   for other land

        REAL, ALLOCATABLE ::  AVLAIS( :, : )      !  average LAI
        REAL, ALLOCATABLE ::  NORNOS( :, :, :)    !  normalized NO emissions

 
        REAL, ALLOCATABLE ::  MLFAC( :, : )           !  mole factors 
        REAL, ALLOCATABLE ::  MSFAC( :, : )           !  mass factors (tons/hour)

C.......   BEIS2 internal, output species

        REAL, ALLOCATABLE :: EMPOL( :, :, : )         ! emissions of biogenic categories
        REAL, ALLOCATABLE :: EMISL( :, :, : )         ! emissions in moles/hour
        REAL, ALLOCATABLE :: EMISS( :, :, : )         ! emissions in tons/hour


        CHARACTER*5      CTZONE     ! string of time zone
        CHARACTER*16     RADNAM     !  string for shortwave radiation reaching ground
        CHARACTER*16     TMPRNAM    !  string for temperature 
        CHARACTER*50  :: METSCEN   !  temporary string for met scenario name
        CHARACTER*50  :: CLOUDSHM  !  temporary string for cloud scheme name
        CHARACTER*50  :: LUSE      !  temporary string for land use description
        CHARACTER*50  :: LUSE2     !  temporary string for 2nd land use desc.

C.......   Name tables for file NNAME

        CHARACTER*16,ALLOCATABLE ::  NORMV( : )   ! names for VOC vbles
        CHARACTER*16,ALLOCATABLE ::  NORMN( : )   ! names for  NO-emission vbles
        CHARACTER*16,ALLOCATABLE ::  EMSPC( : )   ! names of emitting species 

        CHARACTER(LEN=SPNLEN3)       SPPRO        ! speciation profile to use

        CHARACTER*72    PARMENU( 5 )            ! Methods to calc. PAR
        DATA     PARMENU
     &           / 'Use MM5 generated RGND or RSD',
     &             'Use KUO cloud attenuation',
     &             'Use KF  cloud attenuation',   
     &             'Use No deep convection param'  ,
     &             'Assume Clear Skies' /

C...........   Logical names and unit numbers

        INTEGER         LDEV    !  unit number for log device
        INTEGER         RDEV    !  unit number for speciation profiles file
            
        CHARACTER*16    ENAME   !  logical name for emissions output (moles)
        CHARACTER*16    SNAME   !  logical name for emissions output (mass)
        CHARACTER*16    NNAME   !  logical name for normalized-emissions input
        CHARACTER*16    NNAME2  !  logical name for 2nd norm emissions input
        CHARACTER*16    BNAME   !  logical name for frost switch input
        CHARACTER*16    GNAME   !  logical name for GRID_CRO_2D
        CHARACTER*16    M3NAME  !  logical name for MET_FILE1
        CHARACTER*16    M2NAME  !  logical name for MET_FILE2

C...........   Other variables and their descriptions:

        INTEGER         B, M    !  counters for biogenic, model species
        INTEGER         I, II, III, J, JJ, JJJ, K, L, C, R  !  loop counters and subscripts
        INTEGER         HR      !  current simulation hour

        INTEGER         IOS     !  temporay IO status
        INTEGER         JDATE   !  current simulation date (YYYYDDD)
        INTEGER         JTIME   !  current simulation time (HHMMSS)
        INTEGER         LDATE   !  previous simulation date
        INTEGER         MDATE   !  met file 1 start date
        INTEGER         METNCOLS! no. met file columns
        INTEGER         METNROWS! no. met file rows
        INTEGER         MSPCS   ! no. of emitting species
        INTEGER         MTIME   !  met file 1 start time
        INTEGER         MXSTEPS !  maximum number of time steps
        INTEGER         NLINES  ! no. of lines in GSPRO speciation profiles file  
        INTEGER         NSTEPS  !  duration of met file
        INTEGER         BSTEPS  ! no. of hourly time steps for output
        INTEGER         PARTYPE !  method number to calculate PAR
        INTEGER         RDATE   !  met file 2 start date 
        INTEGER         RTIME   !  met file 2 start time
        INTEGER      :: SWNCOLS = 0   !  bioseason 
        INTEGER      :: SWNROWS = 0   !  bioseason 
        INTEGER      :: SWXOFF  = 0   !  bioseason x offset from met grid
        INTEGER      :: SWYOFF  = 0   !  bioseason y offset from met grid
        INTEGER         TZONE   !  output-file time zone ; not used in program

        LOGICAL ::      EFLAG    = .FALSE.  ! error flag
        LOGICAL ::      SAMEFILE = .TRUE.   ! radiation/cld and tmpr data in same file 
        LOGICAL ::      SWITCH_FILE = .TRUE.  ! use frost switch file

        CHARACTER*300   MESG    !  message buffer for M3EXIT()

C.......   Input met and grid variables:
        CHARACTER*16 :: PROGNAME = 'TMPBIO'   !  program name

C***********************************************************************
C   begin body of program TMPBIO

        LDEV = INIT3()
 
C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Evaluate the environment variables...

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.......   Check to see if frost date switch file to be used
        MESG = 'Using a frost date switch file?'
        SWITCH_FILE = ENVYN ( 'BIOSW_YN', MESG, .TRUE., IOS )

C.........  Write time zone to character string
        WRITE( CTZONE,94000 ) TZONE

        MESG = 'Speciation profile to use for biogenics'
        CALL ENVSTR( 'BIOG_SPRO', MESG, '0000', SPPRO, IOS )

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
            END IF
          END DO
        END DO

        ALLOCATE( MLFAC ( MSPCS, BSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LAT', PROGNAME )

        ALLOCATE( MSFAC ( MSPCS, BSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LAT', PROGNAME )
 
 
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

C.......   Get temperature file

        M3NAME = PROMPTMFILE( 
     &          'Enter name for gridded temperature input file',
     &          FSREAD3, 'MET_FILE1', 'TMPBIO' )

C......    Get name of temperature variable to use

        MESG = 'Variable name for temperature'
        CALL ENVSTR( 'TMPR_VAR', MESG, 'TA', TMPRNAM, IOS )

C......    Read description of temperature file 
        
        IF ( .NOT. DESC3( M3NAME ) ) THEN

            MESG = 'Could not get description of file "' //
     &             M3NAME( 1:LEN_TRIM( M3NAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

        END IF

        J = INDEX1( TMPRNAM , NVARS3D, VNAME3D )

        IF ( J .LE. 0 ) THEN

             MESG = 'Could not find ' // TMPRNAM // 'in file ' //
     &                   M3NAME( 1:LEN_TRIM( M2NAME ) ) 

             CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
 
        END IF

C.........  Initialize reference grid with met file
        CALL CHKGRID( M3NAME, 'GRID' , 0 , EFLAG ) 

C.........  Store met grid settings
        METNCOLS = NCOLS
        METNROWS = NROWS

C.........  Get met and cloud scheme descriptions from M3NAME if they exist

        METSCEN  = GETCFDSC( FDESC3D, '/MET SCENARIO/', .FALSE. )
        CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. )

        IF ( PARTYPE .NE. 5 ) THEN

C.........  Open second met file if needed

           IF ( .NOT. SAMEFILE ) THEN

             M2NAME = PROMPTMFILE(
     &           'Enter name for gridded radiation/cloud input file',
     &           FSREAD3, 'MET_FILE2', 'TMPBIO' )

C......    Read description of radiation/cloud file

             IF ( .NOT. DESC3( M2NAME ) ) THEN

                MESG = 'Could not get description of file "' //
     &             M2NAME( 1:LEN_TRIM( M2NAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

             END IF

C..............  Check that all met files have the same grid.
             CALL CHKGRID( M2NAME, 'GRID' , 0 , EFLAG )

C........  If grid definition does not match BGRD file then stop

             IF ( EFLAG ) THEN
              MESG = 'Problems opening input files. See ERROR(S) above.'
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             END IF

             CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. )
        
           ELSE
             
             M2NAME = M3NAME
 
           END IF
        
        END IF

        IF ( PARTYPE .EQ. 1 ) THEN

C............. Get name of radiation variable to use

            MESG = 'Variable name for radiation'
            CALL ENVSTR( 'RAD_VAR', MESG, 'RGRND', RADNAM, IOS )


            J = INDEX1( RADNAM, NVARS3D, VNAME3D )
            IF ( J .LE. 0 ) THEN

                MESG = 'Could not find ' // RADNAM // 'in file ' //
     &                 M2NAME( 1:LEN_TRIM( M2NAME ) )

                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

            END IF
         END IF

C........  if solar zenith angle calculation needed then get lat-lon
C........  coordinates from GRID_CRO_2D file

        IF ( PARTYPE .GT. 1  ) THEN

            GNAME = PROMPTMFILE( 
     &              'Enter name for 2D GRID PARAMETERS input file',
     &              FSREAD3, 'GRID_CRO_2D', 'TMPBIO' )

            IF ( .NOT. DESC3( GNAME ) ) THEN
                CALL M3EXIT( 'TMPBIO', 0, 0,
     &                      'Could not get description of file "'
     &                      // GNAME( 1:LEN_TRIM( GNAME ) ) // '"', 2 )

            END IF 

C.............  Check that all met files have the same grid.
            CALL CHKGRID( GNAME, 'GRID' , 0 , EFLAG ) 

        END IF

C.......    Get bioseason switch file, BIOSEASON
        IF ( SWITCH_FILE ) THEN

           BNAME = PROMPTMFILE( 
     &          'Enter name for season switch input file',
     &          FSREAD3, 'BIOSEASON', 'TMPBIO' )
           
C......    Read description of switch file

           IF ( .NOT. DESC3( BNAME ) ) THEN

              MESG = 'Could not get description of file "' //
     &             NNAME( 1:LEN_TRIM( BNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

           END IF

C............  Compare grid definition and call with temporary subgrid flag
C              (since this subgrid can be different from final subgrid)
           CALL CHKGRID( BNAME, 'GRID' , 2 , EFLAG )
           SWXOFF = XOFF_A
           SWYOFF = YOFF_A

        END IF

        SWNCOLS = NCOLS3D
        SWNROWS = NROWS3D

C.......   Get normalized emissions file, BGRD

        NNAME = PROMPTMFILE( 
     &          'Enter name for NORMALIZED EMISSIONS input file',
     &          FSREAD3, 'BGRD', 'TMPBIO' )

C......    Read description of normalized emissions file

        IF ( .NOT. DESC3( NNAME ) ) THEN

            MESG = 'Could not get description of file "' //
     &             NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

        END IF

C.........  Final grid definition 
        CALL CHKGRID( NNAME, 'GRID' , 1, EFLAG )

        LUSE  = GETCFDSC( FDESC3D, '/LANDUSE/', .FALSE. )

        IF ( SWITCH_FILE ) THEN

C.......   Get winter normalized emissions file, BGRDW
C.......   Note that second normalized file is assumed to be
C.......   the winter file

           NNAME2 = PROMPTMFILE( 
     &          'Enter name for winter NORMALIZED EMISSIONS input file',
     &          FSREAD3, 'BGRDW', 'TMPBIO' )

C......    Read description of second normalized emissions file

           IF ( .NOT. DESC3( NNAME2 ) ) THEN

              MESG = 'Could not get description of file "' //
     &             NNAME( 1:LEN_TRIM( NNAME2 ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

           END IF

C............  Check normalized emissions are consistent with the subgrid.
           CALL CHKGRID( NNAME2, 'GRID' , 1 , EFLAG )

           LUSE2  = GETCFDSC( FDESC3D, '/LANDUSE/', .FALSE. )

        END IF

C........  If grid definitions do not match properly
        IF ( EFLAG ) THEN
          MESG = 'Problems opening input files. See ERROR(S) above.'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.......   Get default time characteristic for output file:
C.......   If we're going to prompt, then set the defaults based on met
C.......      otherwise, use environment variables to set defaults

        JDATE  = SDATE3D
        JTIME  = STIME3D
        NSTEPS = MXREC3D

        CALL GETM3EPI( TZONE, JDATE, JTIME, NSTEPS )

C.......   Build description for, and create/open output file
C.......   (all but variables-table in description is borrowed from M3NAME)

        SDATE3D = JDATE
        STIME3D = JTIME
        BSTEPS  = NSTEPS
        MXREC3D = BSTEPS
        NVARS3D = MSPCS
        NLAYS3D = 1
        TSTEP3D = 10000

        DO  M = 1, MSPCS
            VNAME3D( M ) = EMSPC( M )
            UNITS3D( M ) = 'moles/hr'
            VDESC3D( M ) = 'biogenic emissions of the indicated species'
            VTYPE3D( M ) = M3REAL
        END DO

        FDESC3D = ' '   ! array

        FDESC3D( 1 ) = 'Gridded biogenic emissions from SMOKE-BEIS2'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        FDESC3D( 4 ) = '/TZONE/ '   // CTZONE
        IF ( SWITCH_FILE ) THEN
           FDESC3D( 5 ) = '/LANDUSE/ ' // LUSE // LUSE2
        ELSE
           FDESC3D( 5 ) = '/LANDUSE/ ' // LUSE
        END IF

        FDESC3D( 6 ) = '/MET SCENARIO/ ' // METSCEN
        FDESC3D( 7 ) = '/CLOUD SCHEME/ ' // CLOUDSHM

        ENAME = PROMPTMFILE(
     &          'Enter name for BGTS output file - moles',
     &          FSUNKN3, 'BGTS_L', 'TMPBIO' )

        DO M = 1, MSPCS
            UNITS3D( M ) = 'tons/hr'
        END DO

        SNAME = PROMPTMFILE(
     &          'Enter name for BGTS output file - mass',
     &          FSUNKN3, 'BGTS_S', 'TMPBIO' )


C........  if solar zenith angle calculation needed then get lat-lon
C........  coordinates from GRID_CRO_2D file
        IF ( PARTYPE .GT. 1  ) THEN

            ALLOCATE( LAT ( METNCOLS, METNROWS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LAT', PROGNAME )

            ALLOCATE( LON ( METNCOLS, METNROWS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LON', PROGNAME )

            IF ( .NOT. READ3( GNAME, 'LAT', 1, 0, 0, LAT ) ) THEN
              MESG = 'Could not read LAT from file "' //
     &                GNAME( 1:LEN_TRIM( GNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( GNAME, 'LON', 1, 0, 0, LON ) ) THEN
              MESG = 'Could not read LON from file "' //
     &                GNAME( 1:LEN_TRIM( GNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

        END IF    ! Use cloud to calculate PAR

C.......   Build name table for variables in normalized emissions file"

        ALLOCATE( NORMV( BTYPES * ( BSPCS - 1 ) ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NORMV', PROGNAME )

        ALLOCATE( NORMN( LUSES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NORMN', PROGNAME )

        I = 0

        DO  B = 1, BSPCS - 1
          DO  K = 1, BTYPES

            I = I + 1
            NORMV( I ) = BIOLTYPE( K ) // BIOSPC( B )

          END DO
        END DO

        DO  L = 1, LUSES

            NORMN( L ) = BIOLUSE( L )( 1:LEN_TRIM( BIOLUSE( L )))//'NO'

        END DO


C.......   Allocate memory for normalized emissions
        ALLOCATE( PINE( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PINE', PROGNAME )

        ALLOCATE( DECD( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DECD', PROGNAME )
        ALLOCATE( CONF( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CONF', PROGNAME )
        ALLOCATE( AGRC( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AGRC', PROGNAME )
        ALLOCATE( LEAF( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LEAF', PROGNAME )
        ALLOCATE( OTHR( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OTHR', PROGNAME )

        ALLOCATE( AVLAI( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVLAI', PROGNAME )

        ALLOCATE( NORNO( NCOLS, NROWS, LUSES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NORNO', PROGNAME )

        IF ( SWITCH_FILE ) THEN
           ALLOCATE( PINEW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'PINEW', PROGNAME )
           ALLOCATE( DECDW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'DECDW', PROGNAME )
           ALLOCATE( CONFW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'CONFW', PROGNAME )
           ALLOCATE( AGRCW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'AGRCW', PROGNAME )
           ALLOCATE( LEAFW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'LEAFW', PROGNAME )
           ALLOCATE( OTHRW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'OTHRW', PROGNAME )
           ALLOCATE( AVLAIW( NCOLS, NROWS ), STAT=IOS )
           CALL CHECKMEM( IOS, 'AVLAIW', PROGNAME )
           ALLOCATE( NORNOW( NCOLS, NROWS, LUSES ), STAT=IOS )
           CALL CHECKMEM( IOS, 'NORNOW', PROGNAME )

           ALLOCATE( PINES( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'PINEW', PROGNAME )
           ALLOCATE( DECDS( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'DECDW', PROGNAME )
           ALLOCATE( CONFS( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'CONFW', PROGNAME )
           ALLOCATE( AGRCS( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'AGRCW', PROGNAME )
           ALLOCATE( LEAFS( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'LEAFW', PROGNAME )
           ALLOCATE( OTHRS( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'OTHRW', PROGNAME )
           ALLOCATE( AVLAIS( NCOLS, NROWS ), STAT=IOS )
           CALL CHECKMEM( IOS, 'AVLAIW', PROGNAME )
           ALLOCATE( NORNOS( NCOLS, NROWS, LUSES ), STAT=IOS )
           CALL CHECKMEM( IOS, 'NORNOW', PROGNAME )

           ALLOCATE( SEASON( SWNCOLS, SWNROWS), STAT=IOS )
           CALL CHECKMEM( IOS, 'SEASON', PROGNAME )
           SEASON = 0   ! array

        END IF

C.......   Loops reading the various categories of normalized emissions:
        I = 0

        DO  M = 1, BSPCS - 1

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0, 
     &                      PINE( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      PINEW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      DECD( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      DECDW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      CONF( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      CONFW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      AGRC( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      AGRCW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF 

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      LEAF( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      LEAFW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      OTHR( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      OTHRW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

        END DO        !  end loop to read normalized emissions 

        IF ( .NOT. READ3( NNAME, 'AVLAI', 1, 0, 0, AVLAI ) ) THEN
             MESG = 'Could not read AVLAI from file "' //
     &              NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
        END IF
        IF ( SWITCH_FILE ) THEN
           IF ( .NOT. READ3( NNAME2, 'AVLAI', 1, 0, 0, 
     &                      AVLAIW ) ) THEN
              MESG = 'Could not read AVLAI from file "' //
     &              NNAME2( 1:LEN_TRIM( NNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
           END IF
        END IF

        I = 0
        DO  J = 1, BSPCS

            I = I + 1
            IF ( .NOT. READ3( NNAME, NORMN( I ), 1, 0, 0, 
     &                      NORNO( 1,1,J ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  NORMN( I )( 1 : LEN_TRIM( NORMN( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMN( I ), 1, 0, 0, 
     &                      NORNOW( 1,1,J ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMN( I )( 1 : LEN_TRIM( NORMN( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

        END DO        !  end loop reading normalized NO's


C.......   Allocate memory for met
        ALLOCATE( TASFC( METNCOLS, METNROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TASFC', PROGNAME )
        ALLOCATE( TSOLAR( METNCOLS, METNROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TSOLAR', PROGNAME )
        ALLOCATE( PRSFC( METNCOLS, METNROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRSFC', PROGNAME )

C.......   Allocate memory for emissions 
        ALLOCATE( EMPOL( NCOLS, NROWS, BSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMPOL', PROGNAME )
        ALLOCATE( EMISL( NCOLS, NROWS, MSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISL', PROGNAME )
        ALLOCATE( EMISS( NCOLS, NROWS, MSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISS', PROGNAME )

C.........  Set up gridded met file(s) dates and times for specific time zone

        IF( M3NAME .NE. ' ' ) THEN
          CALL PREBMET( M3NAME, M2NAME, SAMEFILE, TZONE, TSTEP3D, JDATE, 
     &                  JTIME, NSTEPS, MDATE, MTIME, RDATE, RTIME )
        END IF

C.........  Reset no. output steps if needed
        IF( BSTEPS > NSTEPS ) THEN
           BSTEPS = NSTEPS
        END IF

        IF ( M2NAME .EQ. M3NAME ) THEN

           RDATE = MDATE  
           RTIME = MTIME
  
        END IF 

        LDATE = 0     

C.........  Loop thru the number of time steps (hourly)
 
        DO  HR = 1, BSTEPS

           EMISL = 0        ! array
           EMISS = 0        ! array

           IF( JDATE .NE. LDATE ) THEN

               CALL WRDAYMSG( JDATE, MESG )               

C..........  If new date, read season switch 

               IF ( SWITCH_FILE ) THEN
                 
                  MESG = 'Reading gridded season switch data..'  
                  CALL M3MESG( MESG ) 
                  IF ( .NOT. READ3( BNAME, 'SEASON', 1, 
     &                JDATE, 0, SEASON ) ) THEN
                     MESG = 'Could not read SEASON from file ' //
     &               BNAME( 1:LEN_TRIM( BNAME ) )
                     CALL M3EXIT( 'TMPBIO', JDATE, 0, MESG, 2 )
                  END IF

                  MESG = 'Applying gridded season switch data..' 
                  CALL M3MESG( MESG )


                  DO I = 1, METNCOLS
                    DO J = 1, METNROWS

                        II = I - XOFF
                        JJ = J - YOFF

                        IF( II .LE. 0 .OR. II .GT. NCOLS .OR.
     &                      JJ .LE. 0 .OR. JJ .GT. NROWS       ) CYCLE

C........................  If switch equal to 0 use winter normalized emissions
C........................  Allow for subgrid in switch
C........................  If grid is larger subgrid, then extra cells will
C                          behave as the nearest available cell on the boundary.

                       III = MAX( MIN( I - SWXOFF, SWNCOLS ), 1 )
                       JJJ = MAX( MIN( J - SWYOFF, SWNROWS ), 1 )
                       IF ( SEASON ( III, JJJ ) .EQ. 0 ) THEN
                          DO M = 1, BSPCS - 1
                             PINES( II, JJ, M ) = PINEW( II, JJ, M ) 
                             DECDS( II, JJ, M ) = DECDW( II, JJ, M )
                             CONFS( II, JJ, M ) = CONFW( II, JJ, M )
                             LEAFS( II, JJ, M ) = LEAFW( II, JJ, M )
                             AGRCS( II, JJ, M ) = AGRCW( II, JJ, M )
                             OTHRS( II, JJ, M ) = OTHRW( II, JJ, M )
                          END DO

                          AVLAIS( II, JJ ) = AVLAIW ( II, JJ )

                          DO M = 1, BSPCS
                             NORNOS( II, JJ, M ) = NORNOW( II, JJ, M )
                          END DO

                       ELSE

                          DO M = 1, BSPCS - 1
                             PINES( II, JJ, M ) = PINE( II, JJ, M )
                             DECDS( II, JJ, M ) = DECD( II, JJ, M )
                             CONFS( II, JJ, M ) = CONF( II, JJ, M )
                             LEAFS( II, JJ, M ) = LEAF( II, JJ, M )
                             AGRCS( II, JJ, M ) = AGRC( II, JJ, M )
                             OTHRS( II, JJ, M ) = OTHR( II, JJ, M )
                          END DO

                          AVLAIS( II, JJ ) = AVLAI ( II, JJ )

                          DO M = 1, BSPCS
                             NORNOS( II, JJ, M ) = NORNO( II, JJ, M )
                          END DO

                       END IF

                    END DO
                  END DO
 
               END IF

           END IF

C.............  Write to screen because WRITE3 only writes to LDEV
           WRITE( *, 94030 ) HHMMSS( JTIME )

C.............  Read temperature data

           IF ( .NOT. READ3( M3NAME, TMPRNAM, 1, 
     &          MDATE, MTIME, TASFC ) ) THEN
              MESG = 'Could not read ' // TMPRNAM // 'from file ' //
     &                M3NAME( 1:LEN_TRIM( M3NAME ) )
              CALL M3EXIT( 'TMPBIO', MDATE, MTIME, MESG, 2 )
           END IF


C............  If necessary read solar radiation

           IF ( PARTYPE .EQ. 1 ) THEN      ! Use RADNAM

              IF ( .NOT. READ3( M2NAME, RADNAM, ALLAYS3, RDATE,
     &             RTIME, TSOLAR(1,1) ) ) THEN
                 MESG = 'Could not read ' // RADNAM // 'from file ' //
     &                M2NAME( 1:LEN_TRIM( M2NAME ) ) 

                 CALL M3EXIT( 'TMPBIO', RDATE, RTIME, MESG, 2 )
              END IF

C...............  Convert shortwave radiation to PAR  

              DO I = 1, METNCOLS
               DO J = 1, METNROWS
                    TSOLAR( I, J ) = TSOLAR( I, J ) * SOL2PAR ! Calc. PAR
               END DO
              END DO

C...............  Calculate non-speciated emissions
C...............     must pass met date and time here

              IF ( SWITCH _FILE ) THEN
  
                 CALL HRBIOS( MDATE, MTIME, METNCOLS, METNROWS, 
     &                   PINES, DECDS, CONFS, AGRCS, LEAFS, OTHRS,
     &                   AVLAIS, NORNOS, TASFC, TSOLAR, EMPOL )
          
              ELSE 
 
                 CALL HRBIOS( MDATE, MTIME, METNCOLS, METNROWS,
     &                   PINE, DECD, CONF, AGRC, LEAF, OTHR, AVLAI,
     &                   NORNO, TASFC, TSOLAR, EMPOL )

              END IF

           ELSE        ! Use clouds or clear skies to calculate PAR

C..............  Read surface pressure data from M3NAME

              IF ( .NOT. READ3( M3NAME, 'PRES', 1, 
     &                        MDATE, MTIME, PRSFC ) ) THEN
                  MESG = 'Could not read PRES from file "' //
     &                   M3NAME( 1:LEN_TRIM( M3NAME ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', MDATE, MTIME, MESG, 2 )
              END IF

C...............  convert to millibars

              DO  C = 1, METNCOLS
              DO  R = 1, METNROWS
                       PRSFC( C, R ) = PRSFC( C, R ) * 0.010  ! Pa to mb
              END DO
              END DO

C............. Calculate non-speciated emissions
C............. must pass met date and time here

              IF ( SWITCH_FILE ) THEN

                 CALL HRBIO( MDATE, MTIME, METNCOLS, METNROWS, LAT, LON,  
     &                  PRSFC, TASFC, PINES, DECDS, CONFS, AGRCS, 
     &                  LEAFS, OTHRS, AVLAIS, NORNOS, PARTYPE, EMPOL )
              ELSE
 
                 CALL HRBIO( MDATE, MTIME, METNCOLS, METNROWS, LAT, LON,
     &                  PRSFC, TASFC, PINE, DECD, CONF, AGRC, LEAF,
     &                  OTHR, AVLAI, NORNO, PARTYPE, EMPOL )

              END IF

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
                 END DO
               END DO
             END DO
           END DO

C.............  Write out speciated emissions    

           IF ( .NOT. WRITE3( ENAME, 'ALL', 
     &                        JDATE, JTIME, EMISL ) ) THEN
               CALL M3EXIT( 'TMPBIO', JDATE, JTIME, 
     &                      'Error writing BIO OUTPUT file' , 2 )
           END IF                              !  if write3 failed


           IF ( .NOT. WRITE3( SNAME, 'ALL',
     &                        JDATE, JTIME, EMISS ) ) THEN
               CALL M3EXIT( 'TMPBIO', JDATE, JTIME,
     &                      'Error writing BIO OUTPUT file' , 2 )
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

           END IF
      END DO                !  end loop on hours HR


C.........   End of program:

        CALL M3EXIT( 'TMPBIO', 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )

94030   FORMAT( 8X, 'at time ', A8 )

        END PROGRAM TMPBIO  

