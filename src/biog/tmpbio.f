
        PROGRAM TMPBIO

C***********************************************************************
C  program body starts at line  187
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
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
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
        USE MODSPRO

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     !
        INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants


C...........   PARAMETERS and their descriptions:

        CHARACTER*50  SCCSW          ! SCCS string with version number at end
 
        PARAMETER ( SCCSW   = '@(#)$Id$'    )
     
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
        INTEGER         TRIMLEN
        CHARACTER*16    VERCHAR

        EXTERNAL        CRLF, ENVINT, ENVYN, GETDATE, GETFLINE, GETNUM, 
     &                  GETYN, HHMMSS, INDEX1, PROMPTMFILE, 
     &                  PROMPTFFILE, TRIMLEN, VERCHAR

C.........  Gridded meteorology data
                
        REAL, ALLOCATABLE :: LAT  ( :, : )    !  grid lat (deg) -90 <= LAT <= 90
        REAL, ALLOCATABLE :: LON  ( :, : )    !  grid lon (deg) -180 <= LON <= 180 
        REAL, ALLOCATABLE :: TASFC ( :, : )    !  level-1 air  temperature (K)
        REAL, ALLOCATABLE :: PRSFC  ( :, : )    !  pressure (Pa)
        REAL, ALLOCATABLE :: TSOLAR ( :, :)     !  Photosynthetic Active Radiation (PAR)

C.......   Gridded normalized emissions description input from file BGRD

        REAL, ALLOCATABLE ::  PINE( :, :, : )         !   for pine
        REAL, ALLOCATABLE ::  DECD( :, :, : )         !   for deciduous forest
        REAL, ALLOCATABLE ::  CONF( :, :, : )         !   for coniferous forest
        REAL, ALLOCATABLE ::  AGRC( :, :, : )         !   for agricultural land
        REAL, ALLOCATABLE ::  LEAF( :, :, : )         !   for leaf area
        REAL, ALLOCATABLE ::  OTHR( :, :, : )         !   for other land

        REAL, ALLOCATABLE ::  AVLAI( :, : )           !  average LAI
        REAL, ALLOCATABLE ::  NORNO( :, :, : )        !  normalized NO emissions
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
        CHARACTER*16    GNAME   !  logical name for GRID_CRO_2D
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
        INTEGER         NSTEPS  !  run duration (hours)
        INTEGER         PARTYPE !  method number to calculate PAR
        INTEGER         RDATE   !  met file 2 start date 
        INTEGER         RTIME   !  met file 2 start time
        INTEGER         TZONE   !  output-file time zone ; not used in program

        LOGICAL         EFLAG   !  error flag
        LOGICAL ::      SAMEFILE = .TRUE.   ! radiation/cld and tmpr data in same file 

        CHARACTER*300   MESG    !  message buffer for M3EXIT()

C.......   Input met and grid variables:
        CHARACTER*16 :: PROGNAME = 'TMPBIO'   !  program name

C***********************************************************************
C   begin body of program TMPBIO

        LDEV = INIT3()
 
C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Evaluate the environment variables...

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

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
            ENDIF
          ENDDO
        ENDDO

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

C.......   Get normalized emissions file, BGRD

        NNAME = PROMPTMFILE( 
     &          'Enter name for NORMALIZED EMISSIONS input file',
     &          FSREAD3, 'BGRD', 'TMPBIO' )

C......    Read description of temperature file

        IF ( .NOT. DESC3( NNAME ) ) THEN

            MESG = 'Could not get description of file "' //
     &             NNAME( 1:TRIMLEN( NNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

        ENDIF

        CALL CHKGRID( NNAME, 'GRID' , EFLAG )
 
        NCOLS = NCOLS3D 
        NROWS = NROWS3D
        NGRID = NCOLS3D * NROWS3D
        LUSE  = GETCFDSC( FDESC3D, '/LANDUSE/', .FALSE. )
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
     &             M3NAME( 1:TRIMLEN( M3NAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

        ENDIF

        J = INDEX1( TMPRNAM , NVARS3D, VNAME3D )

        IF ( J .LE. 0 ) THEN

             MESG = 'Could not find ' // TMPRNAM // 'in file ' //
     &                   M3NAME( 1:TRIMLEN( M2NAME ) ) 

             CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
 
        ENDIF

        CALL CHKGRID( M3NAME, 'GRID' , EFLAG ) 

C........  Get met and cloud scheme descriptions from M3NAME if they exist

        METSCEN  = GETCFDSC( FDESC3D, '/MET SCENARIO/', .FALSE. )
        CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. )

        IF ( PARTYPE .NE. 5 ) THEN

C.........  Open second met file if needed

           IF ( .NOT. SAMEFILE ) THEN

             M2NAME = PROMPTMFILE(
     &           'Enter name for gridded radiation/cloud input file',
     &           FSREAD3, 'MET_FILE2', 'TMPBIO' )

C......    Read description of temperature file

             IF ( .NOT. DESC3( M2NAME ) ) THEN

                MESG = 'Could not get description of file "' //
     &             M2NAME( 1:TRIMLEN( M2NAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

             ENDIF

             CALL CHKGRID( M2NAME, 'GRID' , EFLAG )
             CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. )
        
           ELSE
             
             M2NAME = M3NAME
 
           ENDIF
        
        ENDIF

        IF ( PARTYPE .EQ. 1 ) THEN

C......    Get name of radiation variable to use

            MESG = 'Variable name for radiation'
            CALL ENVSTR( 'RAD_VAR', MESG, 'RGRND', RADNAM, IOS )


            J = INDEX1( RADNAM, NVARS3D, VNAME3D )
            IF ( J .LE. 0 ) THEN

                MESG = 'Could not find ' // RADNAM // 'in file ' //
     &                 M2NAME( 1:TRIMLEN( M2NAME ) )

                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

            ENDIF
         ENDIF

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
        MXREC3D = NSTEPS
        NVARS3D = MSPCS
        NLAYS3D = 1
        TSTEP3D = 10000

        DO  M = 1, MSPCS
            VNAME3D( M ) = EMSPC( M )
            UNITS3D( M ) = 'moles/hr'
            VDESC3D( M ) = 'biogenic emissions of the indicated species'
            VTYPE3D( M ) = M3REAL
        ENDDO

        FDESC3D = ' '   ! array

        FDESC3D( 1 ) = 'Gridded biogenic emissions from SMOKE-BEIS2'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )
        FDESC3D( 4 ) = '/TZONE/ '   // CTZONE
        FDESC3D( 5 ) = '/LANDUSE/ ' // LUSE
        FDESC3D( 6 ) = '/MET SCENARIO/ ' // METSCEN
        FDESC3D( 7 ) = '/CLOUD SCHEME/ ' // CLOUDSHM

        ENAME = PROMPTMFILE(
     &          'Enter name for BGTS output file - moles',
     &          FSUNKN3, 'BGTS_L', 'TMPBIO' )

        DO M = 1, MSPCS
            UNITS3D( M ) = 'ton/hr'
        ENDDO

        SNAME = PROMPTMFILE(
     &          'Enter name for BGTS output file - mass',
     &          FSUNKN3, 'BGTS_S', 'TMPBIO' )


C........  if solar zenith angle calculation needed then get lat-lon
C........  coordinates from GRID_CRO_2D file

        IF ( PARTYPE .GT. 1  ) THEN

            GNAME = PROMPTMFILE( 
     &              'Enter name for 2D GRID PARAMETERS input file',
     &              FSREAD3, 'GRID_CRO_2D', 'TMPBIO' )

            IF ( .NOT. DESC3( GNAME ) ) THEN
                CALL M3EXIT( 'TMPBIO', 0, 0,
     &                      'Could not get description of file "'
     &                      // GNAME( 1:TRIMLEN( GNAME ) ) // '"', 2 )

            ENDIF 
            CALL CHKGRID( GNAME, 'GRID' , EFLAG ) 

            ALLOCATE( LAT ( NCOLS, NROWS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LAT', PROGNAME )

            ALLOCATE( LON ( NCOLS, NROWS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LON', PROGNAME )

            IF ( .NOT. READ3( GNAME, 'LAT', 1, 0, 0, LAT ) ) THEN
              MESG = 'Could not read LAT from file "' //
     &                GNAME( 1:TRIMLEN( GNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( GNAME, 'LON', 1, 0, 0, LON ) ) THEN
              MESG = 'Could not read LON from file "' //
     &                GNAME( 1:TRIMLEN( GNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            ENDIF

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

          ENDDO
        ENDDO

        DO  L = 1, LUSES

            NORMN( L ) = BIOLUSE( L )( 1:TRIMLEN( BIOLUSE( L )))//'NO'

        ENDDO


C.......   Loops reading the various categories of normalized emissions:

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

        I = 0

        DO  M = 1, BSPCS - 1

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0, 
     &                      PINE( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : TRIMLEN( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      DECD( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : TRIMLEN( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      CONF( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : TRIMLEN( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      AGRC( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : TRIMLEN( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
 
            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      LEAF( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : TRIMLEN( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      OTHR( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : TRIMLEN( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
 

        ENDDO        !  end loop to read normalized emissions 

        IF ( .NOT. READ3( NNAME, 'AVLAI', 1, 0, 0, AVLAI ) ) THEN
             MESG = 'Could not read AVLAI from file "' //
     &              NNAME( 1:TRIMLEN( NNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
        END IF

        I = 0
        DO  J = 1, BSPCS

            I = I + 1
            IF ( .NOT. READ3( NNAME, NORMN( I ), 1, 0, 0, 
     &                      NORNO( 1,1,J ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  NORMN( I )( 1 : TRIMLEN( NORMN( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

        ENDDO        !  end loop reading normalized NO's


C.......   Allocate memory for met and emissions 

        ALLOCATE( TASFC( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TASFC', PROGNAME )

        ALLOCATE( TSOLAR( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TSOLAR', PROGNAME )

        ALLOCATE( EMPOL( NCOLS, NROWS, BSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMPOL', PROGNAME )

        ALLOCATE( PRSFC( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRSFC', PROGNAME )

        ALLOCATE( EMISL( NCOLS, NROWS, MSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISL', PROGNAME )
        ALLOCATE( EMISS( NCOLS, NROWS, MSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISS', PROGNAME )

C.........  Set up gridded met file(s) dates and times for specific time zone

        IF( M3NAME .NE. ' ' ) THEN
          CALL PREBMET( M3NAME, M2NAME, SAMEFILE, TZONE, TSTEP3D, JDATE, 
     &                  JTIME, NSTEPS, MDATE, MTIME, RDATE, RTIME )
        END IF

        IF ( M2NAME .EQ. M3NAME ) THEN

           RDATE = MDATE  
           RTIME = MTIME
  
        ENDIF 

        LDATE = 0     

C.........  Loop thru the number of time steps
 
        DO  HR = 1, NSTEPS

           EMISL = 0
           EMISS = 0

           IF( JDATE .NE. LDATE ) THEN

               CALL WRDAYMSG( JDATE, MESG )               

           ENDIF

C.............  Write to screen because WRITE3 only writes to LDEV
           WRITE( *, 94030 ) HHMMSS( JTIME )

C.............  Read temperature data

           IF ( .NOT. READ3( M3NAME, TMPRNAM, 1, 
     &          MDATE, MTIME, TASFC ) ) THEN
              MESG = 'Could not read ' // TMPRNAM // 'from file ' //
     &                M3NAME( 1:TRIMLEN( M3NAME ) )
              CALL M3EXIT( 'TMPBIO', MDATE, MTIME, MESG, 2 )
           END IF


C............. If necessary read solar radiation

           IF ( PARTYPE .EQ. 1 ) THEN      ! Use RADNAM

              IF ( .NOT. READ3( M2NAME, RADNAM, ALLAYS3, RDATE,
     &             RTIME, TSOLAR(1,1) ) ) THEN
                 MESG = 'Could not read ' // RADNAM // 'from file ' //
     &                M2NAME( 1:TRIMLEN( M2NAME ) ) 

                 CALL M3EXIT( 'TMPBIO', RDATE, RTIME, MESG, 2 )
              END IF

C............. Convert shortwave radiation to PAR  

              DO I = 1, NCOLS
               DO J = 1, NROWS
                    TSOLAR( I, J ) = TSOLAR( I, J ) * SOL2PAR ! Calc. PAR
               ENDDO
              END DO

C............. Calculate non-speciated emissions
C............. must pass met date and time here
 
              CALL HRBIOS( MDATE, MTIME, NCOLS, NROWS,
     &                   PINE, DECD, CONF, AGRC, LEAF, OTHR, AVLAI,
     &                   NORNO, TASFC, TSOLAR, EMPOL )

           ELSE        ! Use clouds or clear skies to calculate PAR

C..............  Read surface pressure data from M3NAME

              IF ( .NOT. READ3( M3NAME, 'PRES', 1, 
     &                        MDATE, MTIME, PRSFC ) ) THEN
                  MESG = 'Could not read PRES from file "' //
     &                   M3NAME( 1:TRIMLEN( M3NAME ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', MDATE, MTIME, MESG, 2 )
              END IF

C...............  convert to millibars

              DO  C = 1, NCOLS
              DO  R = 1, NROWS
                       PRSFC( C, R ) = PRSFC( C, R ) * 0.010  ! Pa to mb
              ENDDO
              ENDDO

C............. Calculate non-speciated emissions
C............. must pass met date and time here

              CALL HRBIO( MDATE, MTIME, NCOLS, NROWS, LAT, LON, PRSFC, 
     &                    TASFC, PINE, DECD, CONF, AGRC, LEAF, OTHR,  
     &                    AVLAI, NORNO, PARTYPE, EMPOL )
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

           ENDIF
      ENDDO                !  end loop on hours HR


C.........   End of program:

        CALL M3EXIT( 'TMPBIO', 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )

94030   FORMAT( 8X, 'at time ', A8 )

        END PROGRAM TMPBIO  

