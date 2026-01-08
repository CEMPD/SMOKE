
        SUBROUTINE OPENLAYOUT( SDATE, STIME, TSTEP, EMLAYS, REP_LAYR, 
     &                         EXPLONLY, INVPROG, INVVERS, METSCEN, 
     &                         CLOUDSHM, VGLVSXG, GFLAG, GRDNM, LNAME, 
     &                         RDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C    Opens the output layer fractions file.  Either the full file or the
C    explicit-only file are opened, depending on input settings. Also, the
C    report file is opened.
C
C  PRECONDITIONS REQUIRED:
C    Input files opened. Input control settings evaluated.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       I/O API 
CC
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***********************************************************************
 
C...........  MODULES for public variables
C.........  This module contains arrays for plume-in-grid and major sources
        USE M3UTILIO

        USE MODELEV, ONLY: NHRSRC

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: VGTYP, VGTOP

C...........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC

        IMPLICIT NONE
 
C...........   INCLUDES:
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'IOSTRG3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'CONST3.EXT'    ! physical and mathematical constants

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(16)   VERCHAR
C       INTEGER         PROMPTFFILE
C       CHARACTER(16)   PROMPTMFILE

C        EXTERNAL        VERCHAR, PROMPTFFILE, PROMPTMFILE
        EXTERNAL     VERCHAR

C...........  LOCAL PARAMETERS and their descriptions:
C       CHARACTER(50), PARAMETER :: 
C    &  CVSW = '$Name SMOKEv5.2.1_Sep2025$' ! CVS release tag

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: SDATE    ! Julian start date
        INTEGER     , INTENT (IN) :: STIME    ! start time (HHMMSS)
        INTEGER     , INTENT (IN) :: TSTEP    ! time step (HHMMSS)
        INTEGER     , INTENT (IN) :: EMLAYS   ! number of emissions layers
        INTEGER     , INTENT (IN) :: REP_LAYR ! layer for reporting
        LOGICAL     , INTENT (IN) :: EXPLONLY ! true: open explicit-only file
        CHARACTER(*), INTENT (IN) :: INVPROG  ! inventory program
        CHARACTER(*), INTENT (IN) :: INVVERS  ! inventory program version
        CHARACTER(*), INTENT (IN) :: METSCEN  ! met scenario name
        CHARACTER(*), INTENT (IN) :: CLOUDSHM ! cloud scheme name
        REAL        , INTENT (IN) :: VGLVSXG( 0:MXLAYS3 ) !  vertical coord values
        LOGICAL     , INTENT (IN) :: GFLAG    ! true: using variable grid
        CHARACTER(*), INTENT (IN) :: GRDNM    ! grid name
        CHARACTER(*), INTENT(OUT) :: LNAME    ! layer fractions logical file nam
        INTEGER     , INTENT(OUT) :: RDEV     ! report unit number

C...........   Local variables

        INTEGER       J

        CHARACTER(300)   MESG      !  buffer for M3EXIT() messages

        CHARACTER(16) :: PROGNAME = 'OPENLAYOUT'   !  program name

C***********************************************************************
C   begin body of subroutine OPENLAYOUT

C.........  Set up and open output file, which will primarily using I/O API 
C           settings from the cross-point met file (XNAME), which are 
C           already retrieved
        CALL HDRMISS3 

        SDATE3D = SDATE
        STIME3D = STIME
        TSTEP3D = TSTEP
        NLAYS3D = EMLAYS

        J = LBOUND( VGLVS3D, 1 )
        VGLVS3D( J:J+EMLAYS ) = VGLVSXG( 0:EMLAYS )  ! array
        VGTYP3D = VGTYP
        VGTOP3D = VGTOP

        FDESC3D = ' '  ! array

        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        FDESC3D( 4 ) = '/MET SCENARIO/ ' // METSCEN
        FDESC3D( 5 ) = '/CLOUD SCHEME/ ' // CLOUDSHM

        FDESC3D( 11 ) = '/PNTS FROM/ ' // INVPROG
        FDESC3D( 12 ) = '/PNTS VERSION/ ' // INVVERS
        
        IF( GFLAG ) THEN
            FDESC3D( 13 ) = '/VARIABLE GRID/ ' // GRDNM
        END IF

C.........  Settings that depend on whether the output file is for all sources
C           or only explicit sources
C.........  For explicit only ...
        IF ( EXPLONLY ) THEN 

            NROWS3D = NHRSRC
            NVARS3D = 2
            NTHIK3D = NSRC

            VNAME3D( 1 ) = 'INDXH'
            VTYPE3D( 1 ) = M3INT
            UNITS3D( 1 ) = 'none'
            VDESC3D( 1 ) = 'Source number'

            VNAME3D( 2 ) = 'LFRAC'
            VTYPE3D( 2 ) = M3REAL
            UNITS3D( 2 ) = 'none'
            VDESC3D( 2 ) = 'Fraction of plume emitted into layer'

            FDESC3D( 1 ) = 'Explicit sources hourly plume rise layer '//
     &                     'fractions'

            MESG = 'Enter logical name for EXPLICIT LAYER FRACTIONS ' //
     &             'MATRIX'
            LNAME = PROMPTMFILE( MESG, FSUNKN3, 'PLAY_EX', PROGNAME )

C.........  For all sources...
        ELSE
            NROWS3D = NSRC
            NVARS3D = 1

            VNAME3D( 1 ) = 'LFRAC'
            VTYPE3D( 1 ) = M3REAL
            UNITS3D( 1 ) = 'none'
            VDESC3D( 1 ) = 'Fraction of plume emitted into layer'

            FDESC3D( 1 ) = 'By-source hourly plume rise layer ' //
     &                     'fractions'

            MESG = 'Enter logical name for LAYER FRACTIONS MATRIX'
            LNAME = PROMPTMFILE( MESG, FSUNKN3, 'PLAY', PROGNAME )

        END IF

C.........  Get file name of report of plume exceeding specified layer
C.........  Write header to the report
        IF( REP_LAYR .GT. 0 ) THEN

            WRITE( MESG,94010 ) 'Enter logical name for report of ' //
     &                          'plumes exceeding layer', REP_LAYR
            RDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 
     &                          'REPRTLAY', PROGNAME )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

        END SUBROUTINE OPENLAYOUT
