
        SUBROUTINE OPENPDOUT( NPDSRC, NVAR, TZONE, SDATE, STIME, TSTEP,
     &                        FILFMT, TYPNAM, PFLAG, EAIDX,  SPSTAT,
     &                        FNAME, RDEV )

C***********************************************************************
C  subroutine body starts at line 96
C
C  DESCRIPTION:
C      This subroutine opens the day- or hour-specific output files
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines
C
C  REVISION  HISTORY:
C      Created 12/99 by M. Houyoux
C
C       Version 10/2016 by C. Coats:  USE M3UTILIO
C
C****************************************************************************
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
C***************************************************************************
        USE M3UTILIO

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATDESC, CRL, BYEAR, NIPOL, NIACT,
     &                     EANAM

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions
        CHARACTER(16), EXTERNAL :: VERCHAR

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NPDSRC    ! no. period-specific sources
        INTEGER     , INTENT (IN) :: NVAR      ! no. output variables
        INTEGER     , INTENT (IN) :: TZONE     ! time zone of date/time
        INTEGER     , INTENT (IN) :: SDATE     ! Julian start date
        INTEGER     , INTENT (IN) :: STIME     ! start time HHMMSS
        INTEGER     , INTENT (IN) :: TSTEP     ! time step HHMMSS
        INTEGER     , INTENT (IN) :: FILFMT    ! format of period-specific data
        CHARACTER(*), INTENT (IN) :: TYPNAM    ! 'day' or 'hour'
        LOGICAL     , INTENT (IN) :: PFLAG     ! true: creating profiles
        INTEGER     , INTENT (IN) :: EAIDX( NVAR ) ! pol/act index
        INTEGER     , INTENT (IN) :: SPSTAT( MXSPDAT ) ! true: special val exists
        CHARACTER(*), INTENT(OUT) :: FNAME     ! logical file name
        INTEGER     , INTENT(IN OUT) :: RDEV   ! report unit number

C...........   LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.8_Oct2020$' ! CVS release tag

C...........   Other local variables

        INTEGER       J, K, V   ! counter and indices

        CHARACTER(5)       CTZONE      ! string of time zone
        CHARACTER(5)       SCRNAM      ! upcase(typnam)
        CHARACTER(NAMLEN3) VARNAM      ! name for integer index
        CHARACTER(IOVLEN3) VBUF        ! tmp buffer for variable names
        CHARACTER(300)     MESG        ! message buffer

        CHARACTER(16) :: PROGNAME = 'OPENPDOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENPDOUT

C.........  Write time zone to character string
        WRITE( CTZONE,94000 ) TZONE

C.........  Determine integer index variable name
        IF( TYPNAM .EQ. 'day' ) THEN
            VARNAM = 'INDXD'
        ELSE IF( TYPNAM .EQ. 'hour' ) THEN
            VARNAM = 'INDXH'
        ELSE
            MESG = 'INTERNAL ERROR: Type name "'// TYPNAM //
     &             '" not recognized'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set header up with missing
        CALL HDRMISS3

C.........  Set of header with actual settings
        SDATE3D = SDATE
        STIME3D = STIME
        TSTEP3D = TSTEP
        NROWS3D = NPDSRC     !  number of rows = # of period sources

C.........  Define variable for source index
        VNAME3D( 1 ) = VARNAM
        VTYPE3D( 1 ) = M3INT
        UNITS3D( 1 ) = 'n/a'
        VDESC3D( 1 ) = 'source IDs'

C.........  Define hour-specific emissions, if any

        J = 1
        DO V = 1, NVAR
            VBUF = EANAM( EAIDX( V ) )
            J = J + 1
            VNAME3D( J ) = VBUF
            VTYPE3D( J ) = M3REAL

C............. If outputs are profiles instead of data values
            IF( PFLAG ) THEN
                UNITS3D( J ) = 'n/a'
                VDESC3D( J ) = TYPNAM // '-specific ' //
     &                             TRIM( VBUF ) // ' diurnal profile'

C............. If outputs are data values...
            ELSE
                UNITS3D( J ) = 'ton/' // TYPNAM
                VDESC3D( J ) = TYPNAM // '-specific ' //
     &                             TRIM( VBUF ) // ' data'
            END IF

        END DO

C.........  Define hour-specific special data values, if any
        K = 0
        DO V = 1, MXSPDAT
            IF( SPSTAT( V ) > 0 ) THEN
                K = K + 1
                J = J + 1
                VNAME3D( J ) = SPDATNAM( V )
                VTYPE3D( J ) = M3REAL
                UNITS3D( J ) = SPDATUNT( V )

                VDESC3D( J ) = TYPNAM // '-specific ' //
     &                         TRIM( SPDATDSC( V ) ) // ' data'

            END IF
        END DO
        NVARS3D = J

C.........  Define general description info
        FDESC3D( 1 ) = CATDESC // TRIM( TYPNAM ) //
     &                '-specific source inventory'
        FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

        IF( NIPOL .GT. 0 ) THEN
            WRITE( FDESC3D( 4 ),94010 ) '/POLLUTANTS/', NIPOL
        END IF

        IF( NIACT .GT. 0 ) THEN
            WRITE( FDESC3D( 5 ),94010 ) '/ACTIVITIES/', NIACT
        END IF

        IF( K .GT. 0 ) THEN
            WRITE( FDESC3D( 6 ),94010 ) '/SPECIAL DATA/', K
        END IF

        WRITE( FDESC3D( 7 ),94010 ) '/BASE YEAR/ '    , BYEAR
        FDESC3D( 8 ) = '/TZONE/ ' // CTZONE

C.........  Set up default file name and prompting message
        SCRNAM = TYPNAM
        CALL UPCASE( SCRNAM )
        IF( PFLAG ) THEN
            MESG = 'Enter logical name for ' // TRIM( SCRNAM ) //
     &              '-SPECIFIC PROFILES output file'
            FNAME = CRL // TYPNAM // 'PRO'

        ELSE
            MESG = 'Enter logical name for ' // TRIM( SCRNAM ) //
     &              '-SPECIFIC output file'
            FNAME = CRL // TYPNAM
        END IF

C.........  Prompt for output file
        CALL UPCASE( FNAME )
        FNAME = PROMPTMFILE( MESG, FSUNKN3, FNAME, PROGNAME )

C.........  If format is CEM format, prompt for report output file name
        IF ( FILFMT .EQ. CEMFMT .AND. RDEV .LE. 0 ) THEN
            MESG = 'Enter logical name for the CEM MATCHING REPORT'
            RDEV = PROMPTFFILE( MESG, .FALSE., .TRUE.,
     &                          'REPINVEN', PROGNAME )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE OPENPDOUT

