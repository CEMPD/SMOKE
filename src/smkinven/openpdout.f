
        SUBROUTINE OPENPDOUT( NPDSRC, NVAR, TZONE, SDATE, STIME, TSTEP, 
     &                        TYPNAM, PFLAG, EAIDX, FNAME )

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
C****************************************************************************/
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL        PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NPDSRC    ! no. period-specific sources
        INTEGER     , INTENT (IN) :: NVAR      ! no. output variables
        INTEGER     , INTENT (IN) :: TZONE     ! time zone of date/time
        INTEGER     , INTENT (IN) :: SDATE     ! Julian start date
        INTEGER     , INTENT (IN) :: STIME     ! start time HHMMSS
        INTEGER     , INTENT (IN) :: TSTEP     ! time step HHMMSS
        CHARACTER(*), INTENT (IN) :: TYPNAM    ! 'day' or 'hour'
        LOGICAL     , INTENT (IN) :: PFLAG     ! true: creating profiles
        INTEGER     , INTENT (IN) :: EAIDX( NVAR ) ! pol/act index
        CHARACTER(*), INTENT(OUT) :: FNAME     ! logical file name

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: SCCSW  = '@(#)$Id$'  ! SCCS string with vers no.

C...........   Other local variables

        INTEGER       L, L2, V      ! counter and indices

        CHARACTER*5            CTZONE      ! string of time zone
        CHARACTER(LEN=NAMLEN3) VARNAM      ! name for integer index
        CHARACTER(LEN=IOVLEN3) VBUF        ! tmp buffer for variable names
        CHARACTER*300          MESG        ! message buffer 

        CHARACTER*16 :: PROGNAME = 'OPENPDOUT' ! program name

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
        NVARS3D = NVAR + 1
        NROWS3D = NPDSRC     !  number of rows = # of period sources

        VNAME3D( 1 ) = VARNAM
        VTYPE3D( 1 ) = M3INT
        UNITS3D( 1 ) = 'n/a'
        VDESC3D( 1 ) = 'source IDs'

        DO V = 1, NVAR 
            VBUF = EANAM( EAIDX( V ) )
            L = LEN_TRIM( VBUF )

            VNAME3D( V + 1 ) = VBUF
            VTYPE3D( V + 1 ) = M3REAL

C............. If outputs are profiles instead of data values
            IF( PFLAG ) THEN
        	UNITS3D( V + 1 ) = 'n/a'
        	VDESC3D( V + 1 ) = TYPNAM // '-specific ' // 
     &                             VBUF( 1:L ) // ' diurnal profile'

C............. If outputs are data values...
            ELSE
        	UNITS3D( V + 1 ) = 'tons/hour'
        	VDESC3D( V + 1 ) = TYPNAM // '-specific ' //
     &                             VBUF( 1:L ) // ' data'
            END IF

   	END DO     

        L2 = LEN_TRIM( TYPNAM )
        FDESC3D( 1 ) = CATDESC // TYPNAM( 1:L2 ) //
     &                '-specific source inventory'
        FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )

        IF( NIPOL .GT. 0 ) THEN
            WRITE( FDESC3D( 4 ),94010 ) '/POLLUTANTS/', NIPOL  
        END IF

        IF( NIACT .GT. 0 ) THEN
            WRITE( FDESC3D( 5 ),94010 ) '/ACTIVITIES/', NIACT
        END IF

        WRITE( FDESC3D( 6 ),94010 ) '/BASE YEAR/ '    , BYEAR 
        FDESC3D( 7 ) = '/TZONE/ ' // CTZONE

C.........  Set up default file name and prompting message
        CALL UPCASE( TYPNAM )
        IF( PFLAG ) THEN
            MESG = 'Enter logical name for ' // TYPNAM( 1:L2 ) //
     &              '-SPECIFIC PROFILES output file'
            FNAME = CRL // TYPNAM // 'PRO'

        ELSE
            MESG = 'Enter logical name for ' // TYPNAM( 1:L2 ) //
     &              '-SPECIFIC output file'
            FNAME = CRL // TYPNAM

        END IF

C.........  Prompt for output file
        FNAME = PROMPTMFILE( MESG, FSUNKN3, FNAME, PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENPDOUT

