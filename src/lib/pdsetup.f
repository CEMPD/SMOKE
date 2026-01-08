
        SUBROUTINE PDSETUP( INFILE, ESDATE, ESTIME, EEDATE, EETIME, 
     &                      TZONE, NIPPA, EANAM, NPOA, NENTRY, 
     &                      EFLAG, PNAME, PDESC )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Checks and sets up to use part-day emission files (day-specific and 
C      hour-specific files) for point sources.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C**************************************************************************
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
C***************************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDE FILES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
C       CHARACTER(2)    CRLF
        INTEGER         GETIFDSC
C       CHARACTER(10)   HHMMSS
C       INTEGER         INDEX1
C       CHARACTER(14)   MMDDYY
C       INTEGER         SECSDIFF

C        EXTERNAL        CRLF, GETIFDSC, HHMMSS, INDEX1, MMDDYY, SECSDIFF
        EXTERNAL     GETIFDSC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT(IN   ) :: INFILE ! name of day- or hour-spec file
        INTEGER     , INTENT(IN   ) :: ESDATE ! episode starting date
        INTEGER     , INTENT(IN   ) :: ESTIME ! episode starting time
        INTEGER     , INTENT(IN   ) :: EEDATE ! episode ending date
        INTEGER     , INTENT(IN   ) :: EETIME ! episode ending time
        INTEGER     , INTENT(IN   ) :: TZONE  ! output time zone
        INTEGER     , INTENT(IN   ) :: NIPPA  ! number of inventory pols/acts
        CHARACTER(*), INTENT(IN   ) :: EANAM( NIPPA ) ! names of pols/acts
        INTEGER     , INTENT(  OUT) :: NPOA   ! no. pols/acts
        INTEGER     , INTENT(  OUT) :: NENTRY ! no. entries in file
        LOGICAL     , INTENT(INOUT) :: EFLAG  ! True if error occurred in routine
        CHARACTER(*), INTENT(  OUT) :: PNAME( * ) ! names of time-spec pol/acts
        CHARACTER(*), INTENT(  OUT) :: PDESC( * ) ! descriptions of pol/acts

C...........   Other local variables
        INTEGER         I, J, LD, L2, N      !  counters and indices

        INTEGER         EDATE            !  ending date
        INTEGER         ETIME            !  ending time
        INTEGER         IOS              !  i/o status
        INTEGER         NSTEPS           !  number of time steps
        INTEGER         NVARS            !  tmp no. of variables in file
        INTEGER         SDATE            !  starting date
        INTEGER         STIME            !  starting time
        INTEGER         ZONE             !  time zone

        CHARACTER(4)    DESCRIBE 
        CHARACTER(10)   CTIMESTR, ETIMESTR   ! inventory and comparison time
        CHARACTER(14)   CDATESTR, EDATESTR   ! inventory and comparison date
        CHARACTER(300)  MESG 

        CHARACTER(IOVLEN3)   VBUF        ! tmp variable name

        CHARACTER(16) :: PROGNAME = 'PDSETUP' ! program name

C***********************************************************************
C   begin body of subroutine PDSETUP

C.........  NOTE- FDESC3 variables passed through I/O API common in include file
C.........  Make settings that depend on day-specific or hour-specific data
        IF( TSTEP3D .EQ. 240000 ) THEN  ! day specific
            NSTEPS = 24 * MXREC3D  
            DESCRIBE = 'Day'
            LD = LEN_TRIM( DESCRIBE )

        ELSE ! hour specific
            NSTEPS= MXREC3D
            DESCRIBE = 'Hour'
            LD = LEN_TRIM( DESCRIBE )

        END IF

C.........  Ensure that the time zone of the file is consistent with the
C           time zone of the output data
        ZONE = GETIFDSC( FDESC3D, '/TZONE/', .TRUE. )

        IF( ZONE .NE. TZONE ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Time zone in ' //
     &             DESCRIBE( 1:LD ) // '-specific file is', ZONE,
     &             CRLF() // BLANK10 // 'and is inconsistent ' //
     &             'with output time zone ', TZONE
            CALL M3MSG2( MESG )

        END IF

C.........  Store header info for memory allocation and flexible usage
C               of pollutant-specific data
        SDATE = SDATE3D
        STIME = STIME3D
        NVARS = NVARS3D
        NENTRY= NROWS3D

C.........  Set ending date/time of file
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, NSTEPS * 10000 )

C.........  Compare day-specific dates/times to epsiode dates/times
        I = SECSDIFF( ESDATE, ESTIME, SDATE, STIME )
        J = SECSDIFF( EEDATE, EETIME, EDATE, ETIME )

        IF( I .GT. 0 ) THEN
            EDATESTR = MMDDYY( ESDATE )
            ETIMESTR = HHMMSS( ESTIME )
            CDATESTR = MMDDYY( SDATE )
            CTIMESTR = HHMMSS( STIME )
            MESG = 'WARNING: ' // DESCRIBE( 1:LD ) // '-specific ' //
     &             'starting date/time of '// CDATESTR// '@ ' //
     &             CTIMESTR // CRLF() // BLANK10 // 
     &             'is later than episode starting date/time of ' //
     &             EDATESTR // '@ ' // ETIMESTR 

            CALL M3MSG2( MESG )
        END IF

        IF( J .LT. 0 ) THEN
            EDATESTR = MMDDYY( EEDATE )
            ETIMESTR = HHMMSS( EETIME )
            CDATESTR = MMDDYY( EDATE )
            CTIMESTR = HHMMSS( ETIME )
            MESG = 'WARNING: ' // DESCRIBE( 1:LD ) // '-specific ' //
     &             'ending date/time of '// CDATESTR// '@ ' //
     &             CTIMESTR // CRLF() // BLANK10 // 
     &             'is earlier than episode ending date/time of ' //
     &             EDATESTR // '@ ' // ETIMESTR 

            CALL M3MSG2( MESG )
        END IF

C.............  Store list of valid pollutant names, or...
C.............  Report which period-specific variables will be ignored
C.............  Skip first variable, which is the source index
        N = 0
        DO I = 2, NVARS

            VBUF = VNAME3D( I )
            J = INDEX1( VBUF, NIPPA, EANAM )
            IF( J .LE. 0 ) THEN
                L2 = LEN_TRIM( VBUF )
                MESG = 'WARNING: ' // DESCRIBE( 1:LD ) // 
     &                 '-specific pollutant "' // VBUF( 1:L2 ) // 
     &                 '" is not in inventory file, so ' //
     &                 'it will be ignored.'
                CALL M3MSG2( MESG )

            ELSE
                N = N + 1
                PNAME( N ) = VBUF
                PDESC( N ) = VDESC3D( I )

            END IF
        END DO

        NPOA = N

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE PDSETUP
