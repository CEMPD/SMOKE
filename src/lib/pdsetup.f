
        SUBROUTINE PDSETUP( INFILE, ESDATE, ESTIME, EEDATE, EETIME, 
     &                      NIPOL, EINAM, NVARS, NENTRY, EFLAG, PNAME )

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
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

        IMPLICIT NONE

C...........   INCLUDE FILES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2     CRLF
        CHARACTER*10    HHMMSS
        CHARACTER*14    MMDDYY
        INTEGER         SECSDIFF

        EXTERNAL        CRLF, HHMMSS, MMDDYY, SECSDIFF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: INFILE ! name of day- or hour-spec file	
        INTEGER     , INTENT (IN) :: ESDATE ! episode starting date
        INTEGER     , INTENT (IN) :: ESTIME ! episode starting time
        INTEGER     , INTENT (IN) :: EEDATE ! episode ending date )
        INTEGER     , INTENT (IN) :: EETIME ! episode ending time (in)
        INTEGER     , INTENT (IN) :: NIPOL  ! number of inventory pollutants
        CHARACTER(*), INTENT (IN) :: EINAM( NIPOL ) ! names of inventory pols
	INTEGER     , INTENT(OUT) :: NVARS  ! no. variables in file
	INTEGER     , INTENT(OUT) :: NENTRY ! no. entries in file
	LOGICAL     , INTENT(OUT) :: EFLAG  ! True if error occurred in routine
        CHARACTER(*), INTENT(OUT) :: PNAME( * ) ! names of time-specific pols
 
C...........   Index of variable names in EINAM
        INTEGER       INDX( NIPOL )

C...........   Other local variables
        INTEGER         I, J, L, L2      !  counters and indices

        INTEGER         EDATE            !  ending date
        INTEGER         ETIME            !  ending time
        INTEGER         IOS              !  i/o status
        INTEGER         NSTEPS           !  number of time steps
        INTEGER         SDATE            !  starting date
        INTEGER         STIME            !  starting time

        CHARACTER*4     DESCRIBE 
        CHARACTER*10    CTIMESTR, ETIMESTR   ! inventory and comparison time
        CHARACTER*14    CDATESTR, EDATESTR   ! inventory and comparison date
        CHARACTER*300   MESG 

        CHARACTER*16 :: PROGNAME = 'PDSETUP' ! program name

C***********************************************************************
C   begin body of subroutine PDSETUP

C.........  Initialize index
        INDX = 0    ! array

C.........  NOTE: FDESC3 variables passed through I/O API common in include file
C.........  Make settings that depend on day-specific or hour-specific data
        IF( TSTEP3D .EQ. 240000 ) THEN  ! day specific
            NSTEPS = 24 * MXREC3D  
            DESCRIBE = 'Day'

        ELSE ! hour specific
            NSTEPS= MXREC3D
            DESCRIBE = 'Hour'

        ENDIF

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
        I = SECSDIFF( ESDATE, ESTIME, SDATE, SDATE )
        J = SECSDIFF( EEDATE, EETIME, EDATE, EDATE )

        L = LEN_TRIM( DESCRIBE )
        IF( I .LT. 0 ) THEN
            EFLAG = .TRUE.
            EDATESTR = MMDDYY( ESDATE )
            ETIMESTR = HHMMSS( ESTIME )
            CDATESTR = MMDDYY( SDATE )
            CTIMESTR = HHMMSS( STIME )
            MESG = 'ERROR: ' // DESCRIBE( 1:L ) // '-specific ' //
     &             'starting date/time of '// CDATESTR// '/' //
     &             CTIMESTR // CRLF() // BLANK10 // 
     &             'is later than episode starting date/time of ' //
     &             EDATESTR // '/' // ETIMESTR // '.'

            CALL M3MSG2( MESG )
        ENDIF

        IF( J .GT. 0 ) THEN
            EFLAG = .TRUE.
            EDATESTR = MMDDYY( EEDATE )
            ETIMESTR = HHMMSS( EETIME )
            CDATESTR = MMDDYY( EDATE )
            CTIMESTR = HHMMSS( ETIME )
            MESG = 'ERROR: ' // DESCRIBE( 1:L ) // '-specific ' //
     &             'ending date/time of '// CDATESTR// '/' //
     &             CTIMESTR // CRLF() // BLANK10 // 
     &             'is earlier than episode ending date/time of ' //
     &             EDATESTR // '/' // ETIMESTR // '.'

            CALL M3MSG2( MESG )
        ENDIF

C.........  Match hour-specific pol names to inventory pol names
        CALL SETVPTR( NIPOL, NVARS, EINAM, VNAME3D, INDX, IOS )

C.............  Store list of valid pollutant names, or...
C.............  Report which day-specific variables will be ignored
        DO I = 1, NVARS

            J = INDX( I ) 
            IF( J .EQ. 0 ) THEN
                L2 = LEN_TRIM( VNAME3D( I ) )
                MESG = 'WARNING: ' // DESCRIBE( 1:L ) // 
     &                 '-specific pollutant "' //
     &                 VNAME3D( I )( 1:L2 ) // 
     &                 '" is not in inventory file, so ' //
     &                 'it will be ignored.'
                CALL M3MSG2( MESG )

            ELSE

                PNAME( I ) = EINAM( J )

            ENDIF

        ENDDO                    

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END

