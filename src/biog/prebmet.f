
        SUBROUTINE PREBMET( MNAME, WNAME, SFILE, TZONE, TSTEP, SDATE, 
     &                      STIME, NSTEPS, MDATE, MTIME, WDATE, WTIME )

C***********************************************************************
C  subroutine body starts at line 104 
C
C  DESCRIPTION:
C       Shifts met data to desired input time zone
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 11/99 by J. Vukovich
C
C****************************************************************************/
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

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)    CRLF
C       INTEGER         SECSDIFF
C       INTEGER         TIME2SEC

C        EXTERNAL  CRLF, SECSDIFF, TIME2SEC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT    (IN) :: MNAME  ! gridded tmpr file
        CHARACTER(*), INTENT    (IN) :: WNAME  ! min/max tmpr file
        LOGICAL     , INTENT    (IN) :: SFILE  ! one or two met files
        INTEGER     , INTENT    (IN) :: TZONE  ! output time zone
        INTEGER     , INTENT    (IN) :: TSTEP  ! processing time step
        INTEGER     , INTENT(IN OUT) :: SDATE  ! output start date
        INTEGER     , INTENT(IN OUT) :: STIME  ! output start time
        INTEGER     , INTENT(IN OUT) :: NSTEPS ! output no. time steps
        INTEGER     , INTENT   (OUT) :: MDATE  ! 1st valid date - MNAME file
        INTEGER     , INTENT   (OUT) :: MTIME  ! 1st valid time - MNAME file
        INTEGER     , INTENT   (OUT) :: WDATE  ! 1st valid date - WNAME file
        INTEGER     , INTENT   (OUT) :: WTIME  ! 1st valid time - WNAME file

C...........   Other local variables
        INTEGER         EDATEZ       ! Sim end date in zone IZONE
        INTEGER         ETIMEZ       ! Sim end time in zone IZONE
        INTEGER         HRS, DYS     ! tmp no. hours & days
        INTEGER         ISDIFF       ! tmp no. seconds difference
        INTEGER         ISECS        ! tmp no. seconds difference
        INTEGER         IZONE        ! time zone of met data
        INTEGER         ML, WL       ! tmp string lengths
        INTEGER         NSTEPSNEW    ! tmp new no. time steps
        INTEGER         SECS         ! tmp value of seconds
        INTEGER         SDATENEW     ! tmp new SDATE
        INTEGER         STIMENEW     ! tmp new STIME
        INTEGER         SDATEZ       ! SDATE in zone IZONE
        INTEGER         STIMEZ       ! STIME in zone IZONE
        INTEGER         MEDATEZ      ! MNAME ending date 
        INTEGER         METIMEZ      ! MNAME ending time
        INTEGER         MSDATEZ      ! MNAME start date 
        INTEGER         MSTIMEZ      ! MNAME start time
        INTEGER         WEDATEZ      ! WNAME ending date 
        INTEGER         WETIMEZ      ! WNAME ending time
        INTEGER         WSDATEZ      ! WNAME start date 
        INTEGER         WSTIMEZ      ! WNAME start time

        LOGICAL      :: EFLAG   = .FALSE. ! true: error occurred
        LOGICAL      :: MADJUST = .FALSE. ! true: date/time adjusted for MNAME
        LOGICAL      :: WADJUST = .FALSE. ! true: date/time adjusted for WNAME

        CHARACTER(300)  MESG     !  message buffer

        CHARACTER(16) :: PROGNAME = 'PREBMET' ! program name

C***********************************************************************
C   begin body of subroutine PRETMPR

C.........  Get lengths of file names
        ML = LEN_TRIM( MNAME )
        WL = LEN_TRIM( WNAME )

C.........  Convert start date and time to time zone of meteorology.
C.........  Met data time zone (IZONE) assume GMT
        IZONE = 0
        SDATEZ = SDATE
        STIMEZ = STIME
        CALL NEXTIME( SDATEZ, STIMEZ, (TZONE - IZONE )*10000 )

C.........  Compute end date and time of simulation (in GMT)
        EDATEZ = SDATEZ
        ETIMEZ = STIMEZ
        CALL NEXTIME( EDATEZ, ETIMEZ, (NSTEPS-1)*10000 )

C.........  Get header information from gridded temperature file
        IF( .NOT. DESC3( MNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             MNAME( 1:ML ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Store start gridded meterology date/time (in GMT)
        MSDATEZ = SDATE3D
        MSTIMEZ = STIME3D

C.........  Calculate end gridded meteorology date/time (in GMT)
        MEDATEZ = MSDATEZ
        METIMEZ = MSTIMEZ
        CALL NEXTIME( MEDATEZ, METIMEZ, (MXREC3D-1)*10000 )

        IF ( .NOT. SFILE ) THEN
C.........  Get header information from meteorology  file
          IF( .NOT. DESC3( WNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             WNAME( 1:WL ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
          END IF

C.........  Store start min/max meterology date/time (in GMT)
          WSDATEZ = SDATE3D
          WSTIMEZ = STIME3D

C.........  Calculate end min/max meteorology date/time (in GMT)
          WEDATEZ = WSDATEZ
          WETIMEZ = WSTIMEZ
          CALL NEXTIME( WEDATEZ, WETIMEZ, (MXREC3D-1)*10000 )

        ENDIF

C.........  Check if simulation start date/time < gridded temperature
C           start date/time.  If so, adjust simulation start date/time to
C           be consistent with the meteorology.

        IF( ( SDATEZ .EQ. MSDATEZ .AND. STIMEZ .LT. MSTIMEZ ) .OR.
     &        SDATEZ .LT. MSDATEZ ) THEN

            MADJUST = .TRUE.
            SDATEZ = MSDATEZ
            STIMEZ = MSTIMEZ

        END IF

C.........  Check if simulation end date/time > gridded temperature ending
C           date/time.  If so, adjust simulation ending date/time to be
C           consistent with the meteorology.
        IF( ( EDATEZ .EQ. MEDATEZ .AND. ETIMEZ .GT. METIMEZ ) .OR.
     &        EDATEZ .GT. MEDATEZ ) THEN

            MADJUST = .TRUE.
            EDATEZ = MEDATEZ
            ETIMEZ = METIMEZ

        END IF

        IF ( .NOT. SFILE ) THEN

C.........  Check if simulation start date/time < gridded temperature
C           start date/time.  If so, adjust simulation start date/time to
C           be consistent with the meteorology.
          IF( ( SDATEZ .EQ. WSDATEZ .AND. STIMEZ .LT. WSTIMEZ ) .OR.
     &         SDATEZ .LT. WSDATEZ ) THEN

            WADJUST = .TRUE.
            SDATEZ = WSDATEZ
            STIMEZ = WSTIMEZ

          END IF

C.........  Check if simulation end date/time > gridded temperature ending
C           date/time.  If so, adjust simulation ending date/time to be
C           consistent with the meteorology.
          IF( ( EDATEZ .EQ. WEDATEZ .AND. ETIMEZ .GT. WETIMEZ ) .OR.
     &        EDATEZ .GT. WEDATEZ ) THEN

            WADJUST = .TRUE.
            EDATEZ = WEDATEZ
            ETIMEZ = WETIMEZ

          END IF

        ENDIF


C.........  Compute start date/time in output time zone
        SDATENEW = SDATEZ
        STIMENEW = STIMEZ
        CALL NEXTIME( SDATENEW, STIMENEW, ( IZONE - TZONE )*10000 )

C.........  Compute new number of time steps
        ISDIFF    = SECSDIFF( SDATEZ, STIMEZ, EDATEZ, ETIMEZ )
        ISECS     = TIME2SEC( TSTEP )
        NSTEPSNEW = ISDIFF / ISECS + 1


C.........  Write messages stating changes because of meteorology files
        IF( MADJUST .OR. WADJUST ) THEN

            IF( MADJUST ) THEN
                MESG = 'WARNING: Adjusting simulation start and/or '//
     &                 'end date/time for '// CRLF()// BLANK10 //
     &                 'file "' // MNAME( 1:ML ) // '" ...'
                CALL M3MSG2( MESG )
            END IF

            IF ( .NOT. SFILE ) THEN
              IF( WADJUST ) THEN
                MESG = 'WARNING: Adjusting simulation start and/or '//
     &                 'end date for '// CRLF()// BLANK10 //
     &                 'file "' // WNAME( 1:WL ) // '" ...'
                CALL M3MSG2( MESG )
              END IF

            ENDIF

            WRITE( MESG,94010 ) BLANK5 //
     &        'For time zone :', TZONE   , CRLF()// BLANK10 //
     &        'Old start date:', SDATE   , CRLF()// BLANK10 //
     &        'Old start time:', STIME   , CRLF()// BLANK10 //
     &        'Old duration  :', NSTEPS  , CRLF()// BLANK10 //
     &        'New start date:', SDATENEW, CRLF()// BLANK10 //
     &        'New start time:', STIMENEW, CRLF()// BLANK10 //
     &        'New duration  :', NSTEPSNEW

            CALL M3MSG2( MESG )

        ENDIF


C.........  Now that date/time in met time zone is consistent with met...
C.........  Compare simulation start date and time to that of gridded met
C           data and give a warning if met file starts sooner.
C.........  Reset output met start date and time based on simulation
        IF( SDATEZ .NE. MSDATEZ .OR. STIMEZ .NE. MSTIMEZ ) THEN

            ISDIFF = SECSDIFF( SDATEZ, STIMEZ, MSDATEZ, MSTIMEZ )
            SECS   = TIME2SEC( TSTEP )
            HRS    = ISDIFF / SECS
            
            WRITE( MESG,94010 ) 'WARNING: Gridded meteorology file '//
     &             'starts', HRS, 'hours before simulation.'
            CALL M3MSG2( MESG )

            MDATE = SDATEZ
            MTIME = STIMEZ

C.........  Otherwise, set met start date and time based on met file
        ELSE
            MDATE = MSDATEZ
            MTIME = MSTIMEZ

        END IF

        IF ( .NOT. SFILE ) THEN 

          IF( SDATEZ .NE. WSDATEZ .OR. STIMEZ .NE. WSTIMEZ ) THEN

            ISDIFF = SECSDIFF( SDATEZ, STIMEZ, WSDATEZ, WSTIMEZ )
            SECS   = TIME2SEC( TSTEP )
            HRS    = ISDIFF / SECS

            WRITE( MESG,94010 ) 'WARNING: Gridded meteorology file '//
     &             'starts', HRS, 'hours before simulation.'
            CALL M3MSG2( MESG )

            WDATE = SDATEZ
            WTIME = STIMEZ

C.........  Otherwise, set met start date and time based on met file
          ELSE
            WDATE = WSDATEZ
            WTIME = WSTIMEZ

          END IF

        ENDIF

C.........   Reset simulation date/time/duration to new values
        SDATE  = SDATENEW
        STIME  = STIMENEW
        NSTEPS = NSTEPSNEW

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE PREBMET
