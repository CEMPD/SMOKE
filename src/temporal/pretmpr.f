
        SUBROUTINE PRETMPR( MNAME, WNAME, TZONE, TSTEP, SDATE, STIME, 
     &                      NSTEPS, MDATE, MTIME, WDATE, WTIME, 
     &                      WEDATEZ )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     This subroutine compares the simulation start date, start time
C     and duration with the dates and times of the temperature input files
C     and modifies the date, time, and duration accordingly. It gives warnings
C     when the temperature file(s) are not completely consistent with the
C     requested simulations period.
C
C  PRECONDITIONS REQUIRED:
C     Temperature file MNAME open
C     Min/Max temperature file WNAME open
C     Episode start date, time, time zone, and duration initialized
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 1/99 by M. Houyoux
C
C****************************************************************************/
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
C***************************************************************************

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        CHARACTER*14    MMDDYY
        INTEGER         SECSDIFF
        INTEGER         TIME2SEC

        EXTERNAL  CRLF, SECSDIFF, TIME2SEC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT    (IN) :: MNAME  ! gridded tmpr file
        CHARACTER(*), INTENT    (IN) :: WNAME  ! min/max tmpr file
        INTEGER     , INTENT    (IN) :: TZONE  ! output time zone
        INTEGER     , INTENT    (IN) :: TSTEP  ! processing time step
        INTEGER     , INTENT(IN OUT) :: SDATE  ! output start date
        INTEGER     , INTENT(IN OUT) :: STIME  ! output start time
        INTEGER     , INTENT(IN OUT) :: NSTEPS ! output no. time steps
        INTEGER     , INTENT   (OUT) :: MDATE  ! 1st valid date - MNAME file
        INTEGER     , INTENT   (OUT) :: MTIME  ! 1st valid time - MNAME file
        INTEGER     , INTENT   (OUT) :: WDATE  ! 1st valid date - WNAME file
        INTEGER     , INTENT   (OUT) :: WTIME  ! 1st valid time - WNAME file
        INTEGER     , INTENT   (OUT) :: WEDATEZ! last valid date- WNAME file

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
        INTEGER         WETIMEZ      ! WNAME ending time
        INTEGER         WSDATEZ      ! WNAME start date 
        INTEGER         WSTIMEZ      ! WNAME start time

        LOGICAL      :: EFLAG   = .FALSE. ! true: error occurred
        LOGICAL      :: MADJUST = .FALSE. ! true: date/time adjusted for MNAME
        LOGICAL      :: WADJUST = .FALSE. ! true: date/time adjusted for WNAME

        CHARACTER*300   MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'PRETMPR' ! program name

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
        CALL NEXTIME( EDATEZ, ETIMEZ, ( NSTEPS-1 ) * 10000 )

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
        CALL NEXTIME( MEDATEZ, METIMEZ, ( MXREC3D-1 )*10000 )

C.........  Get header information from min/max temperature file
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
        CALL NEXTIME( WEDATEZ, WETIMEZ, ( MXREC3D-1 )*240000 )

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

C.........  Check if simulation start date < min/max temperature start date
C           and adjust simulation start date, if needed.
        IF( SDATEZ .LT. WSDATEZ ) THEN

            WADJUST = .TRUE.
            SDATEZ = WSDATEZ
            EDATEZ = 0        ! min/max files are always one full day at a time

        END IF

C.........  Compute start date/time in output time zone
        SDATENEW = SDATEZ
        STIMENEW = STIMEZ
        CALL NEXTIME( SDATENEW, STIMENEW, ( IZONE - TZONE )*10000 )

C.........  Compute new number of time steps
        ISDIFF    = SECSDIFF( SDATEZ, STIMEZ, EDATEZ, ETIMEZ )
        ISECS     = TIME2SEC( TSTEP )
        NSTEPSNEW = ISDIFF / ISECS + 1

C.........  Write messages stating changes because of meteorology files
c        IF( MADJUST .OR. WADJUST ) THEN
        IF( MADJUST ) THEN

            IF( MADJUST ) THEN
        	MESG = 'WARNING: Adjusting simulation start and/or '//
     &                 'end date/time for '// CRLF()// BLANK10 //
     &                 'file "' // MNAME( 1:ML ) // '" ...'
                CALL M3MSG2( MESG )
            END IF

            IF( WADJUST ) THEN
        	MESG = 'WARNING: Adjusting simulation start and/or '//
     &                 'end date for '// CRLF()// BLANK10 //
     &                 'file "' // WNAME( 1:WL ) // '" ...'
                CALL M3MSG2( MESG )
            END IF

            WRITE( MESG,94010 ) BLANK5 //
     &        'For time zone :', TZONE   , CRLF()// BLANK10 //
     &        'Old start date:', SDATE   , CRLF()// BLANK10 //
     &        'Old start time:', STIME   , CRLF()// BLANK10 //
     &        'Old duration  :', NSTEPS  , CRLF()// BLANK10 //
     &        'New start date:', SDATENEW, CRLF()// BLANK10 //
     &        'New start time:', STIMENEW, CRLF()// BLANK10 //
     &        'New duration  :', NSTEPSNEW

            CALL M3MSG2( MESG )

        END IF

C.........  Now that date/time in met time zone is consistent with met...
C.........  Compare simulation start date and time to that of gridded met
C           data and give a warning if met file starts sooner.
C.........  Reset output met start date and time based on simulation
        IF( SDATEZ .NE. MSDATEZ .OR. STIMEZ .NE. MSTIMEZ ) THEN

            ISDIFF = SECSDIFF( MSDATEZ, MSTIMEZ, SDATEZ, STIMEZ )
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

C.........  Compare simulation start date to that of min/max met data and 
C           give a warning if met file starts sooner
C.........  Reset output min/max start date and time based on simulation
        IF( SDATEZ .NE. WSDATEZ ) THEN

            ISDIFF = SECSDIFF( WSDATEZ, 0, SDATEZ, 0 )
            SECS   = TIME2SEC( TSTEP )
            DYS    = ISDIFF / ( SECS * 24 )
            
            WRITE( MESG,94010 ) 'WARNING: Min/max meteorology file '//
     &             'starts', DYS, 'days before simulation.'
            CALL M3MSG2( MESG )

            WDATE = SDATEZ
            WTIME = 0

C.........  Otherwise, set min/max start date and time based on min/max file
        ELSE
            WDATE = WSDATEZ
            WTIME = 0

        END IF

C.........  Check if simulation end date < min/max temperature ending date and
C           give a warning about how many hours will be affected.
        ISDIFF = SECSDIFF( WEDATEZ, 230000, EDATEZ, ETIMEZ )

        IF( ISDIFF .GT. 0 ) THEN

            HRS = ISDIFF / 3600

            WRITE( MESG,'(A,I3,A)' ) 
     &         'WARNING: The last ', HRS, ' hour(s) of the '//
     &         'simulation will reuse min/max ' // CRLF() // BLANK10 // 
     &         'temperatures from ' // MMDDYY( WSDATEZ )

            CALL M3MSG2( MESG )

        END IF

C.........   Reset simulation date/time/duration to new values
        SDATE  = SDATENEW
        STIME  = STIMENEW
        NSTEPS = NSTEPSNEW

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE PRETMPR
