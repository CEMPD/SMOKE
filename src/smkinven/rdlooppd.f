
        SUBROUTINE RDLOOPPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, 
     &                       DAYFLAG, FNAME, SDATE, STIME, NSTEPS, 
     &                       EASTAT )

C***************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine loops through the day-specific or hour-specific input
C      files and reads them. It determines whether these files are in EPS
C      or EMS-95 format, and reads the accordingly
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 12/99 by M. Houyoux
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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

C...........   EXTERNAL FUNCTIONS
        CHARACTER*2  CRLF
        INTEGER      GETFLINE
        INTEGER      JUNIT
        INTEGER      SECSDIFF

        EXTERNAL     CRLF, GETFLINE, JUNIT, SECSDIFF

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  SUBROUTINE ARGUMENTS
        INTEGER,      INTENT (IN) :: FDEV           ! file unit no.
        INTEGER,      INTENT (IN) :: TZONE          ! output time zone
        INTEGER,      INTENT (IN) :: INSTEP         ! expected data time step HHMMSS
        INTEGER,      INTENT (IN) :: OUTSTEP        ! output time step HHMMSS
        INTEGER,      INTENT (IN) :: MXPDSRC        ! max. day- or hr-spec srcs
        LOGICAL,      INTENT (IN) :: DAYFLAG        ! true: day-specific
        CHARACTER(*), INTENT (IN) :: FNAME          ! logical file name
        INTEGER,     INTENT(INOUT):: SDATE          ! Julian start date in TZONE
        INTEGER,     INTENT(INOUT):: STIME          ! data start time in TZONE
        INTEGER,      INTENT(OUT) :: NSTEPS         ! no. time steps
        LOGICAL,      INTENT(OUT) :: EASTAT( NIPPA )! true: pol/act in data

C...........   Character strings of day- or hr-specific list file 
        CHARACTER*300, ALLOCATABLE, SAVE :: NLSTSTR( : )

C...........   Other local variables
        INTEGER, SAVE :: EDATE             ! ending date
        INTEGER, SAVE :: ETIME             ! ending time
        INTEGER, SAVE :: FILFMT            ! file format code
        INTEGER          IFIL              ! file no. counter
        INTEGER          IDEV              ! input file unit no.
        INTEGER          IOS               ! i/o status
        INTEGER          IREC              ! record counter
        INTEGER          LP                ! length of period description
        INTEGER, SAVE :: NFILE             ! no. files in the list

        LOGICAL, SAVE :: EFLAG = .FALSE.   ! true: error found
        LOGICAL          DFLAG             ! true: retrieve date/time & pol/act 
        LOGICAL, SAVE :: DCALLONE = .TRUE. ! true: 1st time called for day-spec
        LOGICAL, SAVE :: HCALLONE = .TRUE. ! true: 1st time called for hr-spec
        LOGICAL          NFLAG             ! true: retrieve no. sources
        LOGICAL          NEWLOOP           ! true: start of a new read loop

        CHARACTER*3      PERIOD            ! tmp period name
        CHARACTER*200    NAMTMP            ! file name buffer
        CHARACTER*300    MESG              ! message buffer

        CHARACTER*16 :: PROGNAME = 'RDLOOPPD' !  program name

C***********************************************************************
C   begin body of program RDLOOPPD

C.........  When start date is zero, set flag for getting dates, times, and
C           pollutants/activities
        IF( SDATE .EQ. 0 ) THEN

            DFLAG = .TRUE.
            NFLAG = .FALSE.

C.........  When start date is non-zero, but maximum count is zero, set flag
C           for getting maximum record count per date.
        ELSE IF( MXPDSRC .EQ. 0 ) THEN

            DFLAG = .FALSE.
            NFLAG = .TRUE.

C.........  Otherwise, set to read all data
        ELSE

            DFLAG = .FALSE.
            NFLAG = .FALSE.

        END IF

C.........  Create label for file name messages
        PERIOD = 'Hour'
        IF( DAYFLAG ) PERIOD = 'Day'
        LP = LEN_TRIM( PERIOD )

C.........  Generate message for GETFLINE and RDLINES calls
        MESG = CATEGORY( 1:CATLEN ) // ' ' // PERIOD( 1:LP ) // 
     &         '-specific inventory file in list format'

C.........  Determine format of day- or hour-specific inputs, and store file
C           names
C.........  Day-specific inputs...
C.........  First time called for day-specific
        IF( DAYFLAG .AND. DCALLONE ) THEN

C.............  Get number of lines of inventory files in list format
            NFILE = GETFLINE( FDEV, MESG )

C.............  Allocate memory for storing contents of list-format'd file
            ALLOCATE( NLSTSTR( NFILE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLSTSTR', PROGNAME )

C.............  Store lines of day-specific list file
            CALL RDLINES( FDEV, MESG, NFILE, NLSTSTR )

C.............  Get file format for files listed in day-specific list file
            CALL CHKLSTFL( NFILE, FNAME, NLSTSTR, FILFMT )

            DCALLONE = .FALSE.

C.........  Hours-specific inputs...
C.........  First time called for hour-specific
        ELSE IF( .NOT. DAYFLAG .AND. HCALLONE ) THEN

C.............  Get number of lines of inventory files in list format
            NFILE = GETFLINE( FDEV, MESG )

C.............  Allocate memory for storing contents of list-format'd file
            IF( ALLOCATED( NLSTSTR ) ) DEALLOCATE( NLSTSTR )
            ALLOCATE( NLSTSTR( NFILE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLSTSTR', PROGNAME )

C.............  Store lines of hour-specific list file
            CALL RDLINES( FDEV, MESG, NFILE, NLSTSTR )

C.............  Get file format for files listed in hour-specific list file
            CALL CHKLSTFL( NFILE, FNAME, NLSTSTR, FILFMT )

            HCALLONE = .FALSE.

        END IF

C.........  Get the number of lines in the input file
        MESG = PERIOD(1:LP) // '-specific emissions file'
        NFILE = GETFLINE( FDEV, MESG )

        NEWLOOP = .TRUE.
        IREC = 0
        DO IFIL = 1, NFILE

            NAMTMP = NLSTSTR( IFIL )

C.............  Open files, and report status
            IDEV = JUNIT()
            OPEN( IDEV, ERR=1003, FILE=NAMTMP, STATUS='OLD' )

            WRITE( MESG,94010 ) 'Successful open ' //
     &             'for emissions file:' // CRLF() // BLANK5 //
     &             NAMTMP( 1:LEN_TRIM( NAMTMP ) )
            CALL M3MSG2( MESG )

C.............  Read EMS-95 day-specific or hour-specific file for EMS-95 format
            IF( FILFMT .EQ. EMSFMT ) THEN

        	CALL RDEMSPD( IDEV, TZONE, OUTSTEP, MXPDSRC, DFLAG, 
     &                        NFLAG, NEWLOOP, DAYFLAG, SDATE, STIME, 
     &                        EDATE, ETIME, EASTAT )

            ELSE IF( FILFMT .EQ. EPSFMT ) THEN

                CALL RDEPSPD( IDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, 
     &                        DFLAG, NFLAG, NEWLOOP, DAYFLAG, 
     &                        SDATE, STIME, EDATE, ETIME, EASTAT )

            ELSE

                WRITE( MESG,94010 ) 'File format', FILFMT, 'of ',
     &                 PERIOD( 1:LP ) // '-specific data not recognized'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

            NEWLOOP = .FALSE.

            CLOSE( IDEV )

        END DO

C.........  Abort if error found
        IF( EFLAG ) THEN

            MESG = 'Problem reading day/hour specific list file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Compute number of time steps based on SDATE and EDATE
        NSTEPS = 1 + SECSDIFF( SDATE, STIME, EDATE, ETIME ) / 
     &               ( OUTSTEP / 10000 * 3600 )

C.........  Rewind input file
        REWIND( FDEV )

        RETURN

999     MESG = 'ERROR: Unexpected end of file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

1003    MESG = 'ERROR: Could not open file ' // CRLF() // BLANK10//
     &         NAMTMP
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDLOOPPD
