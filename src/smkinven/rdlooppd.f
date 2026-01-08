
        SUBROUTINE RDLOOPPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, 
     &                       DAYFLAG, FNAME, SDATE, STIME, NSTEPS, 
     &                       FMTOUT, EASTAT, SPSTAT )

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

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE M3UTILIO

        USE MODINFO, ONLY: NIPPA, BYEAR, CATEGORY, CATLEN

        USE MODLISTS, ONLY: FIREFLAG

        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS
C       CHARACTER(2) CRLF
        INTEGER      GETFLINE
C       INTEGER      JULIAN
C       INTEGER      JUNIT
C       INTEGER      SECSDIFF

C        EXTERNAL     CRLF, GETFLINE, JULIAN, JUNIT, SECSDIFF
        EXTERNAL     GETFLINE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

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
        INTEGER,      INTENT(OUT) :: FMTOUT         ! file format code
        INTEGER,      INTENT(OUT) :: EASTAT( NIPPA )! true: pol/act in data
        INTEGER,      INTENT(OUT) :: SPSTAT( MXSPDAT )! true: special in data

C...........   Local file formats
        INTEGER, ALLOCATABLE, SAVE :: FILFMT( : )  ! file format code

C...........   Character strings of day- or hr-specific list file 
        CHARACTER(300), ALLOCATABLE, SAVE :: LSTSTR( : )

C...........   Other local variables
        INTEGER   I, L                        ! counters and indices

        INTEGER          DD                ! tmp day
        INTEGER, SAVE :: EDATE             ! ending date
        INTEGER, SAVE :: ETIME             ! ending time
        INTEGER          IFIL              ! file no. counter
        INTEGER          IDEV              ! input file unit no.
        INTEGER          IOS               ! i/o status
        INTEGER          IREC              ! record counter
        INTEGER          LP                ! length of period description
        INTEGER          MM                ! tmp month
        INTEGER          MMDD1             ! tmp start month and day 
        INTEGER          MMDD2             ! tmp end month and day 
        INTEGER, SAVE :: NFILE             ! no. files in the list

        LOGICAL, SAVE :: EFLAG = .FALSE.   ! true: error found
        LOGICAL          DFLAG             ! true: retrieve date/time & pol/act 
        LOGICAL, SAVE :: DCALLONE = .TRUE. ! true: 1st time called for day-spec
        LOGICAL, SAVE :: GETRANGE = .TRUE. ! true: get date range if it's available
        LOGICAL, SAVE :: HCALLONE = .TRUE. ! true: 1st time called for hr-spec
        LOGICAL          NFLAG             ! true: retrieve no. sources
        LOGICAL          NEWLOOP           ! true: start of a new read loop

        CHARACTER(4)     PERIOD            ! tmp period name
        CHARACTER(200)   NAMTMP            ! file name buffer
        CHARACTER(300)   MESG              ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDLOOPPD' !  program name

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

C.........  Get number of lines of inventory files in list format
        NFILE = GETFLINE( FDEV, MESG )

C.........  Determine format of day- or hour-specific inputs, and store file
C           names
        IF( ( DAYFLAG .AND. DCALLONE ) .OR. HCALLONE ) THEN

C.............  Allocate memory for storing file formats
            IF( ALLOCATED( FILFMT ) ) DEALLOCATE( FILFMT )
            ALLOCATE( FILFMT( NFILE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FILFMT', PROGNAME )
            FILFMT = 0  ! array

C.............  Allocate memory for storing contents of list-format'd file
            IF( ALLOCATED( LSTSTR ) ) DEALLOCATE( LSTSTR )
            ALLOCATE( LSTSTR( NFILE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LSTSTR', PROGNAME )

C.............  Store lines of day-specific list file
            CALL RDLINES( FDEV, MESG, NFILE, LSTSTR )

C.............  Get file format for files listed in day-specific list file
            CALL CHKLSTFL( NFILE, FNAME, LSTSTR, FILFMT )
        
            IF( DAYFLAG ) THEN
                DCALLONE = .FALSE.
            ELSE
                HCALLONE = .FALSE.
            END IF
            
        END IF

C.........  Get the number of lines in the input file
        MESG = PERIOD(1:LP) // '-specific emissions file'
        NFILE = GETFLINE( FDEV, MESG )

        NEWLOOP = .TRUE.
        IREC = 0
        DO IFIL = 1, NFILE

            NAMTMP = LSTSTR( IFIL )

C.............  If line is LIST header, skip to next line
            IF( INDEX( NAMTMP, 'LIST' ) > 0 ) CYCLE

C.............  If line contains the range of dates, read this packet and
C               skip to next line
            I = INDEX( NAMTMP, 'DATERANGE' )
            IF( I .GT. 0 ) THEN
                I = I + 9
                L = LEN_TRIM( NAMTMP )
                READ( NAMTMP( I:L ), *, IOSTAT=IOS ) MMDD1, MMDD2

C.................  If problem reading dates from string
                IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'ERROR: DATERANGE header ' //
     &                     'misformatted at line', IFIL
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

C.................  Store Julian start and end dates
                IF ( GETRANGE ) THEN
                    MM = MMDD1/100
                    DD = MMDD1 - MM*100
                    SDATE = BYEAR*1000 + JULIAN( BYEAR, MM, DD )

                    MM = MMDD2/100
                    DD = MMDD2 - MM*100
                    EDATE = BYEAR*1000 + JULIAN( BYEAR, MM, DD )
                END IF

                GETRANGE = .FALSE.

                CYCLE

            END IF

C.............  Open files, and report status
            IDEV = JUNIT()
            OPEN( IDEV, ERR=1003, FILE=NAMTMP, STATUS='OLD' )

            WRITE( MESG,94010 ) 'Successful open ' //
     &             'for emissions file:' // CRLF() // BLANK5 //
     &             NAMTMP( 1:LEN_TRIM( NAMTMP ) )
            CALL M3MSG2( MESG )

C.............  Read day-specific or hour-specific file
            IF( FILFMT( IFIL ) .EQ. EMSFMT ) THEN

                CALL RDEMSPD( IDEV, TZONE, OUTSTEP, MXPDSRC, DFLAG, 
     &                        NFLAG, NEWLOOP, DAYFLAG, SDATE, STIME, 
     &                        EDATE, ETIME, EASTAT, SPSTAT )

            ELSE IF( FILFMT( IFIL ) .EQ. FF10FMT .OR.
     &               FILFMT( IFIL ) .EQ. FF10DYFMT .OR.
     &               FILFMT( IFIL ) .EQ. FF10HRFMT ) THEN        
                FIREFLAG = .FALSE.  ! to process ptfire in FF10 format as of SMOKEv3.7
                CALL RDFF10PD( IDEV, TZONE, OUTSTEP, MXPDSRC, DFLAG,
     &                        NFLAG, NEWLOOP, DAYFLAG, SDATE, STIME,
     &                        EDATE, ETIME, EASTAT, SPSTAT )

            ELSE IF ( FILFMT( IFIL ) .EQ. CEMFMT ) THEN

                CALL RDCEMPD( IDEV, TZONE, INSTEP, MXPDSRC, DFLAG,
     &                        NFLAG, NEWLOOP, DAYFLAG, SDATE, STIME,
     &                        EDATE, ETIME, EASTAT, SPSTAT )

            ELSE IF ( FILFMT( IFIL ) .EQ. ORLDYFRFMT ) THEN

                CALL RDORLFR( IDEV, TZONE, OUTSTEP, MXPDSRC, DFLAG,
     &                        NFLAG, NEWLOOP, DAYFLAG, SDATE, STIME,
     &                        EDATE, ETIME, EASTAT, SPSTAT )

            ELSE

                WRITE( MESG,94010 ) 'File format', FILFMT(IFIL), 'of '//
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

C.........  Set output format ( to make sure that any daterange
C           headers get filtered out)
        FMTOUT = MAXVAL( FILFMT ) 

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
