
        SUBROUTINE RDHDAYS( FDEV, SDATE, EDATE )

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for and reads the holiday dates.
C      The current file structure permit holidays to be set globally for the
C      entire inventory.  Future versions will permit setting holidays by
C      region codes.  The data file contains the month and day of the holiday.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 9/2000 by M. Houyoux
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
C**************************************************************************

C...........   Modules for public variables
C...........   This module contains the temporal allocation information
        USE MODTMPRL, ONLY: NHOLIDAY, HOLREGN, HOLJDATE, HOLALTDY

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL        CHKINT
        LOGICAL        BLKORCMT
        INTEGER        INDEX1 
        INTEGER        GETFLINE
        INTEGER        JULIAN
        INTEGER        STR2INT
        INTEGER        YEAR4 

        EXTERNAL       CHKINT, INDEX1, GETFLINE, JULIAN, STR2INT, YEAR4,
     &                 BLKORCMT

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: FDEV           ! file unit number
        INTEGER, INTENT (IN) :: SDATE          ! epsiode start Julian date
        INTEGER, INTENT (IN) :: EDATE          ! episode end Julian date

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: INDEXA  ( : )  ! sorting index
        INTEGER, ALLOCATABLE :: HOLREGNA( : )  ! unsorted holiday region
        INTEGER, ALLOCATABLE :: HDJDATEA( : )  ! unsorted holiday Julian date
        INTEGER, ALLOCATABLE :: HDALTDYA( : )  ! unsorted holiday Julian date

C...........   Local arrays
        CHARACTER(9)    CAPDAYS( 7 )          ! names of days of week in caps
        CHARACTER(20)   SEGMENT( 5 )          ! parsed input line

C...........   Local variables
        INTEGER         I, J, N               ! indices and counters

        INTEGER         IOS                   ! i/o status
        INTEGER         IREC                  ! record number
        INTEGER         DAY                   ! tmp day of week
        INTEGER         DD                    ! tmp day of month
        INTEGER         JDATE                 ! Julian date
        INTEGER         NLINES                ! no. lines in input file
        INTEGER         MM                    ! tmp month
        INTEGER         REGN                  ! tmp region code
        INTEGER         YY                    ! tmp year

        LOGICAL      :: EFLAG = .FALSE.       ! true: error found

        CHARACTER(80)   LINE                  ! Read buffer for a line
        CHARACTER(300)  MESG                  ! Message buffer

        CHARACTER(16) :: PROGNAME = 'RDHDAYS'    !  program name

C***********************************************************************
C   Begin body of subroutine RDHDAYS

C.........  Get the number of lines in the holidays file
        NLINES = GETFLINE( FDEV, 'Holidays file' )

C.........  Allocate memory for the unsorted holidays data
        ALLOCATE( INDEXA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )
        ALLOCATE( HOLREGNA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HOLREGNA', PROGNAME )
        ALLOCATE( HDJDATEA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HDJDATEA', PROGNAME )
        ALLOCATE( HDALTDYA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HDALTDYA', PROGNAME )

C.........  Create array of the days of the week in capital letters
        DO I = 1, 7
            CAPDAYS( I ) = DAYS( I )
            CALL UPCASE( CAPDAYS( I ) )
        END DO

C.........  Read the unsorted, unfiltered holidays data
        I = 0
        IREC = 0 
        DO N = 1, NLINES

            READ ( FDEV, 93000, END=99, IOSTAT=IOS ) LINE

            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading holidays '//
     &                'file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse the line into its 2 segments
            CALL PARSLINE( LINE, 5, SEGMENT )

C.............  Convert first column to region code, if possible
            IF( CHKINT( SEGMENT( 1 ) ) ) THEN
                REGN = STR2INT( SEGMENT( 1 ) ) 

C.................  Temporary message for unsupported region codes
                IF( REGN .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Non-zero region codes in holidays '//
     &                     'file are not yet supported.' 
                    CALL M3MESG( MESG )
                END IF

            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad region code at line',
     &                 IREC, 'of holidays file.'
                CALL M3MESG( MESG )
                CYCLE

            END IF

C.............  Convert second column to Julian date from Gregorian date...
C.............  Convert Month
            IF( CHKINT( SEGMENT( 2 ) ) ) THEN
                MM = STR2INT( SEGMENT( 2 ) ) 

            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad month at line',
     &                 IREC, 'of holidays file.'
                CALL M3MESG( MESG )
                CYCLE

            END IF

C.............  Convert day
            IF( CHKINT( SEGMENT( 3 ) ) ) THEN
                DD = STR2INT( SEGMENT( 3 ) ) 

            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad day at line',
     &                 IREC, 'of holidays file.'
                CALL M3MESG( MESG )
                CYCLE

            END IF

C.............  Convert year
            IF( CHKINT( SEGMENT( 4 ) ) ) THEN
                YY = STR2INT( SEGMENT( 4 ) )

C.................  Convert to 4-digit year, if needed
                IF( YY .LT. 100 ) THEN
                    YY = YEAR4( YY )

                ELSEIF( YY .LT. 1970 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Bad year at line',
     &                     IREC, 'of holidays file.'
                    CALL M3MESG( MESG )
                    CYCLE

                END IF 

            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad year at line',
     &                     IREC, 'of holidays file.'
                CALL M3MESG( MESG )
                CYCLE

            END IF

C.............  Determine Julian date
            JDATE = YY * 1000 + JULIAN( YY, MM, DD )

C.............  If date is outside episode range, ignore
c            IF( JDATE .LT. SDATE .OR. JDATE .GT. EDATE ) CYCLE

C.............  Determine alternative day of the week
            CALL UPCASE( SEGMENT( 5 ) )
            DAY = INDEX1( SEGMENT( 5 ), 7, CAPDAYS )

            IF( DAY .LE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad alternative day of ',
     &                 'week at line', IREC, 'of holidays file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Store holiday records that match episode

            I = I + 1
            INDEXA  ( I ) = I
            HOLREGNA( I ) = REGN
            HDJDATEA( I ) = JDATE
            HDALTDYA( I ) = DAY
            
        END DO

        NHOLIDAY = I

C.........  Abort if error found in holidays file
        IF( EFLAG ) THEN
            MESG = 'Problem reading holidays file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory for the sorted holidays data
        ALLOCATE( HOLREGN( NHOLIDAY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HOLREGN', PROGNAME )
        ALLOCATE( HOLJDATE( NHOLIDAY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HOLJDATE', PROGNAME )
        ALLOCATE( HOLALTDY( NHOLIDAY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HOLALTDY', PROGNAME )

C.........  Sort the index of the dates
        CALL SORTI2( NHOLIDAY, INDEXA, HOLREGNA, HDJDATEA )

C.........  Store the sorted holidays
        DO I = 1, NHOLIDAY

            J = INDEXA( I )
            HOLREGN ( I ) = HOLREGNA( J )
            HOLJDATE( I ) = HDJDATEA( J )
            HOLALTDY( I ) = HDALTDYA( J )

        END DO

C.........  Deallocate local memory
        DEALLOCATE( INDEXA, HOLREGNA, HDJDATEA, HDALTDYA )

C.........  Successful completion
        RETURN

C.........  Unexpected end of file
99      MESG = 'INTERNAL ERROR: Unexpected end of holidays file'
        CALL M3MSG2( MESG )

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, 1X, I8, 1X, A, 1X, F10.6, 1X, A )

94030   FORMAT( A, 1X, I6.6, A, 100( ' SSC(', I2.2, '):', F10.6, : ) )

        END SUBROUTINE RDHDAYS
