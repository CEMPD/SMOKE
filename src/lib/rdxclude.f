
        SUBROUTINE RDXCLUDE( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the NHAPEXCLUDE file that contains a list of country/
C     state/county FIPS codes and SCCs that are to be excluded from
C     calculation of the NONHAPVOC pollutants when combining criteria
C     and toxics inventories.
C
C  PRECONDITIONS REQUIRED:
C     File unit FDEV already is opened.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 11/02 by M. Houyoux
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module is for cross reference tables
        USE MODXREF, ONLY: INDXTA, IFIPTA, CSRCTA, CSCCTA

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         CHKINT
        CHARACTER*2     CRLF
        INTEGER         GETFLINE
        INTEGER         STR2INT

        EXTERNAL        CHKINT, CRLF, GETFLINE, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! NHAPEXCLUDE file unit no.
 
C...........   Local parameters
        INTEGER    , PARAMETER :: MXCOL = 2

C...........   Array of input fields
        CHARACTER(LEN=SCCLEN3)  SEGMENT( MXCOL )
  
C...........   Other local variables
        INTEGER         I, N    !  counters and indices

        INTEGER         IFIP    !  temporary FIPS code
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         NLINES  !  number of lines
        INTEGER         NXREF   !  number of valid x-ref entries

        LOGICAL      :: EFLAG = .FALSE.   !  true: error found

        CHARACTER*10           FIPFMT   ! formt to write co/st/cy to string
        CHARACTER*128          LINE     !  line buffer
        CHARACTER*256          MESG     !  message buffer
        CHARACTER(LEN=FIPLEN3) CFIP     !  buffer for CFIPS code
        CHARACTER(LEN=SCCLEN3) TSCC     !  temporary SCC

        CHARACTER*16 :: PROGNAME = 'RDXCLUDE' ! program name

C***********************************************************************
C   begin body of subroutine RDXCLUDE

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'non-HAP exclusions file' )

C.............  Set up formats
        WRITE( FIPFMT, '("(I",I2.2,".",I2.2,")")' ) FIPLEN3, FIPLEN3

C.........  Allocate memory for unsorted data used in all source categories 
        ALLOCATE( INDXTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
        ALLOCATE( IFIPTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IFIPTA', PROGNAME )
        ALLOCATE( CSCCTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
        ALLOCATE( CSRCTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )
        INDXTA = 0   ! array
        IFIPTA = 0   ! array
        CSCCTA = ' ' ! array
        CSRCTA = ' ' ! array

C.........  Set up constants for loop.

C.........  Second pass through file: read lines and store unsorted data for
C           the source category of interest
        IREC   = 0
        N      = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading non-HAP exclusions file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

C.............  Skip comment lines
            IF( LINE( 1:1 ) .EQ. CINVHDR ) CYCLE

C.............  Depending on source category, transfer line to temporary
C               fields.  In cases where blanks are allowed, do not use
C               STR2INT to prevent warning messages.
            CALL PARSLINE( LINE, MXCOL, SEGMENT )

            CFIP = SEGMENT( 1 )
            TSCC = SEGMENT( 2 )

C.............  Make sure that the co/st/cy code is an integer
            IF( .NOT. CHKINT( CFIP ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Country/state/county ' //
     &                 'code is not an integer at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  If this record is in error, go to next iteration
            IF( EFLAG ) CYCLE

C.............  Convert co/st/cy code to an integer
            CFIP = SEGMENT( 1 )
            IFIP = STR2INT( CFIP )

C.............  Convert integer co/st/cy code back to string, including
C               leading zeros
            WRITE( CFIP ,FIPFMT ) IFIP

            N = N + 1
            IF( N .GT. NLINES ) CYCLE  ! Ensure no overflow

C.............  Store case-indpendent fields
            INDXTA ( N ) = N
            IFIPTA ( N ) = IFIP
            CSCCTA ( N ) = SEGMENT( 2 )

C.............  Store sorting criteria for source.
C.............  NOTE - if point sources are added, make sure that
C               TSCC is justified correctly.
            CSRCTA( N ) = CFIP // TSCC

        END DO      ! End of loop on I for reading in NHAPEXCLUDE file

C.........  Reset number of cross-reference entries in case some were dropped
        NXREF = N

C.........  Write errors for problems with input
        IF( NXREF .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No valid NHAPEXLCUDE entries!'
            CALL M3MSG2( MESG )

        ELSEIF( NXREF .GT. NLINES ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'storing non-HAP exclusions file was', NLINES,
     &             CRLF() // BLANK10 // 'but actually needed', NXREF
            CALL M3MSG2( MESG )

        ENDIF

C.......  Check for errors reading XREF file, and abort
        IF( EFLAG ) THEN
            MESG = 'Problem reading non-HAP exclusions file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        CALL M3MSG2( 'Processing non-HAP exclusions file...' )

        CALL SORTIC( NXREF, INDXTA, CSRCTA )

C.........  Group cross-reference data into tables for different groups
        CALL XREFTBL( 'NONHAP', NXREF )

C.........  Deallocate cross-reference sorting arrays
        DEALLOCATE( INDXTA, IFIPTA, CSCCTA, CSRCTA )

C.........  Rewind file
        REWIND( FDEV )

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of non-HAP' // CRLF() // BLANK5 //
     &         'exclusions file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDXCLUDE
