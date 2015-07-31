
        SUBROUTINE RDSPDREF( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the SPDREF file that contains cross-reference information
C     for assigning speed profiles to mobile sources
C
C  PRECONDITIONS REQUIRED:
C     File unit FDEV already is opened.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 12/02 by C. Seppanen (based off rdgref.f and rdxclude.f)
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
        USE MODXREF, ONLY: INDXTA, CSRCTA, CSCCTA, ISPDCDA
        
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         BLKORCMT
        LOGICAL         CHKINT
        CHARACTER(2)    CRLF
        INTEGER         GETFLINE
        INTEGER         STR2INT

        EXTERNAL        BLKORCMT, CHKINT, CRLF, GETFLINE, STR2INT
        
C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! SPDREF file unit no.

C...........   Local parameters
        INTEGER    , PARAMETER :: MXCOL = 3

C...........   Array of input fields
        CHARACTER(SCCLEN3)  SEGMENT( MXCOL )
        
C...........   Other local variables
        INTEGER         I, N    !  counters and indices

        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         ISPD    !  tmp spd profile code
        INTEGER         NLINES  !  number of lines
        INTEGER         NXREF   !  number of valid x-ref entries
        
        LOGICAL      :: EFLAG = .FALSE.   !  true: error found
        
        CHARACTER(10)      FIPFMT   !  formt to write co/st/cy to string
        CHARACTER(128)     LINE     !  line buffer
        CHARACTER(256)     MESG     !  message buffer
        CHARACTER(FIPLEN3) CFIP     !  buffer for CFIPS code
        CHARACTER(SCCLEN3) TSCC     !  temporary SCC

        CHARACTER(16) :: PROGNAME = 'RDSPDREF' ! program name

C***********************************************************************
C   begin body of subroutine RDSPDREF

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'Speeds cross-reference file' )

C.............  Set up formats
        WRITE( FIPFMT, '("(I",I2.2,".",I2.2,")")' ) FIPLEN3, FIPLEN3
                
C.........  Allocate memory for unsorted data used in all source categories 
        ALLOCATE( ISPDCDA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS,' IPDCDA', PROGNAME )
        ALLOCATE( INDXTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
        ALLOCATE( CSCCTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
        ALLOCATE( CSRCTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )
        ISPDCDA = 0  ! array
        INDXTA = 0   ! array
        CSCCTA = ' ' ! array
        CSRCTA = ' ' ! array

C.........  Set up constants for loop.

C.........  Second pass through file: read lines and store unsorted data for
C           the source category of interest
        IREC = 0
        N    = 0
        
        DO I = 1, NLINES
        
            READ( FDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

C.............  Check for I/O errors            
            IF( IOS > 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'I/O error', IOS, 
     &              'reading speed cross-reference file at line', IREC 
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Check for end of file            
            IF( IOS < 0 ) THEN
                MESG = 'End of file reached unexpectedly. ' //
     &                 'Check format of speed' // CRLF() // BLANK5 //
     &                 'cross-reference file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE
            
C.............  Parse line into segments
            CALL PARSLINE( LINE, MXCOL, SEGMENT ) 
            
C.............  Make sure that the co/st/cy code is an integer
            IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Country/state/county ' //
     &                 'code is not an integer at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Make sure that SCC is full length
            IF( LEN_TRIM( SEGMENT( 2 ) ) /= SCCLEN3 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: SCC value is not ',
     &                 SCCLEN3, 'digits at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Make sure that the speed profile number is an integer
            IF( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Speed profile code ' //
     &                 'is not an integer at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Check that profile number is not too long
            IF( LEN_TRIM( ADJUSTL( SEGMENT( 3 ) ) ) > SPDLEN3 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Speed profile code ' //
     &              'is longer than', SPDLEN3, 'digits at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Convert speed profile code to an integer
            ISPD = STR2INT( SEGMENT( 3 ) )

C.............  Check that profile number is greater than zero
            IF( ISPD <= 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Speed profile code ' //
     &                 'is less than or equal to zero at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Skip rest of loop if there has been an error
            IF( EFLAG ) CYCLE
            
C.............  Create co/st/cy code string with leading zeros
            WRITE( CFIP,FIPFMT ) STR2INT( SEGMENT( 1 ) )

C.............  Save SCC in string
            TSCC = SEGMENT( 2 )
            
            N = N + 1
            IF( N .GT. NLINES ) CYCLE  ! Ensure no overflow
            
C.............  Store fields
            INDXTA ( N ) = N
            ISPDCDA( N ) = ISPD
            CSCCTA ( N ) = TSCC
            CSRCTA ( N ) = CFIP // TSCC
            
        END DO   ! End of loop reading SPDREF file

C.........  Reset number of cross-reference entries in case some were dropped
        NXREF = N

C.........  Write errors for problems with input
        IF( NXREF .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No valid entries in SPDREF file!'
            CALL M3MSG2( MESG )

        ELSEIF( NXREF .GT. NLINES ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'storing SPDREF file was', NLINES,
     &             CRLF() // BLANK10 // 'but actually needed', NXREF
            CALL M3MSG2( MESG )
            
        ENDIF
        
C.........  Check for errors reading SPDREF file and abort
        IF( EFLAG ) THEN
            MESG = 'Problem reading speed cross-reference file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL M3MSG2( 'Processing speed cross-reference file...' )

C.........  Sort cross-reference table
        CALL SORTIC( NXREF, INDXTA, CSRCTA )

C.........  Group cross-reference data into tables for different groups
        CALL XREFTBL( 'SPEED', NXREF )

C.........  Deallocate cross-reference sorting arrays
        DEALLOCATE( ISPDCDA, INDXTA, CSCCTA, CSRCTA )

C.........  Rewind file
        REWIND( FDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDSPDREF
