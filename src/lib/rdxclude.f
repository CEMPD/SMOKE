
        SUBROUTINE RDXCLUDE( FDEV )

C***********************************************************************
C  subroutine body starts at line 98
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
        USE MODXREF, ONLY: INDXTA, CFIPTA, CSRCTA, CSCCTA, NHAP_EXCL

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         BLKORCMT
        LOGICAL         CHKINT
        CHARACTER(2)    CRLF
        INTEGER         GETFLINE
        INTEGER         STR2INT
        LOGICAL         USEEXPGEO

        EXTERNAL        BLKORCMT, CHKINT, CRLF, GETFLINE, STR2INT,
     &                  USEEXPGEO

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! NHAPEXCLUDE file unit no.
 
C...........   Local parameters
        INTEGER    , PARAMETER :: MXCOL = 8

C...........   Array of input fields
        CHARACTER(CHRLEN3)  SEGMENT( MXCOL )
  
C...........   Array of point source plant characeristics
        CHARACTER(CHRLEN3) CHARS( 5 )

C...........   Other local variables
        INTEGER         I, N, L0, L1, L2   !  counters and indices

        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         NLINES  !  number of lines
        INTEGER         NXREF   !  number of valid x-ref entries

        LOGICAL      :: HDRFLAG = .FALSE.
        LOGICAL      :: EFLAG  = .FALSE.      !  true: error found

        CHARACTER(128)     LINE     !  line buffer
        CHARACTER(256)     MESG     !  message buffer
        CHARACTER(FIPLEN3) CFIP     !  buffer for FIPS code
        CHARACTER(PLTLEN3) PLT      !  temporary plant ID
        CHARACTER(SCCLEN3) TSCC     !  temporary SCC
        CHARACTER(ALLLEN3) CSRCALL  !  buffer for source char, incl pol

        CHARACTER(16) :: PROGNAME = 'RDXCLUDE' ! program name

C***********************************************************************
C   begin body of subroutine RDXCLUDE

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'non-HAP inclusion/exclusions file' )

C.........  Allocate memory for unsorted data used in all source categories 
        ALLOCATE( INDXTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
        ALLOCATE( CFIPTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFIPTA', PROGNAME )
        ALLOCATE( CSCCTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
        ALLOCATE( CSRCTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )
        INDXTA = 0   ! array
        CFIPTA = ' ' ! array
        CSCCTA = ' ' ! array
        CSRCTA = ' ' ! array

C.........  Set up constants for loop.

C.........  Second pass through file: read lines and store unsorted data for
C           the source category of interest
        NHAP_EXCL = .TRUE.     ! true: default non-integrate (exclusion)
        IREC   = 0
        N      = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading non-HAP exclusions/inclusions file at line'
     &              , IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Search for header to define either exclusions or inclusion
            L0 = INDEX( LINE, '/EXCLUDE/' )
            L1 = INDEX( LINE, '/INCLUDE/' )
            L2 = INDEX( LINE, '/END/' )
            
            IF( L0 > 0 ) CYCLE         ! process non-HAP exclusions
 
            IF( L1 > 0 ) THEN
                NHAP_EXCL = .FALSE.    ! process non-HAP inclusions
                CYCLE
            ELSE IF( L2 > 0 ) THEN
                HDRFLAG = .TRUE.
                CYCLE 
            END IF

C.............  Skip blank lines
            IF( HDRFLAG ) CYCLE        ! Skip after /END/
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Depending on source category, transfer line to temporary
C               fields.  In cases where blanks are allowed, do not use
C               STR2INT to prevent warning messages.
            CALL PARSLINE( LINE, MXCOL, SEGMENT )

            CFIP = SEGMENT( 1 )
            TSCC = SEGMENT( 2 )

C.............  Smart interpretation of SCC
            CALL FLTRNEG( TSCC )     ! Filter 0 and -9 to blank
            CALL PADZERO( TSCC )     ! Pad with zeros

C.............  Make sure that the co/st/cy code is an integer
            IF( .NOT. USEEXPGEO() .AND.
     &          .NOT. CHKINT( CFIP ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Country/state/county ' //
     &                 'code is not an integer at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  If this record is in error, go to next iteration
            IF( EFLAG ) CYCLE

C.............  Adjust co/st/cy code
            CFIP = ADJUSTR( CFIP )
            CALL PADZERO( CFIP )

            N = N + 1
            IF( N .GT. NLINES ) CYCLE  ! Ensure no overflow

C.............  Store case-indpendent fields
            INDXTA ( N ) = N
            CFIPTA ( N ) = CFIP
            CSCCTA ( N ) = TSCC
            
C.............  Store sorting criteria for source.
C.............  NOTE - if point sources are added, make sure that
C               TSCC is justified correctly.

C.............  For point sources
            IF ( CATEGORY == 'POINT' ) THEN

                PLT      = SEGMENT( 3 )
                CHARS(1) = SEGMENT( 4 )
                CHARS(2) = SEGMENT( 5 )
                CHARS(3) = SEGMENT( 6 )
                CHARS(4) = SEGMENT( 7 )
                CHARS(5) = SEGMENT( 8 )

C.................  Store sorting criteria as right-justified in fields
                CALL BLDCSRC( CFIP, PLT, CHARS(1),
     &                        CHARS(2), CHARS(3), CHARS(4),
     &                        CHARS(5), POLBLNK3, CSRCALL   )

                CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // TSCC 

C.............  For area and mobile sources
            ELSE

                CSRCTA( N ) = CFIP // TSCC

            END IF

        END DO      ! End of loop on I for reading in NHAPEXCLUDE file

C.........  Reset number of cross-reference entries in case some were dropped
        NXREF = N

C.........  Write errors for problems with input
        IF( NXREF .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No valid non-HAP inclusions or exclusions entries!'
            CALL M3MSG2( MESG )

        ELSEIF( NXREF .GT. NLINES ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &          'storing non-HAP inclusions/exclusions file was',NLINES,
     &          CRLF() // BLANK10 // 'but actually needed', NXREF
            CALL M3MSG2( MESG )

        ENDIF

C.......  Check for errors reading XREF file, and abort
        IF( EFLAG ) THEN
            MESG = 'Problem reading non-HAP exclusions/inclusions file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        IF( NHAP_EXCL ) THEN
            MESG = 'Processing non-HAP exclusions file...'
        ELSE
            MESG = 'Processing non-HAP inclusions file...'
        ENDIF

        CALL M3MSG2( MESG )

        CALL SORTIC( NXREF, INDXTA, CSRCTA )

C.........  Group cross-reference data into tables for different groups
        CALL XREFTBL( 'NONHAP', NXREF )

C.........  Deallocate cross-reference sorting arrays
        DEALLOCATE( INDXTA, CFIPTA, CSCCTA, CSRCTA )

C.........  Rewind file
        REWIND( FDEV )

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of non-HAP' // CRLF() // BLANK5 //
     &         'inclusions/exclusions file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDXCLUDE
