
        SUBROUTINE RDSPROF( FDEV, POLNAM, NMSPC )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the speciation profile, sorts it, and returns
C      the sorted data. Also checks for the optional NONHAP<pollutant>
C      headers, allocates memory for these, and reads these if
C      they exist.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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

C.........  MODULES for public variables
C.........  This module contains the speciation profiles
        USE MODSPRO

        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2   CRLF
        REAL          STR2REAL

        EXTERNAL      CRLF, STR2REAL

C...........   Subroutine arguments

        INTEGER     , INTENT    (IN) :: FDEV        ! file unit number
        CHARACTER(*), INTENT    (IN) :: POLNAM      ! pol name of interest
        INTEGER     , INTENT   (OUT) :: NMSPC       ! no. unique species IDs
                                
C.........  Local parameters
        INTEGER    , PARAMETER :: MXSEG    = 6         ! # of potential segments

C...........   Local unsorted arrays

        INTEGER        INDXA ( MXSPFUL )   ! sorting index
        REAL           DIVISA( MXSPFUL )   ! unsorted divisors
        REAL           FACTRA( MXSPFUL )   ! unsorted split factors
        REAL           XMFA  ( MXSPFUL )   ! unsorted mass fraction
        CHARACTER*21   INPSPA( MXSPFUL )   ! unsorted profile no. // species ID
        CHARACTER*16   SPCIDA( MXSPFUL )   ! unsorted species IDs
        
C...........   Other arrays
        CHARACTER*32 SEGMENT( MXSEG )          ! Segments of parsed lines

C...........   Other local variables

        INTEGER         I, J, L1, L2, N   ! counters and indices
        INTEGER, SAVE:: HDRSTLEN  ! header start field width
        INTEGER         IPS       ! counter for records in which POLID = POLNAM
        INTEGER         IOS       ! i/o error status
        INTEGER         IREC      ! record counter
        INTEGER         LASTHDR   ! line number of last header line

        REAL            DIVISATP          ! tmp divisor
        REAL            FACTRATP          ! tmp split factor
        REAL            XMFATP            ! tmp mass fraction

        LOGICAL       :: DUPFLAG  = .FALSE.   ! true: duplicate entries found
        LOGICAL       :: EFLAG    = .FALSE.   ! true: error found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.    ! true: first time routine called
        LOGICAL       :: ZFLAG    = .FALSE.   ! true: divisor of zero found

        CHARACTER*256   LINE              ! read buffer for a line
        CHARACTER*256   MESG              ! text for M3EXIT()

        CHARACTER(LEN=SPNLEN3)  PPRF      ! previous profile code
        CHARACTER(LEN=SPNLEN3)  TMPPRF    ! tmp profile code
        CHARACTER(LEN=IOVLEN3)  POLID     ! tmp pollutant name
        CHARACTER(LEN=IOVLEN3)  SPECNM    ! tmp species name
        CHARACTER(LEN=IOVLEN3)  PSPCNM    ! previous species name

        CHARACTER*16 :: PROGNAME = 'RDSPROF' ! program name
       
C***********************************************************************
C   Begin body of subroutine RDSPROF

C...........  If firstime the routine is called...
        IF( FIRSTIME ) THEN

            HDRSTLEN = LEN_TRIM( HDRSTART )

            CALL READ_GSPRO_HEADER( FDEV, 'SCAN' )
    
C.............  Allocate memory for arrays, even if no headers found
            ALLOCATE( NSPLST( NSPDEF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NSPLST', PROGNAME )  
            ALLOCATE( SPCDEFPOL( NSPDEF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPCDEFPOL', PROGNAME )  
            ALLOCATE( SPCDEFLST( MXSPLST, NSPDEF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPCDEFLST', PROGNAME )  
            NSPLST    = 0    ! array
            SPCDEFPOL = ' '  ! array
            SPCDEFLST = ' '  ! array

            CALL READ_GSPRO_HEADER( FDEV, 'STORE' )

            FIRSTIME = .FALSE.

        END IF

C...........  Read in speciation profile
        IREC  = 0
        IPS   = 0
        SEGMENT = ' '  ! array
        DO
        
            READ( FDEV, 93000, END=120, IOSTAT=IOS ) LINE
            IREC = IREC + 1
             
            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010)
     &              'I/O error', IOS, 'reading speciation profile '//
     &              'file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  If we haven't already passed the headers, skip all header
C               lines
            IF( IREC .LE. LASTHDR ) CYCLE

C.............  Skip blank and comment lines
            IF( LINE .EQ. ' ' ) CYCLE
            IF( LINE(1:1) .EQ. CINVHDR ) CYCLE

C.............  Separate the line of data into each part
            CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.............  Check for current pollutant of interest
            POLID = ADJUSTL( SEGMENT( 2 ) )           
            IF ( POLID .NE. POLNAM ) CYCLE
            
            IPS = IPS + 1
            
            IF ( IPS .LE. MXSPFUL ) THEN
            
                INDXA ( IPS ) = IPS
                TMPPRF        = ADJUSTL ( SEGMENT( 1 ) ) 
                SPECNM        = ADJUSTL ( SEGMENT( 3 ) )
                INPSPA( IPS ) = ADJUSTR ( TMPPRF ) // SPECNM
                SPCIDA( IPS ) = SPECNM
                FACTRA( IPS ) = STR2REAL( SEGMENT( 4 ) )
                DIVISA( IPS ) = STR2REAL( SEGMENT( 5 ) )
                XMFA  ( IPS ) = STR2REAL( SEGMENT( 6 ) )

            END IF

        END DO 
            
120     CONTINUE    ! End of read on input file
     
        NSPFUL = IPS

        IF( NSPFUL .GT. MXSPFUL ) THEN  ! Check for memory overflow

            WRITE( MESG, 94010 )
     &        'INTERNAL ERROR: Number of profiles ' //
     &        'encountered: ', NSPFUL, CRLF() // BLANK5 //
     &        'Number of profiles expected: ', MXSPFUL

            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )      

        END IF
       
        CALL SORTIC( NSPFUL, INDXA, INPSPA )    ! Sort on INPSPA
        
        PPRF   = ' '
        PSPCNM = ' '
        N = 0
        DO I = 1, NSPFUL

            J = INDXA( I )        
            IF ( DIVISA( J ) .EQ. 0 ) THEN
                ZFLAG = .TRUE.
                CYCLE
            END IF
        
            J = INDXA( I )
            
            TMPPRF = INPSPA( J )( 1:SPNLEN3 )
            SPECNM = SPCIDA( J )

C.............  Make sure duplicates are not used
            IF( TMPPRF .EQ. PPRF   .AND. 
     &          SPECNM .EQ. PSPCNM       ) THEN

                DUPFLAG = .TRUE.
                EFLAG   = .TRUE.
                L1 = LEN_TRIM( TMPPRF )
                L2 = LEN_TRIM( SPECNM )
                MESG = 'ERROR: Duplicate entries in speciation ' //
     &                 'profiles file for profile ' // CRLF() //
     &                 BLANK10 // TMPPRF( 1:L1 ) // ', species ' //
     &                 SPECNM( 1:L2 ) // '.'
                CALL M3MESG( MESG ) 

            ELSE
                N = N + 1

                INPRF   ( N ) = TMPPRF
                SPECID  ( N ) = SPECNM
                MOLEFACT( N ) = TON2GM * FACTRA( J ) / DIVISA( J )
                MASSFACT( N ) = TON2GM * XMFA( J )

            END IF

C...........  Count profiles used in this call of subroutine
            IF( TMPPRF .NE. PPRF ) NPOLSPRO = NPOLSPRO + 1

            PPRF   = TMPPRF
            PSPCNM = SPECNM

        END DO

        IF( DUPFLAG ) THEN
            MESG = 'ERROR: Duplicate speciation profile entries ' //
     &             'found. ' //CRLF()// BLANK10 // 
     &             'Remove duplicate entries and try again.'
            CALL M3MSG2( MESG )
        END IF
        
        IF( ZFLAG ) THEN
            MESG = 'ERROR: At least one of the divisors was zero ' //
     &             'in the speciation profiles.' // CRLF()// BLANK10 // 
     &             'Correct column 5 of input file and try again.'
            CALL M3MSG2( MESG )
        END IF

        IF( EFLAG ) THEN
            MESG = 'Problem(s) found in speciation profiles.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
 
C........... Sort the species names in order to loop through and count
C            the number of times the name changes. This will give us the
C            number of unique species names, NMSPC 
        
        CALL SORTIC( IPS, INDXA, SPECID )       ! Sort on SPECID
        
        SPECNM = EMCMISS3  ! Initialize temporary 
        NMSPC   = 0
        DO I = 1, NSPFUL
            
            J = INDXA ( I )
            IF ( SPECID( J ) .NE. SPECNM ) THEN
                NMSPC = NMSPC + 1
                SPECNM = SPECID( J )
            END IF
            
        END DO
        
C......... Rewind file

        REWIND( FDEV )
       
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS
        
C.............  This internal subprogram reads the header of the speciation
C               profile file
            SUBROUTINE READ_GSPRO_HEADER( FDEV, STATUS )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: FDEV   ! unit number
            CHARACTER(*), INTENT (IN) :: STATUS ! Read status 'SCAN' or 'STORE'

C.............  Local variables
            INTEGER    L
            INTEGER :: NCNT = 0

            LOGICAL :: EFLAG    = .FALSE.  ! true: error found
            LOGICAL :: INHEADER = .FALSE.  ! true: in header section
            LOGICAL :: ONSTART  = .FALSE.  ! true: on header start line

C----------------------------------------------------------------------

C.............  Rewind file
            REWIND( FDEV )

C.............  If status is "SCAN" then Scan for headers, count number 
C               of headers, and get the maximum number of entries per 
C               header.
C.............  If status is "STORE"  then store the info in the allocated
C               header arrays.
            IREC  = 0
            NSPDEF = 0
            MXSPLST = 0
            DO
        
                READ( FDEV, 93000, END=500, IOSTAT=IOS ) LINE
                IREC = IREC + 1

C.................  Check error status
                IF ( IOS .GT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010)
     &                  'I/O error', IOS, 'reading speciation profile'//
     &                  ' file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Skip blank and comment lines
                IF( LINE .EQ. ' ' ) CYCLE
                IF( LINE(1:1) .EQ. CINVHDR ) CYCLE

                ONSTART = .FALSE.

C.................  Check for header start
                L = INDEX( LINE, HDRSTART ) 
                IF( L .GT. 0 ) THEN
                    INHEADER = .TRUE.
                    NCNT = 0
                    NSPDEF = NSPDEF + 1
                    ONSTART = .TRUE.
                END IF

C.................  Check for end of header
                L = INDEX( LINE, HDREND ) 
                IF( L .GT. 0 ) THEN
                    INHEADER = .FALSE.
                    LASTHDR = IREC
                    CYCLE
                END IF

C.................  Checking on header end...
C.................  Separate the line of data into each part
                CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.................  If reach the first nonheader line...
                IF( SEGMENT( 2 ) .NE. ' ' ) THEN

C....................  Give warning for automatically setting the end of
C                      the last header
                    IF( INHEADER ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Setting end of '//
     &                         'header in speciation profile file at '//
     &                         'line', IREC - 1
                        INHEADER = .FALSE.
                        LASTHDR = IREC - 1
                    END IF

C....................  Exit from loop because end of header section
                    EXIT

                END IF

C.................  Keep track of maximum number of entries
                IF( INHEADER .AND. .NOT. ONSTART ) THEN
                    NCNT = NCNT + 1
                    MXSPLST = MAX( MXSPLST, NCNT )
                END IF

                SELECT CASE( STATUS )

C................  If scanning, just go to next interation
                CASE( 'SCAN' )
                    CYCLE

C................  If storing, then store fields
                CASE( 'STORE' )

                    LINE = ADJUSTL( LINE )

C....................  Store name of pollutant
                    IF( ONSTART ) THEN
                        L = LEN_TRIM( LINE )

C.......................  Check if bad header
                        IF( L .LE. HDRSTLEN ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 
     &                        'ERROR: Header does not specific a ' //
     &                        'pollutant in profiles at line', IREC
                            CALL M3MSG2( MESG )
                            CYCLE
                        END IF

                        SPCDEFPOL( NSPDEF ) = LINE( 2:L-1 )

C....................  Store definition of that pollutant and count
                    ELSE
                        NSPLST( NSPDEF ) = NCNT
                        SPCDEFLST( NCNT,NSPDEF ) = TRIM( LINE )

                    END IF

                CASE DEFAULT
                    MESG = 'INTERNAL ERROR: Subprogram ' //
     &                     'READ_GSPRO_HEADER called with unknown '//
     &                     'STATUS argument "'// TRIM( STATUS ) // '"'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END SELECT

                CYCLE

            END DO

C.............  Abort if error
            IF( EFLAG ) THEN
                MESG = 'Problem reading speciation profile file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Rewind input file
            REWIND( FDEV )

            RETURN

500         MESG = 'Unexpected end of file prior to end of header'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
       
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000       FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE READ_GSPRO_HEADER

        END SUBROUTINE RDSPROF                                                                            
