
        SUBROUTINE RDSPROF( FDEV, POLNAM, NPROFMAX, NPROF, NMSPC,
     &                      INPRF, SPECID, MOLEFACT, MASSFACT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the speciation profile, sorts it, and returns
C      the sorted data
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
        INTEGER     , INTENT    (IN) :: NPROFMAX    ! max no. expected profiles
        INTEGER     , INTENT   (OUT) :: NPROF       ! no. profiles found
        INTEGER     , INTENT   (OUT) :: NMSPC       ! no. unique species IDs
        CHARACTER(*), INTENT   (OUT) :: INPRF   ( NPROFMAX ) ! spec prof numbers
        CHARACTER(*), INTENT(IN OUT) :: SPECID  ( NPROFMAX ) ! names of species
        REAL        , INTENT(IN OUT) :: MOLEFACT( NPROFMAX ) ! mole-based facs
        REAL        , INTENT(IN OUT) :: MASSFACT( NPROFMAX ) ! mass-based facs
                                
C.........  Local parameters
        INTEGER, PARAMETER :: MXSEG = 6        ! # of potential line segments

C...........   Local unsorted arrays

        INTEGER        INDXA ( NPROFMAX )   ! sorting index
        REAL           DIVISA( NPROFMAX )   ! unsorted divisors
        REAL           FACTRA( NPROFMAX )   ! unsorted split factors
        REAL           XMFA  ( NPROFMAX )   ! unsorted mass fraction
        CHARACTER*21   INPSPA( NPROFMAX )   ! unsorted profile no. // species ID
        CHARACTER*16   SPCIDA( NPROFMAX )   ! unsorted species IDs
        
C...........   Other arrays
        CHARACTER*20 SEGMENT( MXSEG )             ! Segments of parsed lines

C...........   Other local variables

        INTEGER         I, J, L1, L2, N   ! counters and indices
        INTEGER         IPS       ! counter for records in which POLID = POLNAM
        INTEGER         IOS       ! i/o error status
        INTEGER         IREC      ! record counter

        REAL            DIVISATP          ! tmp divisor
        REAL            FACTRATP          ! tmp split factor
        REAL            XMFATP            ! tmp mass fraction

        LOGICAL      :: DUPFLAG = .FALSE.   ! true: duplicate entries found
        LOGICAL      :: EFLAG   = .FALSE.   ! true: error found
        LOGICAL      :: ZFLAG   = .FALSE.   ! true: divisor of zero found

        CHARACTER*200   LINE              ! read buffer for a line
        CHARACTER*300   MESG              ! text for M3EXIT()

        CHARACTER(LEN=SPNLEN3)  PPRF      ! previous profile code
        CHARACTER(LEN=SPNLEN3)  TMPPRF    ! tmp profile code
        CHARACTER(LEN=IOVLEN3)  POLID     ! tmp pollutant name
        CHARACTER(LEN=IOVLEN3)  SPECNM    ! tmp species name
        CHARACTER(LEN=IOVLEN3)  PSPCNM    ! previous species name

        CHARACTER*16 :: PROGNAME = 'RDSPROF' ! program name
       
C***********************************************************************
C   Begin body of subroutine RDSPROF

C...........  Read in speciation profile
        IREC  = 0
        IPS   = 0
        DO
        
            READ( FDEV, 93000, END=12, IOSTAT=IOS ) LINE
     
            IREC = IREC + 1
             
            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010)
     &              'I/O error', IOS, 'reading speciation profile '//
     &              'file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Separate the line of data into each part
            CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.............  Check for current pollutant of interest
            POLID = ADJUSTL( SEGMENT( 2 ) )           
            IF ( POLID .NE. POLNAM ) CYCLE
            
            IPS = IPS + 1
            
            IF ( IPS .LE. NPROFMAX ) THEN
            
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
            
12      CONTINUE    !End of read on input file
     
        NPROF = IPS

        IF( NPROF .GT. NPROFMAX ) THEN  ! Check for memory overflow

            WRITE( MESG, 94010 )
     &        'INTERNAL ERROR: Number of profiles ' //
     &        'encountered: ', NPROF, CRLF() // BLANK5 //
     &        'Number of profiles expected: ', NPROFMAX

            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )      

        END IF
       
        CALL SORTIC( NPROF, INDXA, INPSPA )    ! Sort on INPSPA
        
        PPRF   = ' '
        PSPCNM = ' '
        N = 0
        DO I = 1, NPROF

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
        DO I = 1, NPROF
            
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

94010  FORMAT( 10( A, :, I8, :, 1X ) )
        
       END SUBROUTINE RDSPROF                                                                            
