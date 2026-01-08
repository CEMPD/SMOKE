
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
C       Created ??/???? by ???
C       09/2025 by HT UNC-IE:  Use M3UTILIO
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
        USE M3UTILIO

C.........  MODULES for public variables
C.........  This module contains the speciation profiles
        USE MODSPRO, ONLY: MXSPFUL, HDRSTART, NSPLST, NSPFUL,
     &                     SPCDEFPOL, SPCDEFLST, NPOLSPRO, SPECID,
     &                     NSPDEF, MXSPLST, INPRF,MOLEFACT, MASSFACT 

        USE MODLISTS, ONLY: MXIDAT, INVDNAM

        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
c       LOGICAL       BLKORCMT
c       CHARACTER(2)  CRLF
c       INTEGER       INDEX1
c       REAL          STR2REAL

c       EXTERNAL      BLKORCMT, CRLF, INDEX1, STR2REAL
        LOGICAL, EXTERNAL :: BLKORCMT

C...........   Subroutine arguments

        INTEGER     , INTENT    (IN) :: FDEV          ! file unit number
        CHARACTER(IOVLEN3), INTENT    (IN) :: POLNAM  ! pol name of interest; HT: change size from * to IOVLEN3
        INTEGER     , INTENT   (OUT) :: NMSPC         ! no. unique species IDs
                                
C.........  Local parameters
        INTEGER    , PARAMETER :: MXSEG = 6         ! # of potential segments

C...........   Local unsorted arrays

        INTEGER        INDXA ( MXSPFUL )   ! sorting index
        REAL           DIVISA( MXSPFUL )   ! unsorted divisors
        REAL           FACTRA( MXSPFUL )   ! unsorted split factors
        REAL           XMFA  ( MXSPFUL )   ! unsorted mass fraction
        CHARACTER(21)  INPSPA( MXSPFUL )   ! unsorted profile no. // species ID
        CHARACTER(IOVLEN3)  SPCIDA( MXSPFUL )   ! unsorted species IDs
        
C...........   Other arrays
        CHARACTER(32) SEGMENT( MXSEG )          ! Segments of parsed lines

C...........   Other local variables

        INTEGER         I, J, L1, L2, N   ! counters and indices
        INTEGER, SAVE:: HDRSTLEN  ! header start field width
        INTEGER         IPS       ! counter for records in which POLID = POLNAM
        INTEGER         IOS       ! i/o error status
        INTEGER         IREC      ! record counter
        INTEGER      :: LASTHDR = 0  ! line number of last header line

        REAL            DIVISATP          ! tmp divisor
        REAL            FACTRATP          ! tmp split factor
        REAL            XMFATP            ! tmp mass fraction

        LOGICAL       :: DUPFLAG  = .FALSE.   ! true: duplicate entries found
        LOGICAL       :: EFLAG    = .FALSE.   ! true: error found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.    ! true: first time routine called
        LOGICAL       :: ZFLAG    = .FALSE.   ! true: divisor of zero found

        CHARACTER(256)  LINE              ! read buffer for a line
        CHARACTER(256)  MESG              ! text for M3EXIT()

        CHARACTER(SPNLEN3)  PPRF      ! previous profile code
        CHARACTER(SPNLEN3)  TMPPRF    ! tmp profile code
        CHARACTER(IOVLEN3)  POLID     ! tmp pollutant name
        CHARACTER(IOVLEN3)  SPECNM    ! tmp species name
        CHARACTER(IOVLEN3)  PSPCNM    ! previous species name

        CHARACTER(16) :: PROGNAME = 'RDSPROF' ! program name
       
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
            IF( BLKORCMT( LINE ) ) CYCLE

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
        NPOLSPRO = 0
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
            MESG = 'WARNING: At least one of the divisors was zero ' //
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

C.............  Local allocatable variables
            INTEGER, ALLOCATABLE, SAVE :: NCNT( : )  ! Count of entries for NONHAP definition
            INTEGER, ALLOCATABLE, SAVE :: NHAPINDX( : ) ! Position of count of NHAP<pol> in list
            LOGICAL, ALLOCATABLE, SAVE :: NHAPFLAG( : ) ! true: flags INVDNAM with NHAP<pol> pollutant
            CHARACTER(10),ALLOCATABLE,SAVE :: NHAPUSPC( : )  ! a list of unique #NHAP species

C.............   Local arrays
            CHARACTER(64) SEGMENT( 3 )          ! Segments of parsed lines

C.............  Local variables
            INTEGER    L, J, K, N, IREC, IOS            

            LOGICAL, SAVE :: FIRSTIME = .TRUE.    ! true: first time routine called
            LOGICAL :: EFLAG    = .FALSE.  ! true: error found
            LOGICAL :: INHEADER = .FALSE.  ! true: in header section

            CHARACTER(10)   J_K              ! read buffer for a line
            CHARACTER(256)  LINE              ! read buffer for a line
            CHARACTER(256)  MESG              ! text for M3EXIT()

C----------------------------------------------------------------------

C.............  Rewind file
            REWIND( FDEV )

C.............  The firsttime this subprogram is called, allocate memory
C               for local counters.
            IF ( FIRSTIME ) THEN 
                ALLOCATE( NCNT( MXIDAT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'NCNT', PROGNAME )  
                ALLOCATE( NHAPINDX( MXIDAT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'NHAPINDX', PROGNAME )  
                ALLOCATE( NHAPFLAG( MXIDAT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'NHAPFLAG', PROGNAME )
                ALLOCATE( NHAPUSPC( MXIDAT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'NHAPUSPC', PROGNAME )
                FIRSTIME = .FALSE.
            END IF

C.............  Initialize local counters (each time subprogram called)
            NCNT     = 0        ! array
            NHAPINDX = 0        ! array
            NHAPFLAG = .FALSE.  ! array
            NHAPUSPC = ' '      ! array

C.............  If status is "SCAN" then Scan for headers, count number 
C               of headers, and get the maximum number of entries per 
C               header.
C.............  If status is "STORE"  then store the info in the allocated
C               header arrays.
            IREC  = 0
            NSPDEF = 0   ! global variable from module
            MXSPLST = 0  ! global variable from module
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

                INHEADER = .FALSE.

C.................  Check if blank or comment line (or could be a header)
                IF( BLKORCMT( LINE ) ) THEN

C.....................  Check for NHAP header line
                    IF( LINE( 1:HDRSTLEN ) .EQ. HDRSTART ) THEN

                        INHEADER = .TRUE.

                        CALL PARSLINE( LINE, 3, SEGMENT )

C........................  Figure out if this is a new NHAP pollutant or one
C                          that we've already had in a previous header line
                        J = INDEX1( SEGMENT( 2 ), MXIDAT, INVDNAM )
                        K = INDEX1( SEGMENT( 3 ), MXIDAT, INVDNAM )

                        IF ( J < 1 .AND. STATUS == 'SCAN' ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Pollutant "'//
     &                        TRIM( SEGMENT(2) )//'" at line', IREC,
     &                        'in Speciation Profile header is not '//
     &                        'in the Inventory Table.'
                            CALL M3MSG2( MESG )
                            CYCLE
                        END IF

C.........................  Skip duplicate #NHAP entries
                        WRITE( J_K, '(2I5)' ) J, K
                        L = INDEX1( J_K, MXIDAT, NHAPUSPC ) 

                        IF ( L > 0 ) THEN
                             CYCLE    ! skip duplicate #NHAP entries
                        ELSE
                             NHAPUSPC( J ) = J_K
                        END IF

C.........................  If pollutant previously found, then set counters
                        IF ( NHAPFLAG( J ) ) THEN
                            NCNT( J ) = NCNT( J ) + 1
                            NSPDEF = NHAPINDX( J )

C........................  If pollutant not previously found, update/initialize counters
                        ELSE
                            NSPDEF = NSPDEF + 1
                            NCNT    ( J ) = 1
                            NHAPINDX( J ) = NSPDEF
                            NHAPFLAG( J ) = .TRUE.
                        END IF

C.....................  Otherwise, skip plain blank or comment
                    ELSE
                        CYCLE
                    END IF      ! End NHAP header or not

                END IF          ! End blank or comment line

C.................  If an actual profile data line, skip line
                IF( .NOT. INHEADER ) CYCLE

C.................  Keep track of maximum number of entries for any pollutant
                MXSPLST = MAXVAL( NCNT )
                NSPDEF  = MAXVAL( NHAPINDX )

                SELECT CASE( STATUS )

C................  If scanning, just go to next interation
                CASE( 'SCAN' )
                    CYCLE

C................  If storing, then store fields
                CASE( 'STORE' )

                    N = NCNT( J )
                    SPCDEFPOL(   NSPDEF ) = SEGMENT( 2 ) 
                    SPCDEFLST( N,NSPDEF ) = SEGMENT( 3 )
                    NSPLST   (   NSPDEF ) = N  ! Final value will be max value

                CASE DEFAULT
                    MESG = 'INTERNAL ERROR: Subprogram ' //
     &                     'READ_GSPRO_HEADER called with unknown '//
     &                     'STATUS argument "'// TRIM( STATUS ) // '"'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END SELECT

            END DO

500         CONTINUE  ! continue after ending read loop
       
C.............  Abort if error
            IF( EFLAG ) THEN
                MESG = 'Problem reading speciation profile file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Rewind input file
            REWIND( FDEV )

            RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000       FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE READ_GSPRO_HEADER

        END SUBROUTINE RDSPROF                                                                            
