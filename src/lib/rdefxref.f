
        SUBROUTINE RDEFXREF( FDEV, PROCESS ) 

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:  
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Copied from rdmplist.F 6/99 by M Houyoux
C
C***********************************************************************
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module is for mobile-specific data
        USE MODMOBIL

C...........   This module is for cross reference tables
        USE MODXREF

C.........  This module contains emission factor tables and related
        USE MODEMFAC

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         CHKINT
        CHARACTER*2     CRLF
        INTEGER         FIND1
        INTEGER         FINDC
        INTEGER         STR2INT
        INTEGER         GETFLINE

        EXTERNAL  CHKINT, CRLF, FIND1, STR2INT, GETFLINE


C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV       ! x-ref file unit no.
        LOGICAL, INTENT (IN) :: PROCESS    ! true: call XREFTBL

C...........   LOCAL PARAMETERS and their descriptions

        INTEGER, PARAMETER :: MXFLDS = 29  ! Max. no. columns in MPLIST
        INTEGER, PARAMETER :: MXPSIVAL = 24 * 3600

C...........   Sorted activities
        INTEGER                   INDXA  ( NIACT ) !  sorting index for actvtys
        CHARACTER(LEN=IOVLEN3) :: SRTACT ( NIACT ) !  sorted activity names

C...........   Array for reading temporal x-ref fields
        CHARACTER*20            SEGMENT( MXFLDS )

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: IDXPSIA( : )  ! 
        INTEGER, ALLOCATABLE :: ALLPSIA( : )
        INTEGER, ALLOCATABLE :: ALLPOSA( : )

C...........    LOCAL VARIABLES and their descriptions:

        INTEGER   C, I, J, K, L, M, N, P   ! indices and counters

        INTEGER   IDUM      ! dummy integer
        INTEGER   IOS       ! i/o status
        INTEGER   IREC      ! record counter
        INTEGER   JACT      ! index to unsorted activities
        INTEGER   LPOS      ! activity position from previous iteration
        INTEGER   LPSI      ! PSI from previous iteration
        INTEGER   NALLPSI   ! no. all PSIs for creating unique list
        INTEGER   NLINES    ! no. lines in xref file
        INTEGER   NREF      ! no. of PSIs in a segment of LINE
        INTEGER   NCOL      ! no. of column in LINE
        INTEGER   NP        ! no. PSIs in group in xref file
        INTEGER   NXREF     ! no. of 
        INTEGER   POS       ! tmp position of activity
        INTEGER   PSI       ! tmp PSI
        INTEGER   PSIARR( 24 ) ! tmp parameter scheme indices per xref entry
        INTEGER   RWT       ! tmp roadway type

        LOGICAL	 :: EFLAG = .FALSE.  ! true: error found
        LOGICAL	 :: AFLAG = .FALSE.  ! true: activity-specific entry was skipped
        LOGICAL	 :: SKIPREC = .FALSE.! true: skip current x-ref entry

        CHARACTER*5     CPOS      ! tmp sorted position of pol
        CHARACTER*10    RWTFMT    ! format for writing roadway type to string
        CHARACTER*300	LINE      ! line of input from MPLIST
        CHARACTER*300	MESG      ! message buffer

        CHARACTER(LEN=IOVLEN3) ACT      !  temporary activity name
        CHARACTER(LEN=SICLEN3) CDUM     !  dummy character field for SIC
        CHARACTER(LEN=LNKLEN3) CLNK     !  temporary link code
        CHARACTER(LEN=ALLLEN3) CSRCALL  !  buffer for source char, incl pol
        CHARACTER(LEN=FIPLEN3) CFIP     !  buffer for CFIPS code
        CHARACTER(LEN=FIPLEN3) FIPZERO  !  buffer for zero FIPS code
        CHARACTER(LEN=RWTLEN3) CRWT     !  buffer for roadway type
        CHARACTER(LEN=VIDLEN3) CVID     !  buffer for vehicle type number
        CHARACTER(LEN=VIDLEN3) VIDHOLD  !  buffer to get around XREFTBL check
        CHARACTER(LEN=VIDLEN3) VIDZERO  !  buffer for zero vehicle type number
        CHARACTER(LEN=SCCLEN3) TSCC     !  temporary SCC
        CHARACTER(LEN=SCCLEN3) SCCBLNK  !  blank SCC
        CHARACTER(LEN=SCCLEN3) SCCZERO  !  buffer for zero SCC
        CHARACTER(LEN=LNKLEN3) LNKZERO  !  buffer for zero link

        CHARACTER*16 :: PROGNAME = 'RDEFXREF' ! program name

C***********************************************************************
C   begin body of subroutine  RDEFXREF

C.........  Set up zero strings for FIPS code of zero and SCC code of zero
        FIPZERO = REPEAT( '0', FIPLEN3 )
        SCCZERO = REPEAT( '0', SCCLEN3 )
        SCCBLNK = REPEAT( ' ', SCCLEN3 )
        LNKZERO = REPEAT( '0', LNKLEN3 )
        VIDHOLD = REPEAT( '-', VIDLEN3 )
        VIDZERO = REPEAT( '0', VIDLEN3 )

C.........  Set up roadway type format
        WRITE( RWTFMT, '("(I",I2.2,")")' ) RWTLEN3

C.........  Sort the actual list of activity names and store it
        DO I = 1, NIACT
            INDXA( I ) = I
        ENDDO

        CALL SORTIC( NIACT, INDXA, ACTVTY )

        DO I = 1, NIACT
            J = INDXA( I )
            SRTACT( I ) = ACTVTY( J )
        END DO

        MESG = 'Reading emission factors cross-reference file...'
        CALL M3MSG2( MESG )

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 
     &                     'Emission factors cross reference file' )

C.........   First pass through file.  Count the number of lines matching 
C            pollutant list.  Do this so that we know how much memory to
C            allocate for the unsorted, unprocessed arrays.

        NREF   = 0
        IREC   = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END = 999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

C.............  Check for line read errors

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading emission factor xref file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

C.............  Parse the line of the file 
            CALL PARSLINE( LINE, MXFLDS, SEGMENT )

C.............  Search for activity name in list and skip entry if it does not
C               match.
            ACT = SEGMENT( 5 )

            K = FINDC( ACT, NIACT, SRTACT )

            IF( K .LE. 0 ) CYCLE

            NREF = NREF + 1

        END DO

C.........  Reset file pointer to start of file
        REWIND( FDEV )

C.........  Allocate memory for unsorted emission factor xref table
        ALLOCATE( INDXTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
        ALLOCATE( ISPTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPTA', PROGNAME )
        ALLOCATE( CSCCTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
        ALLOCATE( CSRCTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )
        ALLOCATE( IPSIA( NREF,24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPSIA', PROGNAME )

C.........  Allocate memory for arrays to use in creating PSILIST
        ALLOCATE( ALLPSIA( NREF*24 ),STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLPSIA', PROGNAME )
        ALLOCATE( IDXPSIA( NREF*24 ),STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXPSIA', PROGNAME )
        ALLOCATE( ALLPOSA( NREF*24 ),STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLPOSA', PROGNAME )
        ALLOCATE( NPSI( NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NPSI', PROGNAME )
       
C.........  Second pass through file: read lines and store unsorted data for
C           the source category of interest
        IREC = 0
        C    = 0
        M    = 0
        CDUM = ' '
        DO I = 1, NLINES

            READ( FDEV, 93000, END = 999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

C.............  Check for line read errors

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading emission factor xref file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

            SEGMENT = ' '  ! array

C.............  Parse the line of the file 
            CALL PARSLINE( LINE, MXFLDS, SEGMENT )

C.............  Make sure that no delimiters were surrounding an asterisk
C.............  Count the actual number of columns in the line
            NCOL = 0
            DO J = MXFLDS, 1, -1
                L = LEN_TRIM( SEGMENT( J ) )
 
C.................  If field is blank, skip it 
                IF( L .EQ. 0 ) CYCLE

C.................  Otherwise, process entry
                IF( SEGMENT( J )( 1:1 ) .EQ. '*' .OR.
     &              SEGMENT( J )( L:L ) .EQ. '*'      ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: delimiters used '//
     &                     'around asterisks in emission factor xref'//
     &                     CRLF() // BLANK10 // 'file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE

                END IF

                IF( NCOL .EQ. 0 ) NCOL = J

            END DO

C.............  Process segments of LINE into known parts
C.............  For now this is configured for mobile source only

            CFIP = SEGMENT( 1 )
            CRWT = ADJUSTL( SEGMENT( 2 ) )
            CLNK = SEGMENT( 3 )
            CVID = SEGMENT( 4 )
            ACT  = SEGMENT( 5 )

C.............  Filter -9 and zeros for roadway type, veh type, link ID
            CALL FLTRNEG( CRWT )
            CALL FLTRNEG( CLNK )
            CALL FLTRNEG( CVID )
            CALL PADZERO( CVID )

C.............  Convert roadway type to integer
            RWT = 0
            IF( CRWT .NE. ' ' ) THEN

                IF( .NOT. CHKINT( CRWT ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: non-number found ' //
     &                     'for roadway type in emission factor xref' //
     &                     CRLF() // BLANK10 // 'file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE

                ELSE
                    RWT = STR2INT( CRWT )
 
                END IF
            END IF

C.............  Make sure roadway type is not road class
            IF( RWT .NE. 0 ) THEN
                K = FIND1( RWT, NRCLAS, AMSRDCLS )
                IF( K .GT. 0 ) THEN
                    RWT = RDWAYTYP( K )
                    WRITE( CRWT, RWTFMT ) RWT
                END IF 
            ELSE
                CRWT = REPEAT( '0', RWTLEN3 )
 
            END IF 

C.............  Build internal SCC (can't use MBSCCADJ because input file
C               does not use SCCs). Use VIDHOLD to work around check in 
C               Xreftbl for no SCC when more than road class is given.  This
C               will allow use of Link-specific without vehicle type.
            IF( RWT  .NE. 0       .AND. 
     &          CVID .EQ. VIDZERO .AND. 
     &          CLNK .NE. ' '           ) THEN
                CVID = VIDHOLD
            END IF
            TSCC = CRWT // CVID
            CALL PADZERO( TSCC )

C.................  Post-process country/state/county code for blank and
C                   missing values
            CALL FLTRNEG( CFIP )
            CALL PADZERO( CFIP )

C.................  Compare activity name with master list.
            JACT = FINDC( ACT, NIACT, SRTACT )

            IF( JACT .LE. 0 ) CYCLE

C.................  Skip lines that are not valid for this inven and src cat
            IF( SKIPREC ) CYCLE

            N = 0 
            DO J = 6, NCOL  !  head of segment processing loop

                P = INDEX( SEGMENT( J ), '*' )

C.................  For only one PSI given
                IF( P .LE. 0 ) THEN 
             
                    N        = N + 1
                    PSIARR( N ) = STR2INT( SEGMENT( J ) )

C.................  For multiple PSIs given
                ELSE

                    NP = STR2INT( SEGMENT( J )( 1:P-1 ) )

                    DO K = 1, NP
                        N = N + 1
                        L = LEN_TRIM( SEGMENT( J ) )
                        PSIARR( N ) = STR2INT( SEGMENT( J )( P+1:L ) )
                    END DO

                END IF

            END DO

C.............  Check for too many or too few records in line
            IF( N .NE. 24 ) THEN  

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &               'ERROR: Need PSIs for 24 hours at line', IREC,
     &               'of EF xref (', N, 'found now).'
                CALL M3MESG( MESG )
                CYCLE         ! To head of file read loop

            END IF

C.................  Write activity name position to a character string
            WRITE( CPOS, '(I5)' ) JACT  

            C = C + 1

            IF ( C .LE. NREF ) THEN

C.................  If link ID is not set, then set road class to blank
                IF( CLNK .EQ. ' ' ) CRWT = ' '

C.................  Store sorting criteria as right-justified in fields
C.................  For mobile, we are using "CSCCTA" to store the roadway
C                   type so that we can use the XREFTBL subroutine.  
                CSRCALL = ' '
                CALL BLDCSRC( CFIP, CRWT, CLNK, CHRBLNK3, CHRBLNK3, 
     &                        CHRBLNK3, CHRBLNK3, POLBLNK3, CSRCALL )

                INDXTA( C ) = C
                ISPTA ( C ) = JACT
                CSCCTA( C ) = TSCC
                CSRCTA( C ) = CSRCALL( 1:SC_ENDP( NCHARS ) ) // 
     &                        TSCC // CPOS
                                
                DO J = 1, 24

                    M = M + 1
                    PSI            = PSIARR( J )
                    IPSIA  ( C,J ) = PSI
                    IDXPSIA( M )   = M
                    ALLPSIA( M )   = PSIARR( J )
                    ALLPOSA( M )   = JACT
 
                    IF( PSI .GT. MXPSIVAL ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 ) 
     &                         'Value of PSIs cannot exceed ', MXPSIVAL
                        CALL M3MSG2( MESG )
                    END IF

                END DO
   
            END IF               !  end dimensioning check

        END DO

C.........  Set actual number of cross-reference entries and total number of
C           all PSIs in file
        NXREF   = C
        NALLPSI = M

        IF( NXREF .GT. NREF ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &           'INTERNAL ERROR: EF xref table dimensioned ', NREF,
     &           'but needed ', NXREF
            CALL M3MSG2( MESG )

        ELSE IF( NXREF .EQ. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &           'ERROR: Emission factors xref table had 0 entries '//
     &           'that matched inventory'
            CALL M3MSG2( MESG )
 
        ENDIF

        IF ( EFLAG ) THEN
            MESG = 'Problem reading emission factor ' //
     &             'cross-reference file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MESG = 'Processing emission factor cross-reference file...'
        CALL M3MSG2( MESG )

C.........  Generate list of unique PSIs in the inventory
        CALL SORTI2( NALLPSI, IDXPSIA, ALLPOSA, ALLPSIA )

C.........  Count up actual number of PSIs
        NPSI = 0  ! array
        LPOS = -9
        DO I = 1, NALLPSI

            J   = IDXPSIA( I )
            POS = ALLPOSA( J )
            PSI = ALLPSIA( J )

            IF( POS .NE. LPOS ) LPSI = -9

            IF( PSI .NE. LPSI ) NPSI( POS ) = NPSI( POS ) + 1

            LPOS = POS
            LPSI = PSI

        END DO
        MXXACTV = POS
        MXXNPSI = MAXVAL( NPSI )

C.........  Allocate memory for unique PSI list
        ALLOCATE( PSILIST( MXXNPSI,MXXACTV ),STAT=IOS )
        CALL CHECKMEM( IOS, 'PSILIST', PROGNAME )

C.........  Initialize unique PSI list
        PSILIST = 0    ! array

C.........  Store unique PSI list
        LPOS = -9
        DO I = 1, NALLPSI

            J = IDXPSIA( I )
            POS = ALLPOSA( J )
            PSI = ALLPSIA( J )

            IF( POS .NE. LPOS ) THEN
                N = 0
                LPSI = -9
                LPOS = POS
            END IF

            IF( PSI .NE. LPSI ) THEN
                N = N + 1
                LPSI = PSI
                PSILIST( N,POS ) = PSI
            END IF

        END DO
        
C.........  Sort emission factor cross-reference entries. Since CPOS was used
C           in building CSRCTA, and CPOS will equal "0" when the x-ref entry is
C           not pollutant-specific, the non-pollutant-specific entries will
C           always appear first.  This is necessary for the table-generating
C           subroutines.

        CALL SORTIC( NXREF, INDXTA, CSRCTA )

C.........  Exit the routine now if no processing requested
        IF ( .NOT. PROCESS ) RETURN

C.........  Generate cross-reference tables and assign x-ref entry to the
C           sources
        CALL XREFTBL( 'EMISFACS', NXREF )

C.........  Deallocate temporary unsorted arrays
        DEALLOCATE( INDXTA, ISPTA, CSCCTA, CSRCTA )
        DEALLOCATE( IDXPSIA, ALLPOSA, ALLPSIA )

C.........  Rewind file
        REWIND( FDEV )

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of temporal' // CRLF() // BLANK5 //
     &         'cross reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   I/O formats........................... 93xxx

93000	FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010	FORMAT( 10( A, :, I7, :, 2X ) )

        END SUBROUTINE RDEFXREF

