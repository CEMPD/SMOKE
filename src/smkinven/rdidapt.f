
        SUBROUTINE RDIDAPT( FDEV, NRAWIN, NRAWBV, WKSET,
     &                      NRAWOUT, EFLAG, NDROP, EDROP )

C***********************************************************************
C  subroutine body starts at line 187
C
C  DESCRIPTION:
C      This subroutine reads the IDA format point-source inventory
C      files.  It call be called multiple times for multiple files.
C
C  PRECONDITIONS REQUIRED:
C      Files must be opened and their unit numbers stored in FDEV() in the
C      order listed in the description.  OTHER?
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      created by M. Houyoux (04/99)
C      allow for activities by A. Holland (02/02) 
C
C*************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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

C...........   MODULES for public variables
C...........   This module is the point source inventory arrays
        USE MODSOURC

C.........  This module contains the lists of unique inventory information
        USE MODLISTS

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
         INCLUDE 'CONST3.EXT'    !  physical and mathematical constants
         INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
      	LOGICAL                CHKINT
        LOGICAL                CHKREAL
        CHARACTER*2            CRLF
        INTEGER                ENVINT
        LOGICAL                ENVYN
        INTEGER                GETFLINE
        INTEGER                GETNLIST
        INTEGER                INDEX1
        INTEGER                STR2INT
        REAL                   STR2REAL
        REAL                   YR2DAY 

        EXTERNAL    CRLF, ENVYN, GETFLINE, GETNLIST, INDEX1, 
     &              STR2INT, STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
C...........   NOTE that NDROP and EDROP are not used at present
        INTEGER     , INTENT (IN) :: FDEV   ! unit number of input file
        INTEGER     , INTENT (IN) :: NRAWIN ! total raw record-count 
        INTEGER     , INTENT (IN) :: NRAWBV ! total raw record times pols
        INTEGER     , INTENT (IN) :: WKSET  ! weekly profile interpretation
        INTEGER     , INTENT(OUT) :: NRAWOUT! outgoing source * pollutants
        LOGICAL     , INTENT(OUT) :: EFLAG  ! outgoing error flag
        INTEGER     ,INTENT(INOUT):: NDROP  !  number of records dropped
        REAL        ,INTENT(INOUT):: EDROP( MXIDAT )  ! emis dropped per pol

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: BLIDLEN  = 6   ! width of boiler field
        INTEGER, PARAMETER :: DESCLEN  = 40  ! width of plant description field
        INTEGER, PARAMETER :: MXVARFIL = 53  ! maximum pollutants in file
        INTEGER, PARAMETER :: PTOTWIDE = 52  ! total width of all pol fields
        INTEGER, PARAMETER :: PTNONPWD = 249 ! width of non-pol fields

C...........   Local parameters, dependent
        INTEGER, PARAMETER :: LINSIZ  = PTNONPWD + MXVARFIL * PTOTWIDE

C...........   Local parameter arrays...
C...........   Start and end positions in the file format of the first set
C              of pollutant fields.
        INTEGER, PARAMETER :: ISINIT( NPTPPOL3 ) = 
     &                              ( / 250,263,276,283,286,296,299 / )

        INTEGER, PARAMETER :: IEINIT( NPTPPOL3 ) = 
     &                              ( / 262,275,282,285,295,298,301 / )

C...........   Local arrays
        INTEGER          IS( NPTPPOL3 )  ! start position for each pol char
        INTEGER          IE( NPTPPOL3 )  ! end position for each pol char
        
C...........   Local allocatable arrays
	CHARACTER*50, ALLOCATABLE  ::  SEGMENT( : ) ! list-formatted strings

C...........   Counters of total number of input records
        INTEGER, SAVE :: NSRCSAV = 0 ! cumulative source count
        INTEGER, SAVE :: NSRCVAR = 0 ! cumulative source x pollutants count

C.........  Temporary variables for storing source characteristics.  These
C           variables must be the width of the fields for global source
C           characteristics definition for use in BLDCSRC.
        CHARACTER(LEN=PLTLEN3) FCID  ! tmp plant ID
        CHARACTER(LEN=CHRLEN3) PTID  ! tmp point ID
        CHARACTER(LEN=CHRLEN3) SKID  ! tmp stack ID
        CHARACTER(LEN=CHRLEN3) SGID  ! tmp segment ID

C...........   Other local variables
        INTEGER         I, J, K, L, N, V  ! counters and indices

        INTEGER         CNY     !  county code
        INTEGER         COD     !  tmp pollutant position in INVDNAM
        INTEGER         ES      !  counter for source x pollutants
        INTEGER         FIP, SCC, SIC  ! tmp fip, scc, sic
        INTEGER         ICC     !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  line counter
        INTEGER, SAVE:: MXWARN  !  maximum number of warnings
        INTEGER, SAVE:: NNOTE =0!  no. of notes in this routine
        INTEGER         NPVAR   !  no. of variables per data
        INTEGER         NSEG    !  no. of input segments
        INTEGER         NVAR    !  number of pollutants in file
        INTEGER, SAVE:: NWARN =0!  number of warnings in this routine
        INTEGER         SS      !  counter for sources
        INTEGER         STA     !  state code
        INTEGER         TPF     !  tmp temporal adjustments setting

        REAL            DAY2YR  ! factor to convert from daily data to annual
        REAL            VAL     !  tmp data value
        REAL            VANN    !  tmp annual data value
        REAL            RBUF    !  tmp real value
        REAL            LAT     !  tmp Y-coordinate
        REAL            LON     !  tmp X-coordinate
        REAL            DM, HT, FL, TK, VL  ! Temporary stack parms

        LOGICAL, SAVE :: CFLAG              ! true: recalc vel w/ flow & diam
        LOGICAL, SAVE :: FFLAG    = .FALSE. ! true: fill in 0. annual with seasonal
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: 1st time routine called
        LOGICAL, SAVE :: WFLAG    = .FALSE. ! true: all lat-lons to western hemi
        LOGICAL          FIXED              ! true: input file is fixed-format
        LOGICAL          MISSFLAG           ! true:  all emissions missing on line


        CHARACTER*300   MESG    !  message buffer

        CHARACTER(LEN=POLLEN3)  CCOD  ! character pollutant index to INVDNAM
        CHARACTER(LEN=FIPLEN3)  CFIP  ! character FIP code
        CHARACTER(LEN=IOVLEN3)  CPOL  ! tmp pollutant code
        CHARACTER(LEN=LINSIZ)   LINE  ! input line from inventory file
        CHARACTER(LEN=SCCLEN3)  TSCC  ! tmp scc

        CHARACTER*16 :: PROGNAME = 'RDIDAPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDIDAPT

C.........  Set up settings the first time the subroutine is called
        IF( FIRSTIME ) THEN

C.............  Get settings from the environment
            CFLAG = ENVYN( 'VELOC_RECALC', 
     &                 'Flag for recalculating velocity', .FALSE., IOS )

            WFLAG = ENVYN( 'WEST_HSPHERE',
     &                 'Western hemisphere flag', .TRUE., IOS )

            MESG = 'Fill in 0. annual data based on seasonal data.'
            FFLAG = ENVYN( 'FILL_ANN_WSEAS', MESG, .FALSE., IOS )

            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

        ENDIF

C.........  Reinitialize for multiple subroutine calls
        EFLAG = .FALSE.
        ICC   = -9
        INY   = 0
        NVAR  = 0
        FIXED = .TRUE.

C........................................................................
C.............  Head of the main read loop  .............................
C........................................................................

        SS   = NSRCSAV
        ES   = NSRCVAR
        IREC = 0
        DO

C.............  Read a line of IDA file as a character string
C.............  If line is a header line or blank, it will advance anyway
            READ( FDEV, 93010, END=199, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading inventory file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            L = LEN_TRIM( LINE )  ! store width of line and check

C.............  Skip blank lines
            IF( L .EQ. 0 ) CYCLE

C.............  Scan for header lines and check to ensure all are set 
C               properly
            CALL GETHDR( MXVARFIL, .TRUE., .TRUE., .TRUE., 
     &                   LINE, ICC, INY, NVAR, IOS )

C.............  Interpret error status
            IF( IOS .EQ. 4 ) THEN
                WRITE( MESG,94010 ) 
     &                 'Maximum allowed data variables ' //
     &                 '(MXVARFIL=', MXVARFIL, CRLF() // BLANK10 //
     &                 ') exceeded in input file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.

            END IF

C.............  If a header line was encountered...
            IF( IOS .GE. 0 ) THEN
            
C................  Check if types is set to determine free or fixed format
	        IF ( LINE(2:5) .EQ. 'TYPE' ) THEN
                
C...................  Try to find activity in name of data, otherwise,
C                     assume emissions
		    MESG = LINE
                    CALL UPCASE( MESG )
                    I = INDEX( MESG, 'ACTIVITY' )
                    IF ( I .GT. 0 ) THEN
                        FIXED = .FALSE.
                    END IF
                END IF
                
C................  Go to next line
		CYCLE
                
            END IF
            
            IF ( FIRSTIME ) DAY2YR = 1. / YR2DAY( INY )
            
C.............   Allocate memory for line segments if not already done
	    IF ( .NOT. ALLOCATED( SEGMENT ) ) THEN
            
C................  Compute the number of segments - different depending on 
C                  whether pollutant or activity
		NSEG = 37
                DO V = 1, NVAR
                    J = DATPOS( V )
                    IF( INVSTAT( J ) .EQ. 1 ) NSEG = NSEG + NPPOL
                    IF( INVSTAT( J ) .EQ. -1 ) NSEG = NSEG + NPACT
                END DO
                
                ALLOCATE( SEGMENT( NSEG ), STAT=IOS )
                CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )
                SEGMENT = ' '   ! array
                
            END IF
            
C.............  Separate lines into parts for fixed format
	    IF ( FIXED ) THEN
            
C................  Initialize fixed-format field positions
                IS = ISINIT - PTOTWIDE   ! array
                IE = IEINIT - PTOTWIDE   ! array
                
                SEGMENT(  1 ) = ADJUSTL( LINE (   1:2   ) )  ! state
                SEGMENT(  2 ) = ADJUSTL( LINE (   3:5   ) )  ! county
                SEGMENT(  3 ) = ADJUSTL( LINE (   6:20  ) )  ! plant ID
                SEGMENT(  4 ) = ADJUSTL( LINE (  21:35  ) )  ! point ID
                SEGMENT(  5 ) = ADJUSTL( LINE (  36:47  ) )  ! stack ID
                SEGMENT(  6 ) = ADJUSTL( LINE (  48:53  ) )  ! DOE plant ID
                SEGMENT(  7 ) = ADJUSTL( LINE (  54:59  ) )  ! boiler ID
                SEGMENT(  8 ) = ADJUSTL( LINE (  60:61  ) )  ! segment ID
                SEGMENT(  9 ) = ADJUSTL( LINE (  62:101 ) )  ! plant descr.
                SEGMENT( 10 ) = ADJUSTL( LINE ( 102:111 ) )  ! SCC code
                
                FCID = SEGMENT( 3 )
                PTID = SEGMENT( 4 )
                SKID = SEGMENT( 5 )
                SGID = SEGMENT( 8 )
                
                K = 37
                DO V = 1, NVAR
                
                    DO I = 1, NPPOL
                        K = K + 1
                        IS( I ) = IS( I ) + PTOTWIDE
                        IE( I ) = IS( I ) + PTOTWIDE
                        
                        SEGMENT( K ) = LINE( IS( I ):IE( I ) )

                    END DO
                    
                END DO
                
C.............  Separate line into parts for list format
	    ELSE
            
                CALL PARSLINE( LINE, NSEG, SEGMENT )
                
            END IF
            
C.............  Check state/county codes, error for missing
	    IF( .NOT. CHKINT( SEGMENT( 1 ) ) .OR.
     &          .NOT. CHKINT( SEGMENT( 2 ) ) .OR.
     &          SEGMENT( 1 ) .EQ. ' '        .OR.
     &          SEGMENT( 2 ) .EQ. ' '             ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: State and/or county ' //
     &                   'code is non-integer or missing at line', IREC
                CALL M3MESG( MESG )
            END IF
            
C..............  Set the default temporal resolution of the data
	    TPF = MTPRFAC * WKSET
            
C..............  Perform unit conversions on input data

C.............  Read stack height and convert units
            CALL READ_REAL( 4, IREC, .TRUE., LINE( 120:123 ), 
     &                      'stack height' , HT, EFLAG )
            IF( HT .LT. 0. ) HT = 0.
            HT   = HT * FT2M                 ! ft to m
            
C.............  Read stack diameter and convert units
            CALL READ_REAL( 6, IREC, .TRUE., LINE( 124:129 ), 
     &                      'stack diameter' , DM, EFLAG )
            IF( DM .LT. 0. ) DM = 0.
            DM   = DM * FT2M                 ! ft to m

C.............  Read exit temperature and convert units
            CALL READ_REAL( 4, IREC, .TRUE., LINE( 130:133 ), 
     &                      'stack exit temperature' , TK, EFLAG )
            IF( TK .LT. 0. ) THEN
                TK = 0.
            ELSE
                TK   = ( TK - 32 ) * FTOC + CTOK ! F to K
            END IF
            
C.............  Read exit flow rate and convert units
            CALL READ_REAL( 10, IREC, .TRUE., LINE( 134:143 ), 
     &                      'stack exit temperature' , FL, EFLAG )
            IF( FL .LT. 0. ) FL = 0.
            FL   = FL * FT2M3                ! ft^3/s to m^3/s

            
C.............  Read exit velocity and convert units
            CALL READ_REAL( 7, IREC, .TRUE., LINE( 144:152 ), 
     &                      'stack exit velocity' , VL, EFLAG )

C.............  Recompute velocity if that option has been set
            IF ( CFLAG .OR. VL .EQ. 0. ) THEN
                RBUF = ( 0.25 * PI * DM * DM )
                IF ( RBUF .GT. 0 ) VL = FL / RBUF

C.............  Otherwise, convert units
            ELSE
                VL = VL * FT2M   ! ft/s to m/s
            END IF
            
C.............  Read Standard Industrial Code
            CALL READ_INTEGER( 4, IREC, .TRUE., LINE( 227:230 ), 
     &                         'SIC', SIC, EFLAG )

C.............  If SIC is missing, set to zero
            SIC = MAX( SIC, 0 )

C.............  Read stack latitude
            CALL READ_REAL( 9, IREC, .FALSE., LINE( 231:239 ), 
     &                      'latitude' , LAT, EFLAG )
            
C.............  Read stack longitude and correct hemisphere if necessary
            CALL READ_REAL( 9, IREC, .FALSE., LINE( 240:248 ), 
     &                      'longtiude' , LON, EFLAG )
            IF( WFLAG .AND. LON .GT. 0 ) LON = -LON

C.............  Make sure that all of the needed real values are real...

C.............  Emissions or activity and associated data
	    K = 37
            DO V = 1, NVAR
            
              J = DATPOS( V )
              IF( INVSTAT( J ) .EQ. 1 ) NPVAR = NPPOL
              IF( INVSTAT( J ) .EQ. -1 ) NPVAR = NPACT
              
              MISSFLAG = .TRUE.
              DO N = 1, NPVAR
              
                  K = K + 1
                  
C.................  Ensure field is a real
                CALL READ_REAL( 50, IREC, .TRUE., SEGMENT( K ),
     &                          'inventory data', VAL, EFLAG )
              
C.................  Convert field to a numeric value and store annual
C                   value for next iteration
	        IF( N .EQ. 1 ) VANN = VAL
              
                IF( SEGMENT( K ) .NE. ' ' ) THEN
              
C....................  If reading emissions data (fixed) and annual data is
C                      zero (VANN=0) but ozone-season data is not missing
C                      (SEGMENT(K)=' ',N=2), fill in annual data with ozone
C                      season data, if user has requested it

	          IF( FIXED .AND. VANN .EQ. 0. .AND.
     &                FFLAG .AND. N    .EQ. 2        ) THEN
     
                      CALL READ_REAL( 50, IREC, .TRUE., SEGMENT( K ),
     &                                'seasonal data', VAL, EFLAG )
                                  
                      VAL = VAL * DAY2YR
                  
                      WRITE( SEGMENT( K-1 ), '(E10.3)' ) VAL
                  
                      NNOTE = NNOTE + 1
                      IF( NNOTE .LE. MXWARN ) THEN
                          WRITE(MESG,94010) 'NOTE: Using ' //
     &                      'seasonal data to fill in annual ' //
     &                      'data'// CRLF()// BLANK10// 'at line',
     &                      IREC, 'for'// TMPNAM( V )
                          CALL M3MESG( MESG )
                      END IF
                  
                  
C...........................  Remove monthly factors for this source.
C                             Note that this will impact ALL pollutants,
C                             even if only one pollutant gets filled.
C                             This is necessary unless TPF is changed to
C                             be pollutant-dependent.

		      TPF = WKSET
                  
                  END IF
              
                  MISSFLAG = .FALSE.     ! data not missing for this V
              
                END IF
          
            END DO
        
              IF( NWARN .LT. MXWARN .AND. MISSFLAG ) THEN
        
                L = LEN_TRIM( TMPNAM( V ) )
                WRITE( MESG, 94010 ) 'WARNING: All columns of ' //
     &                   'data for '// TMPNAM( V )( 1:L ) //
     &                   'are missing at line', IREC
                CALL M3MESG( MESG )
                NWARN = NWARN + 1
                SEGMENT( K ) = '0.'
            
              END IF
        
            END DO
            
            
C.............  If there has been an error, do not try to store any of the
C               records.  Instead go to next line of file.
            IF( EFLAG ) CYCLE
            
            
C.............  Now use the file format definition to parse the line into
C               the various data fields...

	    CALL READ_INTEGER( 50, IREC, .FALSE., SEGMENT( 1 ),
     &                         'state code', STA, EFLAG )
            
            CALL READ_INTEGER( 50, IREC, .FALSE., SEGMENT( 2 ),
     &                         'county code', CNY, EFLAG )
     
     	    FIP = ICC * 100000 + 1000 * STA + CNY
            
            TSCC = SEGMENT( 10 )     ! SCC code
            
C.............  Make adjustments to pad with zeros, if needed
	    WRITE( CFIP,94120 ) FIP
            CALL PADZERO( CFIP )
            CALL PADZERO( TSCC )

C.............  Store source characteristics if dimension is okay
            SS = SS + 1
            
            IF( SS .LE. NRAWIN ) THEN
                IFIPA   ( SS ) = FIP
                ISICA   ( SS ) = SIC
                TPFLGA  ( SS ) = TPF
                INVYRA  ( SS ) = INY
                STKHTA  ( SS ) = HT
                STKDMA  ( SS ) = DM
                STKTKA  ( SS ) = TK
                STKVEA  ( SS ) = VL
                XLOCAA  ( SS ) = LON
                YLOCAA  ( SS ) = LAT
                CSCCA   ( SS ) = TSCC
                CORISA  ( SS ) = SEGMENT( 6 )
                CBLRIDA ( SS ) = SEGMENT( 7 )
                CPDESCA ( SS ) = SEGMENT( 9 )
            END IF
            
            K = 37
            DO V = 1, NVAR
            
                ES = ES + 1
                
                J = DATPOS( V )
                IF( INVSTAT( J ) .EQ. 1 ) NPVAR = NPPOL
                IF( INVSTAT( J ) .EQ. -1 ) NPVAR = NPACT
                
                IF( ES .LE. NRAWBV ) THEN
                
                    INDEXA ( ES ) = ES
                    INRECA ( ES ) = SS
                    
                    WRITE( CCOD,94125 ) J
                    
                    CALL BLDCSRC( CFIP, FCID, PTID, SKID, SGID,
     &                            TSCC, CHRBLNK3, CCOD, CSOURCA( ES ) )
     
C.............  Store main data value
C.............  Units conversion only applies for activities, where
C               NPVAR will be = 1 (no seasonal activity data)
	            K = K + 1
                    CALL READ_REAL( 50, IREC, .TRUE., SEGMENT( K ),
     &                              'inventory data', VAL, EFLAG )
                    POLVLA ( ES, 1 ) = INVDCNV( J ) * VAL
                    
C................  Store related values (if any)
                    DO N = 2, NPVAR
                    
                        K = K + 1
                        CALL READ_REAL( 50 , IREC, .TRUE., SEGMENT( K ),
     &                                  'inventory data', VAL, EFLAG )
                        POLVLA ( ES, N ) = VAL
                        
                    END DO
                    
                END IF
                
            END DO             ! loop through data variables
            
          END DO               ! end of loop for reading input file
          
199	CONTINUE               ! exit from the FDEV-read loop

 	CLOSE( FDEV )
        
        WRITE( MESG,94010 ) 
     &         'IDA FILE processed:'  // CRLF() // BLANK10 //
     &              'This-file source-count', SS - NSRCSAV,
     &         CRLF() // BLANK10 //
     &              'Cumulative source-count', SS,
     &         CRLF() // BLANK10 //
     &              'This-file source*pollutant-count', ES - NSRCVAR,
     &         CRLF() // BLANK10 //
     &              'Cumulative source*pollutant-count', ES

        CALL M3MSG2( MESG )

C.........  Update saved cumulative counts
        NSRCSAV = SS        !  source
        NSRCVAR = ES        !  source*pollutant

C.........  Write message if overflow occurred
        IF( NSRCSAV .GT. NRAWIN ) THEN

            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Source memory allocation ' //
     &             'insufficient for IDA inventory'
            CALL M3MSG2( MESG )

        END IF

        IF( NSRCVAR .GT. NRAWBV ) THEN

            WRITE( MESG, 94010 )
     &	      'INTERNAL ERROR: Number of valid src x variables ' //
     &        'encountered: ', NRAWOUT, CRLF() // BLANK5 //
     &        'Maximum number of raw records allowed: ', NRAWBV
     
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, '', 2 )

        ELSE
            NRAWOUT = NSRCVAR

        END IF

C.........  Deallocate local allocatable arrays 
        DEALLOCATE( TMPNAM, DATPOS, SEGMENT )

C.........  Make sure routine knows it's been called already
        FIRSTIME = .FALSE.

C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A )  

93020   FORMAT( A52 )    ! must match value of PTOTWIDE

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 ) 

        CONTAINS

            SUBROUTINE READ_INTEGER( LENGTH, IREC, OFLAG, STRING, DESC, 
     &                               INTVAL, EFLAG )

            INTEGER      , INTENT ( IN )  :: LENGTH
            INTEGER      , INTENT ( IN )  :: IREC
            LOGICAL      , INTENT ( IN )  :: OFLAG  !  true: field is optional 
            CHARACTER(LEN=LENGTH), INTENT( IN ) :: STRING
            CHARACTER*(*), INTENT ( IN )  :: DESC
            INTEGER      , INTENT ( OUT ) :: INTVAL
            LOGICAL      , INTENT ( OUT ) :: EFLAG 

C.............  Check to see if the field is blank
            IF ( STRING .EQ. ' ' ) THEN

C.................  If field is blank and optional, then set to missing
                IF ( OFLAG ) THEN
                    INTVAL = IMISS3

C.................  If field is blank and not optional, then error
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: required ' // DESC //
     &                     ' is blank at line', IREC
                    CALL M3MESG( MESG )
                    INTVAL = IMISS3
                    RETURN
                END IF
            END IF

C.............  Try to read value and see what the error status is
            READ( STRING, *, IOSTAT = IOS ) INTVAL

C.............  If error, then write message and continue
            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: ' // DESC //
     &                 ' has non-integer value "' // STRING //
     &                 '"' // CRLF() // BLANK 10 // 'at line', IREC
                CALL M3MESG( MESG )
                INTVAL = IMISS3

C.............  Check if missing value
            ELSE IF ( INTVAL .EQ. -9 ) THEN

C.................  If field is missing and optional, then set to missing
                IF ( OFLAG ) THEN
                    INTVAL = IMISS3

C.................  If field is missing and not optional, then error
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: required ' // DESC //
     &                     ' has -9 missing value at line', IREC
                    CALL M3MESG( MESG )
                    INTVAL = IMISS3

                END IF

            END IF

            RETURN

C............................................................................

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE READ_INTEGER

C--------------------------------------------------------------------------
C--------------------------------------------------------------------------

            SUBROUTINE READ_REAL( LENGTH, IREC, OFLAG, STRING, DESC, 
     &                            REALVAL, EFLAG )

            INTEGER      , INTENT ( IN )  :: LENGTH
            INTEGER      , INTENT ( IN )  :: IREC
            LOGICAL      , INTENT ( IN )  :: OFLAG  !  true: field is optional 
            CHARACTER(LEN=LENGTH), INTENT( IN ) :: STRING
            CHARACTER*(*), INTENT ( IN )  :: DESC
            REAL         , INTENT ( OUT ) :: REALVAL
            LOGICAL      , INTENT ( OUT ) :: EFLAG 

C.............  Check to see if the field is blank
            IF ( STRING .EQ. ' ' ) THEN

C.................  If field is blank and optional, then set to missing
                IF ( OFLAG ) THEN
                    REALVAL = BADVAL3

C.................  If field is blank and not optional, then error
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: required ' // DESC //
     &                     ' is blank at line', IREC
                    CALL M3MESG( MESG )
                    REALVAL = BADVAL3
                END IF

                RETURN

            END IF

C.............  Try to read value and see what the error status is
            READ( STRING, *, IOSTAT = IOS ) REALVAL

C.............  If error, then write message and continue
            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: ' // DESC //
     &                 ' has non-readable value "' // STRING //
     &                 '"' // CRLF() // BLANK 10 // 'at line', IREC
                CALL M3MESG( MESG )
                REALVAL = BADVAL3

C.............  Check if missing value
            ELSE IF ( REALVAL .EQ. -9 ) THEN

C.................  If field is missing and optional, then set to zero
                IF ( OFLAG ) THEN
                    REALVAL = BADVAL3

C.................  If field is missing and not optional, then error
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: required ' // DESC //
     &                     ' has -9 missing value at line', IREC
                    CALL M3MESG( MESG )
                    REALVAL = BADVAL3

                END IF

            END IF

            RETURN

C............................................................................

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE READ_REAL


        END SUBROUTINE RDIDAPT
