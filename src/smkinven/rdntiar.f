
        SUBROUTINE RDNTIAR( FDEV, NRAWBP, WKSET, CURREC, EFLAG, 
     &                      NDROP, EDROP )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine reads the NTI format area-source inventory
C      files.  It can read multiple NTI files, if needed.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      copied from rdemsar.f by C. Seppanen (11/02) 
C
C**************************************************************************
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

C...........   MODULES for public variables
C...........   This module is the point source inventory arrays
        USE MODSOURC

C.........  This module contains the lists of unique inventory information
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL                CHKINT
        LOGICAL                CHKREAL
        CHARACTER*2            CRLF
        INTEGER                ENVINT
        LOGICAL                ENVYN
        INTEGER                FINDC
        INTEGER                INDEX1
        INTEGER                STR2INT
        REAL                   STR2REAL
        REAL                   YR2DAY 

        EXTERNAL    CHKINT, CHKREAL, CRLF, ENVINT, ENVYN, FINDC, 
     &              INDEX1, STR2INT, STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
C...........   NOTE that NDROP and EDROP are not used at present
        INTEGER     , INTENT (IN) :: FDEV   ! unit number of input file
        INTEGER     , INTENT (IN) :: NRAWBP ! total raw record times pols
        INTEGER     , INTENT (IN) :: WKSET  ! weekly profile interpretation
        INTEGER     ,INTENT(INOUT):: CURREC ! current no. source * pollutants
        LOGICAL     , INTENT(OUT) :: EFLAG  ! outgoing error flag
        INTEGER     ,INTENT(INOUT):: NDROP  !  number of records dropped
        REAL        ,INTENT(INOUT):: EDROP( MXIDAT )  ! emis dropped per pol

C...........   Local parameters
        INTEGER, PARAMETER :: MXDATFIL = 60  ! arbitrary max no. data variables

C...........   Other local variables
        INTEGER         I       ! counters and indices

        INTEGER         COD     !  pollutant position in INVDNAM
        INTEGER         ES      !  counter for source x pollutants
        INTEGER         FIP     !  tmp FIPS code
        INTEGER         ICC     !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  line counter
        INTEGER         LDEV    !  unit no. for log file
        INTEGER, SAVE:: MXWARN  !  maximum number of warnings
        INTEGER         NCASPOLS!  number of pollutants per CAS
        INTEGER      :: NSEG = 9!  number of segments in line
        INTEGER         NPOA    !  number of pollutants in file
        INTEGER, SAVE:: NWARN =0!  number of warnings in this routine
        INTEGER         NWRLINE !  number of lines of file written to log
        INTEGER         SCCLEN         ! length of SCC string 
        INTEGER         SS      !  counter for sources
        INTEGER         TPF     !  tmp temporal adjustments setting
        INTEGER         SCASPOS !  position of CAS number in sorted array
        INTEGER         UCASPOS !  position of CAS number in unique array

        REAL            CEFF    !  tmp control effectiveness
        REAL            DAY2YR  !  factor to convert from daily data to annual
        REAL            EANN    !  tmp annual-ave emission value
        REAL            EOZN    !  tmp ozone-season-ave emission value
        REAL            REFF    !  tmp rule effectiveness
        REAL            RPEN    !  tmp rule penetration
        REAL            POLANN  !  annual emission for pollutant
        REAL            POLOZN  !  ozone-season emission for pollutant
        REAL            POLFAC  !  factor to convert from CAS data to pollutant

        LOGICAL, SAVE:: FFLAG    = .FALSE. ! true: fill in 0. annual with seasonal
        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called

        CHARACTER*300   MESG    !  message buffer

        CHARACTER(LEN=50)      SEGMENT( 9 ) ! segments of line
        CHARACTER(LEN=POLLEN3) CCOD  ! character pollutant index to INVDNAM
        CHARACTER(LEN=FIPLEN3) CFIP  ! character FIP code
        CHARACTER(LEN=IOVLEN3) CPOL  ! tmp pollutant name
        CHARACTER(LEN=300)     LINE  ! input line from inventory file
        CHARACTER(LEN=SCCLEN3) TSCC  ! tmp scc
        CHARACTER(LEN=CASLEN3) TCAS  ! tmp cas number
        
        CHARACTER(LEN=300)     TENLINES( 10 ) ! first ten lines of inventory file

        CHARACTER*16 :: PROGNAME = 'RDNTIAR' ! Program name

C***********************************************************************
C   begin body of subroutine RDNTIAR

        IF ( FIRSTIME ) THEN

            MESG = 'Fill in 0. annual data based on seasonal data.'
            FFLAG = ENVYN( 'FILL_ANN_WSEAS', MESG, .FALSE., IOS )

            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

        END IF

C.........  Reinitialize for multiple subroutine calls
        EFLAG = .FALSE.
        ICC   = -9
        INY   = 0
        NPOA  = 0
        
C.........  Get log file number for reports
        LDEV = INIT3()        
        NWRLINE = 0

C.........  Make sure the file is at the beginning
        REWIND( FDEV )

C........................................................................
C.............  Head of the main read loop  .............................
C........................................................................

        ES   = CURREC
        IREC = 0
        DO

C.............  Read a line of NTI file as a character string
            READ( FDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

C.............  Check I/O error status
            IF ( IOS > 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading inventory file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            ELSE IF ( IOS < 0 ) THEN  ! reached end of file
                EXIT
            END IF

            IF ( LINE == ' ' ) CYCLE      ! skip if line is blank

C.............  Scan for header lines and check to ensure all are set 
C               properly (country and year required)
            CALL GETHDR( MXDATFIL, .TRUE., .TRUE., .FALSE., 
     &                   LINE, ICC, INY, NPOA, IOS )

C.............  Interpret error status
            IF( IOS == 4 ) THEN
                WRITE( MESG,94010 ) 
     &                 'Maximum allowed data variables ' //
     &                 '(MXDATFIL=', MXDATFIL, CRLF() // BLANK10 //
     &                 ') exceeded in input file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE IF( IOS > 0 ) THEN
                EFLAG = .TRUE.

            END IF

C.............  If a header line was encountered, go to next line
            IF( IOS >= 0 ) CYCLE

C.............  Write first ten lines to log file
            IF( NWRLINE < 10 ) THEN
            	NWRLINE = NWRLINE + 1
            	TENLINES( NWRLINE ) = TRIM( LINE )
            
                IF( NWRLINE == 10 ) THEN
                    MESG = 'First 10 lines of NTI area inventory:'
                    WRITE( LDEV,* ) TRIM( MESG )
             
                    DO I = 1,NWRLINE
                        WRITE( LDEV,* ) TRIM( TENLINES( I ) )
                    END DO
                END IF
            END IF
            
C.............  Separate line into segments
            CALL PARSLINE( LINE, NSEG, SEGMENT )

C.............  Make sure that all of the needed integer values are integers...

C.............  Check state/county codes, error for missing
            IF( .NOT. CHKINT( SEGMENT( 1 ) ) .OR. 
     &          .NOT. CHKINT( SEGMENT( 2 ) ) .OR.
     &          SEGMENT( 1 ) .EQ. ' '        .OR.
     &          SEGMENT( 2 ) .EQ. ' '             ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: State and/or county ' //
     &                 'code is non-integer or missing at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Determine if SCC has proper length
            SCCLEN = LEN_TRIM( SEGMENT( 3 ) )

            IF ( SCCLEN /= 10 ) THEN   ! check if SCC is 10 char. wide
                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &	                   'SCC not 10 characters wide on line ', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Make sure that all of the needed real values are real...

C.............  Emissions and associated data
            IF( .NOT. CHKREAL( SEGMENT( 5 ) ) .OR.
     &          .NOT. CHKREAL( SEGMENT( 6 ) ) .OR.
     &          .NOT. CHKREAL( SEGMENT( 7 ) ) .OR.
     &          .NOT. CHKREAL( SEGMENT( 8 ) ) .OR.
     &          .NOT. CHKREAL( SEGMENT( 9 ) )      ) THEN
     	
     	        EFLAG = .TRUE.
     	        WRITE( MESG,94010 ) 'ERROR: Emissions data, ' //
     &                 'control percentage, and/or rule ' //
     &                 'percentages are not a number or have ' //
     &                 'bad formatting at line', IREC
                CALL M3MESG( MESG )
     	    END IF
            
            IF( SEGMENT( 5 ) == ' ' .AND.
     &          SEGMENT( 6 ) == ' '       ) THEN
     	
                WRITE( MESG,94010 ) 'WARNING: All emissions ' //
     &                 'data are missing at line', IREC
                CALL M3MESG( MESG )
                SEGMENT( 5 ) = '0.'
                SEGMENT( 6 ) = '0.'
            END IF

C.............  Find CAS number in list of unique CAS's from INVTABLE
            TCAS = SEGMENT( 4 )
            UCASPOS = FINDC( TCAS, NUNIQCAS, UNIQCAS )
            IF( UCASPOS < 1 ) THEN
                WRITE( MESG,94010 ) 'Source dropped: ' //
     &                 'CAS number "' // TRIM( TCAS ) //
     &                 '" in emission file at line', IREC, CRLF() //
     &                 BLANK5 // 'is not in inventory pollutants list'
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Get total number of pollutants per CAS number and
C               position of CAS number in sorted array
            NCASPOLS = UCASNPOL( UCASPOS )
            SCASPOS  = UCASIDX( UCASPOS )
            
C.............  If there has been an error, do not try to store any of the
C               records.  Instead go to next line of file.
            IF( EFLAG ) CYCLE

C.............  Define day to year conversion factor
            DAY2YR  = 1. / YR2DAY( INY )

C.............  Set country/state/county code
            FIP  = ICC  * 100000 +
     &             1000 * STR2INT( SEGMENT( 1 ) ) +
     &                    STR2INT( SEGMENT( 2 ) )
            WRITE( CFIP,94120 ) FIP
            CALL PADZERO( CFIP )

C.............  Save SCC value
            TSCC = SEGMENT( 3 )

C.............  Set the default temporal resolution of the data
            TPF  = MTPRFAC * WKSET

C.............  Read emissions data as real values
            CALL READ_REAL( 50, IREC, .FALSE., SEGMENT( 5 ),
     &                      'annual emissions', EANN, EFLAG )
            CALL READ_REAL( 50, IREC, .TRUE., SEGMENT( 6 ),
     &                      'average day emissions', EOZN, EFLAG )
            CALL READ_REAL( 50, IREC, .TRUE., SEGMENT( 7 ),
     &                      'control efficiency', CEFF, EFLAG )
            CALL READ_REAL( 50, IREC, .TRUE., SEGMENT( 8 ),
     &                      'rule effectiveness', REFF, EFLAG )
            CALL READ_REAL( 50, IREC, .TRUE., SEGMENT( 9 ),
     &                      'rule penetration', RPEN, EFLAG )

            IF( CEFF .LT. AMISS3 ) THEN
                WRITE( MESG,94010 ) 'WARNING: Missing control ' //
     &                 'efficiency at line', IREC 
                IF( NWARN .LT. MXWARN ) CALL M3MESG( MESG )
                NWARN = NWARN + 1
            END IF
            
            IF( REFF .LT. AMISS3 ) THEN
                WRITE( MESG,94010 ) 'WARNING: Missing rule ' //
     &                 'effectiveness at line', IREC
                IF( NWARN .LT. MXWARN ) CALL M3MESG( MESG )
                NWARN = NWARN + 1
            END IF

            IF( RPEN .LT. AMISS3 ) THEN
                WRITE( MESG,94010 ) 'WARNING: Missing rule ' //
     &                 'pentration at line', IREC 
                IF( NWARN .LT. MXWARN ) CALL M3MESG( MESG )
                NWARN = NWARN + 1
            END IF

C.............  Replace annual data with ozone-season information if
C               user option is set
            IF( FFLAG      .AND. 
     &          EANN <= 0. .AND.
     &          EOZN >  0.       ) THEN

                WRITE( MESG,94010 ) 'NOTE: Using seasonal ' //
     &                 'emissions to fill in annual emissions' //
     &                 CRLF() // BLANK10 // 'at line', IREC
                CALL M3MESG( MESG )

                EANN = EOZN * DAY2YR
                
C.................  Remove monthly factors for this source
                TPF = WKSET
                
            END IF

C.............  Store emissions by CAS number for reporting
            EMISBYCAS( UCASPOS ) = EMISBYCAS( UCASPOS ) + EANN
            RECSBYCAS( UCASPOS ) = RECSBYCAS( UCASPOS ) + 1
                    
C.............  Loop through pollutants for this CAS number            
            DO I = SCASPOS, SCASPOS + NCASPOLS - 1

C.................  Set pollutant name
                CPOL = ITNAMA( SCASIDX( I ) )
                CALL UPCASE( CPOL )

C.................  Find location of pollutant in master list    
                COD = INDEX1( CPOL, MXIDAT, INVDNAM )
                WRITE( CCOD,94125 ) COD
                
C.................  Apply factor for this pollutant (only if valid)
                POLFAC = ITFACA( SCASIDX( I ) )
                
                IF( EANN > AMISS3 ) THEN
                    POLANN = EANN * POLFAC
                ELSE
                    POLANN = EANN
                END IF
                
                IF( EOZN > AMISS3 ) THEN
                    POLOZN = EOZN * POLFAC
                ELSE
                    POLOZN = EOZN
                END IF
                
C.................  Store emissions by pollutant for reporting
                EMISBYPOL( I ) = EMISBYPOL( I ) + POLANN

C.................  Check if current pollutant is kept
                IF( .NOT. ITKEEPA( SCASIDX( I ) ) ) THEN
                    CYCLE
                END IF
                
C.................  Increment source number
                ES = ES + 1

C.................  Store source information
                IF( ES <= NRAWBP ) THEN
                    IFIPA ( ES     ) = FIP
                    TPFLGA( ES     ) = TPF
                    INVYRA( ES     ) = INY
                    CSCCA ( ES     ) = TSCC
                    POLVLA( ES,NEM ) = POLANN * INVDCNV( COD )
                    POLVLA( ES,NOZ ) = POLOZN
                    POLVLA( ES,NCE ) = CEFF
                    POLVLA( ES,NRE ) = REFF
                    POLVLA( ES,NRP ) = RPEN
                    
                    CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3,
     &                            CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                            CCOD, CSOURCA( ES ) )
     
                    CSOURCA( ES )( ALLLEN3+1:ALLCAS3 ) = ADJUSTR( TCAS )
                    
                 END IF

            END DO

        END DO          !  to head of FDEV-read loop

        CLOSE( FDEV )

        WRITE( MESG,94010 ) 
     &         'NTI FILE processed:'  // CRLF() // BLANK10 //
     &              'This-file source-count', ES - CURREC,
     &         CRLF() // BLANK10 //
     &              'Cumulative source-count', ES

        CALL M3MSG2( MESG )

C.........  Write message if overflow occurred
        IF( ES > NRAWBP ) THEN

            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Source memory allocation ' //
     &             'insufficient for NTI inventory'
            CALL M3MSG2( MESG )

        ELSE
            CURREC = ES
        
        END IF

C.........  Make sure routine knows it's been called already
        FIRSTIME = .FALSE.

C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        CONTAINS

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
                    RETURN
                END IF
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


        END SUBROUTINE RDNTIAR
