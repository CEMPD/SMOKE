
        SUBROUTINE RDNTIMB( FDEV, NRAWIN, WKSET, NRAWOUT, EFLAG, 
     &                      NDROP, EDROP )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine reads the NTI format mobile-source inventory
C      files.  It can read multiple NTI files, if needed.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      copied from rdntimb.f by C. Seppanen (11/02) 
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

C.........  This module is for mobile-specific data
        USE MODMOBIL
        
C.........  This module contains the lists of unique inventory information
        USE MODLISTS

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
        INTEGER                FIND1
        INTEGER                FINDC
        INTEGER                INDEX1
        INTEGER                STR2INT
        REAL                   STR2REAL
        REAL                   YR2DAY 

        EXTERNAL    CHKINT, CHKREAL, CRLF, ENVINT, ENVYN, FIND1, FINDC, 
     &              INDEX1, STR2INT, STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
C...........   NOTE that NDROP and EDROP are not used at present
        INTEGER     , INTENT (IN) :: FDEV   ! unit number of input file
        INTEGER     , INTENT (IN) :: NRAWIN ! total raw record-count 
        INTEGER     , INTENT (IN) :: WKSET  ! weekly profile interpretation
        INTEGER     , INTENT(OUT) :: NRAWOUT! outgoing source * pollutants
        LOGICAL     , INTENT(OUT) :: EFLAG  ! outgoing error flag
        INTEGER     ,INTENT(INOUT):: NDROP  !  number of records dropped
        REAL        ,INTENT(INOUT):: EDROP( MXIDAT )  ! emis dropped per pol

C...........   Local parameters
        INTEGER, PARAMETER :: MXDATFIL = 60  ! arbitrary max no. data variables

C...........   Counters of total number of input records
        INTEGER, SAVE :: NSRCSAV = 0 ! cumulative source count

C...........   Other local variables
        INTEGER         I, J    ! counters and indices

        INTEGER         COD     !  pollutant position in INVDNAM
        INTEGER         FIP     !  tmp FIPS code
        INTEGER         ICC     !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  line counter
        INTEGER         IVT     ! tmp vehicle type code
        INTEGER, SAVE:: MXWARN  !  maximum number of warnings
        INTEGER         NCASPOLS!  number of pollutants per CAS
        INTEGER      :: NSEG = 6!  number of segments in line
        INTEGER         NPOA    !  number of pollutants in file
        INTEGER, SAVE:: NWARN =0!  number of warnings in this routine
        INTEGER         SCCLEN  ! length of SCC string 
        INTEGER         SS      !  counter for sources
        INTEGER         RWT     ! roadway type
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

        CHARACTER*300   MESG         !  message buffer
        CHARACTER*20  VIDFMT         ! vehicle type ID format
        CHARACTER*20  RWTFMT         ! roadway type number format
        
        CHARACTER(LEN=50)      SEGMENT( 6 ) ! segments of line
        CHARACTER(LEN=POLLEN3) CCOD  ! character pollutant index to INVDNAM
        CHARACTER(LEN=FIPLEN3) CFIP  ! character FIPS code
        CHARACTER(LEN=VIDLEN3) CIVT  ! tmp vehicle type ID
        CHARACTER(LEN=LNKLEN3) CLNK  ! tmp link ID 
        CHARACTER(LEN=RWTLEN3) CRWT  ! tmp roadway type
        CHARACTER(LEN=IOVLEN3) CPOL  ! tmp pollutant name
        CHARACTER(LEN=300)     LINE  ! input line from inventory file
        CHARACTER(LEN=SCCLEN3) TSCC  ! tmp scc
        CHARACTER(LEN=CASLEN3) TCAS  ! tmp cas number
        CHARACTER(LEN=VTPLEN3) VTYPE ! tmp vehicle type

        CHARACTER*16 :: PROGNAME = 'RDNTIMB' ! Program name

C***********************************************************************
C   begin body of subroutine RDNTIMB

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

C.........  Create formats
        WRITE( VIDFMT, '("(I",I2.2,")")' ) VIDLEN3
        WRITE( RWTFMT, '("(I",I2.2,")")' ) RWTLEN3
        
C.........  Make sure the file is at the beginning
        REWIND( FDEV )

C........................................................................
C.............  Head of the main read loop  .............................
C........................................................................

        SS   = NSRCSAV
        IREC = 0
        CLNK = ' '
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
     &          .NOT. CHKREAL( SEGMENT( 6 ) )      ) THEN
     	
     	        EFLAG = .TRUE.
     	        WRITE( MESG,94010 ) 'ERROR: Emissions data, ' //
     &                 'are not a number or have ' //
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

C.............  Define day to year conversion factor
            DAY2YR  = 1. / YR2DAY( INY )

C.............  Set country/state/county code
            FIP  = ICC  * 100000 +
     &             1000 * STR2INT( SEGMENT( 1 ) ) +
     &                    STR2INT( SEGMENT( 2 ) )

C.............  Save SCC value
            TSCC = SEGMENT( 3 )

C.............  Check that vehicle type is an integer and in list of valid types
            IF( .NOT. CHKINT( TSCC( 3:6 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Vehicle type "' //
     &                 TSCC( 3:6 ) // '" at line', IREC, 'is invalid'
                CALL M3MESG( MESG )
            ELSE
                IVT = STR2INT( TSCC( 3:6 ) )
                DO J = 1, NVTYPE
                    IF( IVT == IVTIDLST( J ) ) EXIT
                END DO
                	
                IF( J > NVTYPE ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'ERROR: Vehicle type "' //
     &                     TSCC( 3:6 ) // '" at line', IREC,
     &                     'was not found in list of valid types'
                    CALL M3MESG( MESG )
                ELSE
                    VTYPE = CVTYPLST( J )
                END IF
            END IF

C.............  Check that road class is an integer and in list of valid types
            IF( .NOT. CHKINT( TSCC( 8:10 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Road class "' // 
     &                 TSCC( 8:10 ) // '" at line', IREC, 'is invalid'
                CALL M3MESG( MESG )
            ELSE
                RWT = STR2INT( TSCC( 8:10 ) )
                J = FIND1( RWT, NRCLAS, AMSRDCLS )
                
                IF( J <= 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Road class "' //
     &                     TSCC( 8:10 ) // '" at line', IREC,
     &                     'was not found in list of valid classes'
                    CALL M3MESG( MESG )
                ELSE
                    RWT = RDWAYTYP( J )  
                END IF
            END IF

C.............  If there has been an error, do not try to store any of the
C               records.  Instead go to next line of file.
            IF( EFLAG ) CYCLE

C.............  Set the default temporal resolution of the data
            TPF  = MTPRFAC * WKSET

C.............  Create string source characteristics and pad with zeroes
            WRITE( CFIP,94120 ) FIP
            CALL PADZERO( CFIP )
            CALL PADZERO( TSCC )
            WRITE( CRWT,RWTFMT ) RWT
            WRITE( CIVT,VIDFMT ) IVT

C.............  Read emissions data as real values
            CALL READ_REAL( 50, IREC, .FALSE., SEGMENT( 5 ),
     &                      'annual emissions', EANN, EFLAG )
            CALL READ_REAL( 50, IREC, .TRUE., SEGMENT( 6 ),
     &                      'average day emissions', EOZN, EFLAG )

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
                EMISBYPOL( I ) = EMISBYPOL( I ) + EANN
                
C.................  Check if current pollutant is kept
                IF( .NOT. ITKEEPA( SCASIDX( I ) ) ) THEN
                    CYCLE
                END IF
                                
C.................  Increment source number
                SS = SS + 1

C.................  Store source information
                IF( SS <= NRAWIN ) THEN
                    IFIPA  ( SS ) = FIP
                    IRCLASA( SS ) = RWT
                    IVTYPEA( SS ) = IVT
                    CLINKA ( SS ) = CLNK
                    CVTYPEA( SS ) = VTYPE
                    TPFLGA ( SS ) = TPF
                    INVYRA ( SS ) = INY
                    CSCCA  ( SS ) = TSCC
                    XLOC1A ( SS ) = BADVAL3
                    YLOC1A ( SS ) = BADVAL3
                    XLOC2A ( SS ) = BADVAL3
                    YLOC2A ( SS ) = BADVAL3
                    POLVLA ( SS,NEM ) = POLANN * INVDCNV( COD )
                    POLVLA ( SS,NOZ ) = POLOZN
                    
                    CALL BLDCSRC( CFIP, CRWT, CLNK, CIVT, TSCC, 
     &                            CHRBLNK3, CHRBLNK3, CCOD, 
     &                            CSOURCA( SS ) )
                 END IF

            END DO

        END DO          !  to head of FDEV-read loop

        CLOSE( FDEV )

        WRITE( MESG,94010 ) 
     &         'NTI FILE processed:'  // CRLF() // BLANK10 //
     &              'This-file source-count', SS - NSRCSAV,
     &         CRLF() // BLANK10 //
     &              'Cumulative source-count', SS

        CALL M3MSG2( MESG )

C.........  Update saved cumulative counts
        NSRCSAV = SS        !  source

C.........  Write message if overflow occurred
        IF( NSRCSAV > NRAWIN ) THEN

            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Source memory allocation ' //
     &             'insufficient for NTI inventory'
            CALL M3MSG2( MESG )

        ELSE
            NRAWOUT = NSRCSAV
        
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


        END SUBROUTINE RDNTIMB
