
        SUBROUTINE RDINVDATA( FDEV, FNAME, NRAWBP, TFLAG )

C***********************************************************************
C  subroutine body starts at line 133
C
C  DESCRIPTION:
C      This subroutine controls reading an ASCII inventory file for any source 
C      category from one of many formats.  It determines the format and 
C      calls the appropriate reader subroutines. It controls the looping 
C      through multiple files when a list-formatted file is used as input.
C      This routine only reads the data (emissions and activities) from the
C      inventories.
C
C  PRECONDITIONS REQUIRED:
C      Input file unit FDEV opened
C      Inventory pollutant list created: MXIDAT, INVDCOD, and INVDNAM
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: 
C      Functions: 
C
C  REVISION  HISTORY:
C      Created 1/03 by C. Seppanen (based on rdinven.f)
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
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: INRECA, POLVLA, TPFLGA, INDEXA,
     &                      NSTRECS, SRCSBYREC, RECIDX, IPOSCODA, 
     &                      SRCIDA, INVYR, ICASCODA
        
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NEM, NOZ, NEF, NCE, NRE, NRP, 
     &                     NPPOL, NSRC 
        
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: FILFMT, LSTSTR, MXIDAT, INVDCNV, INVDNAM,
     &                      NUNIQCAS, UNIQCAS, UCASNKEP, ITNAMA, 
     &                      SCASIDX, UCASIDX, UCASNPOL, ITKEEPA, ITFACA

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         CHKREAL
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        LOGICAL         ENVYN
        INTEGER         FINDC
        INTEGER         GETINVYR
        INTEGER         INDEX1
        REAL            STR2REAL
        REAL            YR2DAY

        EXTERNAL        CHKREAL, CRLF, ENVINT, ENVYN, FINDC, 
     &                  GETINVYR, INDEX1, STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
        INTEGER,          INTENT (IN) :: FDEV         ! unit no. of inv file
        CHARACTER(LEN=*), INTENT (IN) :: FNAME        ! logical name of file
        INTEGER,          INTENT (IN) :: NRAWBP       ! no. sources with pollutants
        LOGICAL,          INTENT(OUT) :: TFLAG        ! true: PTREF output

C...........   Local parameters
        INTEGER, PARAMETER :: DATALEN3 = 25  ! length of data field
        
C...........   Dropped emissions
        INTEGER         NDROP             !  number of records dropped
        REAL            EDROP  ( MXIDAT ) !  total dropped for each pol/activity

C...........   File units and logical/physical names
        INTEGER         EDEV( 5 )   !  up to 5 EMS-95 emissions files
        INTEGER         TDEV        !  file listed in list formatted input file

C...........   Output from individual reader routines
        CHARACTER(LEN=DATALEN3), ALLOCATABLE :: READDATA( :,: )  ! data values
        CHARACTER(LEN=IOVLEN3),  ALLOCATABLE :: READPOL ( : )    ! pollutant names

C...........   Other local variables
        INTEGER         I, J, SP    !  counters and indices

        INTEGER         CURFIL      !  current file from list formatted inventory
        INTEGER         CURFMT      !  format of current inventory file
        INTEGER         CURSRC      !  current source number
        INTEGER      :: INY = 0     !  tmp inventory year
        INTEGER         IOS         !  i/o status
        INTEGER         INVYEAR     !  inventory year
        INTEGER         IREC        !  no. of records read
        INTEGER         ISTREC      !  no. of records stored
        INTEGER         MXWARN      !  maximum number of warnings
        INTEGER         NPOLPERCAS  !  no. of pollutants per CAS number
        INTEGER         NPOLPERLN   !  no. of pollutants per line of inventory file
        INTEGER      :: NWARN = 0   !  current number of warnings
        INTEGER         POLCOD      !  pollutant code
        INTEGER         TPF         !  temporal adjustments setting
        INTEGER         SCASPOS     !  position of CAS number in sorted array
        INTEGER         UCASPOS     !  position of CAS number in unique array
        INTEGER         WKSET       !  setting for wkly profile TPFLAG component

        REAL            CEFF        !  tmp control effectiveness
        REAL            DAY2YR      !  factor to convert from daily data to annual
        REAL            EANN        !  annual-ave emission value
        REAL            EMFC        !  emission factor
        REAL            EOZN        !  ozone-season-ave emission value
        REAL            REFF        !  rule effectiveness
        REAL            RPEN        !  rule penetration
        REAL            POLFAC      !  factor for current pollutant

        LOGICAL      :: EFLAG  = .FALSE. ! true: error occured
        LOGICAL      :: DFLAG  = .FALSE. ! true: weekday (not full week) nrmlizr 
        LOGICAL      :: FFLAG  = .FALSE. ! true: fill annual data with seasonal
        LOGICAL      :: HDRFLAG          ! true: current line is part of header
        LOGICAL      :: LSTFLG = .FALSE. ! true: using list-fmt inventory file
        LOGICAL      :: LSTTIME = .FALSE. ! true: last time through 

        CHARACTER(LEN=IOVLEN3) POLNAM    !  tmp pollutant name
        CHARACTER(LEN=300)    INFILE     !  input file line buffer
        CHARACTER(LEN=3000)   LINE       !  input file line buffer
        CHARACTER(LEN=300)    MESG       !  message buffer

        CHARACTER*16 :: PROGNAME =  'RDINVDATA' ! program name

C***********************************************************************
C   begin body of subroutine RDINVDATA

C.........  Check if inventory file is list format
        IF( SIZE( FILFMT ) > 1 ) LSTFLG = .TRUE.
        
C.........   Initialize variables for keeping track of dropped emissions
        NDROP = 0
        EDROP = 0.  ! array

C.........  Get setting for interpreting weekly temporal profiles from the
C           environment. Default is false for non-EMS-95 and true for EMS-95
C           inventory inputs.
        DFLAG = .FALSE.
        
        DO I = 1, SIZE( FILFMT )
            IF ( FILFMT( I ) .EQ. EMSFMT ) DFLAG = .TRUE.
        END DO
        
        MESG = 'Use weekdays only to normalize weekly profiles'
        DFLAG = ENVYN( 'WKDAY_NORMALIZE', MESG, DFLAG, IOS )

C.........  Set weekly profile interpretation flag...
C.........  Weekday normalized
        IF( DFLAG ) THEN
            WKSET = WDTPFAC
            MESG = 'NOTE: Setting inventory to use weekday '//
     &             'normalizer for weekly profiles'

C.........  Full-week normalized
        ELSE
            WKSET = WTPRFAC
            MESG = 'NOTE: Setting inventory to use full-week '//
     &             'normalizer for weekly profiles'

        END IF

C.........  Write message
        CALL M3MSG2( MESG )

C.........  If EMS-95 format, check the setting for the interpretation of
C           the weekly profiles
        DO I = 1, SIZE( FILFMT )
            IF( FILFMT( I ) == EMSFMT .AND. 
     &          WKSET /= WDTPFAC ) THEN

                MESG = 'WARNING: EMS-95 format files will be using ' //
     &             'non-standard approach of ' // CRLF() // BLANK10 //
     &             'full-week normalized weekly profiles.  Can ' //
     &             'correct by setting ' // CRLF() // BLANK10 //
     &             'WKDAY_NORMALIZE to Y and rerunning.'
                CALL M3MSG2( MESG )

            ELSE IF( FILFMT( I ) == EPSFMT .AND. 
     &               WKSET /= WTPRFAC ) THEN

                MESG = 'WARNING: EPS2.0 format files will be using ' //
     &             'non-standard approach of ' // CRLF() // BLANK10 //
     &             'weekday normalized weekly profiles.  Can ' //
     &             'correct by setting ' // CRLF() // BLANK10 //
     &             'WKDAY_NORMALIZE to N and rerunning.'
                CALL M3MSG2( MESG )

            END IF
        END DO

C.........  Get annual data setting from environment
        MESG = 'Fill in 0. annual data based on seasonal data.'
        FFLAG = ENVYN( 'FILL_ANN_WSEAS', MESG, .FALSE., IOS )

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

C.........  Set default inventory characteristics (declared in MODINFO) used
C           by the IDA and EPS formats, including NPPOL
        CALL INITINFO( FILFMT( 1 ) )
        
C.........  Allocate memory for storing inventory data
        ALLOCATE( INDEXA( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )
        ALLOCATE( POLVLA( NRAWBP,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLVLA', PROGNAME )
        ALLOCATE( INRECA( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INRECA', PROGNAME )
        ALLOCATE( IPOSCODA( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPOSCODA', PROGNAME )
        ALLOCATE( ICASCODA( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICASCODA', PROGNAME )
        ALLOCATE( TPFLGA( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TPFLGA', PROGNAME )
        ALLOCATE( INVYR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVYR', PROGNAME )

C.........  Initialize pollutant-specific values as missing
        POLVLA  = BADVAL3  ! array

C.........  If inventory is list format, open first file for reading
        CURFIL = 1

        IF( LSTFLG ) THEN
            LINE = LSTSTR( CURFIL )

C.............  Check if line is year packet
            INVYEAR = GETINVYR( LINE )
            
            IF( INVYEAR > 0 ) THEN
                CURFIL = CURFIL + 1
            END IF

C.............  Store path of file name            
            INFILE = LSTSTR( CURFIL )
            
C.............  Open current file
            OPEN( FDEV, FILE=INFILE, STATUS='OLD', IOSTAT=IOS )

C.............  Check for errors while opening file
            IF( IOS /= 0 ) THEN
            
                WRITE( MESG,94010 ) 'Problem at line ', CURFIL, 'of ' //
     &             TRIM( FNAME ) // '.' // ' Could not open file:' //
     &             CRLF() // BLANK5 // TRIM( INFILE ) 
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE
                WRITE( MESG,94010 ) 'Successful OPEN for ' //
     &             'inventory file:' // CRLF() // BLANK5 //
     &             TRIM( INFILE )
                CALL M3MSG2( MESG ) 

            END IF

C.............  Set default inventory characteristics that depend on file format
            CALL INITINFO( FILFMT( CURFIL ) )
        
        END IF

C.........  Allocate memory to store emissions and pollutant from a single line
C.........  For now, set number of pollutants per line to 1 
C           (will change if format is IDA)
        ALLOCATE( READDATA( 1,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'READDATA', PROGNAME )
        ALLOCATE( READPOL( 1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'READPOL', PROGNAME )
        
        CURFMT = FILFMT( CURFIL )

        IREC = 0    ! current record number
        ISTREC = 0  ! current stored record
        SP = 0      ! current source with pollutant index

C.........  Loop through inventory files and read data
        DO

C.............  If reached end of SRCSBYREC array, make sure we finished the file
            IF( ISTREC == NSTRECS ) THEN
!                READ( FDEV, 93000, IOSTAT=IOS ) LINE
!                IF( IOS == 0 ) THEN   ! successful read -> not the end of the file
!                    WRITE( MESG,94010 ) 'INTERNAL ERROR:' //
!     &                  'reached end of records before end ' //
!     &                  'of file.'
!                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
!                ELSE
                    EXIT
!                END IF
            END IF
        
            READ( FDEV, 93000, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS > 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'I/O error', IOS,
     &             'reading inventory file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Check if we've reached the end of the file            
            IF( IOS < 0 ) THEN

C.................  If list format, try to open next file
                IF( LSTFLG ) THEN

C.....................  Close current file and reset counter
                    CLOSE( FDEV )
                    IREC = 0
                
                    CURFIL = CURFIL + 1

C.....................  Check if there are more files to read
                    IF( CURFIL <= SIZE( FILFMT ) ) THEN 
                        INFILE = LSTSTR( CURFIL )
                
                        OPEN( FDEV, FILE=INFILE, STATUS='OLD', 
     &                        IOSTAT=IOS )
                
C.........................  Check for errors while opening file
		     	        IF( IOS /= 0 ) THEN
					
				            WRITE( MESG,94010 ) 'Problem at line ', 
     &  		               CURFIL, 'of ' // TRIM( FNAME ) // 
     &  		               '.' // ' Could not open file:' //
     &		                   CRLF() // BLANK5 // TRIM( INFILE ) 
				            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
			
			            ELSE
				            WRITE( MESG,94010 ) 
     &  		              'Successful OPEN for ' //
     &		                  'inventory file(s):' // CRLF() // 
     &                        BLANK5 // TRIM( INFILE )
				            CALL M3MSG2( MESG ) 
			
				        END IF

C.........................  Set default inventory characteristics that depend on file format
			            CALL INITINFO( FILFMT( CURFIL ) )
			            CURFMT = FILFMT( CURFIL )

C.........................  Reallocate memory to store emissions from a single line
                        DEALLOCATE( READDATA, READPOL )
                        ALLOCATE( READDATA( 1,NPPOL ), STAT=IOS )
                        CALL CHECKMEM( IOS, 'READDATA', PROGNAME )
                        ALLOCATE( READPOL( 1 ), STAT=IOS )
                        CALL CHECKMEM( IOS, 'READPOL', PROGNAME )
                        
C.........................  Skip back to the beginning of the loop
                        CYCLE
			  
C.....................  Otherwise, no more files to read, so exit
			        ELSE
			            LSTTIME = .TRUE.
			            EXIT
			        END IF

C.................  Otherwise, not a list file, so exit
                ELSE
                    LSTTIME = .TRUE.
                    EXIT
                END IF
             
            END IF   ! end check for end of file
            
C.............  Skip blank lines
            IF( LINE == ' ' ) CYCLE

C.............  Process line depending on file format and source category
            SELECT CASE( CURFMT )
            CASE( IDAFMT )
                SELECT CASE( CATEGORY )
                CASE( 'AREA' )
                    CALL RDDATAIDAAR( LINE, READDATA, READPOL, 
     &                                NPOLPERLN, INVYEAR, HDRFLAG, 
     &                                EFLAG )
                END SELECT
            CASE( NTIFMT )
                SELECT CASE( CATEGORY )
                CASE( 'AREA' )
                    CALL RDDATANTIAR( LINE, READDATA, READPOL, 
     &                                INVYEAR, HDRFLAG, EFLAG )
                    NPOLPERLN = 1
                END SELECT
            END SELECT
            
C.............  Check for header lines
            IF( HDRFLAG ) THEN 

C.................  If IDA format, reallocate emissions memory with proper number
C                   of pollutants per line
                IF( CURFMT == IDAFMT .AND. NPOLPERLN /= 0 ) THEN
                    DEALLOCATE( READDATA, READPOL )
                    ALLOCATE( READDATA( NPOLPERLN,NPPOL ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'READDATA', PROGNAME )
                    ALLOCATE( READPOL( NPOLPERLN ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'READPOL', PROGNAME )
                END IF

C.................  Calculate day to year conversion factor
                IF( INVYEAR /= 0 ) THEN
                    DAY2YR = 1. / YR2DAY( INVYEAR )
                END IF
                
                CYCLE
            END IF

C.............  Check that data values are numbers
            DO I = 1, NPOLPERLN
                POLNAM = READPOL( I )
            
                DO J = 1, NPPOL 
                    IF( .NOT. CHKREAL( READDATA( I,J ) ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Emission data, ' //
     &                     'control percentages, and/or emission ' //
     &                     CRLF() // BLANK10 // 'factor for ' //
     &                     TRIM( POLNAM ) // ' are not a number ' //
     &                     'or have bad formatting at line', IREC
                        CALL M3MESG( MESG )
                        EXIT
                    END IF
                END DO  ! end loop over data values
                
                IF( READDATA( I,1 ) == ' ' .AND.
     &              READDATA( I,2 ) == ' '       ) THEN
                    IF( NWARN < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: All emissions ' //
     &                     'data for ' // TRIM( POLNAM ) //
     &                     ' are missing at line', IREC
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                    READDATA( I,1 ) = '0.'
                    READDATA( I,2 ) = '0.'
                END IF            
            END DO  ! end loop over pollutants per line

C.............  Skip rest of loop if an error has occured
            IF( EFLAG ) CYCLE

C.............  Make sure some emissions are kept for this source
            IF( CURFMT == NTIFMT ) THEN
                POLNAM = READPOL( 1 )
                UCASPOS = FINDC( POLNAM, NUNIQCAS, UNIQCAS )
                IF( UCASPOS < 1 ) THEN
!                    WRITE( MESG,94010 ) 'No valid pollutants found ' //
!     &                 'at line', IREC, '. The source will be ' //
!     &                 'dropped.'
!                    CALL M3MESG( MESG )
                    CYCLE
                ELSEIF( UCASNKEP( UCASPOS ) == 0 ) THEN
                    CYCLE
                END IF

            ELSE
C.................  For non-toxic sources, set UCASPOS to use in ICASCODA
                UCASPOS = 0
            END IF

C.............  Increment number of stored records and double check that we are
C               where we're supposed to be
            ISTREC = ISTREC + 1
            IF( SRCSBYREC( RECIDX( ISTREC ),1 ) /= CURFIL .OR.
     &          SRCSBYREC( RECIDX( ISTREC ),2 ) /= IREC        ) THEN
                MESG = 'INTERNAL ERROR: Current record does ' //
     &                 'not match expected record'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF      

C.............  Store inventory year in sorted order
            CURSRC = SRCIDA( SRCSBYREC( RECIDX( ISTREC ),3 ) )
            INVYR( CURSRC ) = INVYEAR

C.............  Convert data to real numbers and check for missing values
            DO I = 1, NPOLPERLN
                POLNAM = READPOL( I )
                
                EANN = STR2REAL( READDATA( I,1 ) )
                IF( EANN == -9 ) EANN = IMISS3
                
                IF( EANN == IMISS3 ) THEN
                    IF( NWARN < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                     'annual emissions for ' // 
     &                     TRIM( POLNAM ) // ' at line', IREC
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                END IF
                
                EOZN = STR2REAL( READDATA( I,2 ) )
                IF( EOZN == -9 ) EOZN = IMISS3
                
                IF( EOZN == IMISS3 ) THEN
                    IF( NWARN < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                     'seasonal emissions for ' //
     &                     TRIM( POLNAM ) // ' at line', IREC
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                END IF

C.................  For area and point, convert control percentages
                IF( CATEGORY == 'AREA' .OR. CATEGORY == 'POINT' ) THEN
                    EMFC = STR2REAL( READDATA( I,3 ) )
                    IF( EMFC == -9 ) EMFC = IMISS3
                
                    CEFF = STR2REAL( READDATA( I,4 ) )
                    IF( CEFF == -9 ) CEFF = IMISS3
                    
                    IF( CEFF == IMISS3 ) THEN
                        IF( NWARN < MXWARN ) THEN
                            WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                         'control efficiency for ' //
     &                         TRIM( POLNAM ) // ' at line', IREC
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF
                    END IF
                    
                    REFF = STR2REAL( READDATA( I,5 ) )
                    IF( REFF == -9 ) REFF = IMISS3
                    
                    IF( REFF == IMISS3 ) THEN
                        IF( NWARN < MXWARN ) THEN
                            WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                         'rule effectiveness for ' //
     &                         TRIM( POLNAM ) // ' at line', IREC
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF
                    END IF
                    
                    RPEN = STR2REAL( READDATA( I,6 ) )
                    IF( RPEN == -9 ) RPEN = IMISS3
                    
                    IF( RPEN == IMISS3 ) THEN
                        IF( NWARN < MXWARN ) THEN
                            WRITE( MESG,94010 ) 'WARNING: Missing ' //
     &                         'rule penetration for ' //
     &                         TRIM( POLNAM ) // ' at line', IREC
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF
                    END IF
                END IF
                
                IF( CATEGORY == 'POINT' ) THEN
                
                END IF

C.................  Set temporal adjustment for this source and pollutant
                TPF = MTPRFAC * WKSET

C.................  Replace annual data with ozone-season information
                IF( FFLAG .AND. EANN >= 0. .AND. EOZN < 0. ) THEN
                    WRITE( MESG,94010 ) 'NOTE: Using seasonal ' //
     &                 'emissions to fill in annual emissions' //
     &                 CRLF() // BLANK10 // 'for ' // TRIM( POLNAM ) //
     &                 ' at line', IREC
                    CALL M3MESG( MESG )
                    
                    EANN = EOZN * DAY2YR

C.....................  Remove monthly factors for this source and pollutant
                    TPF = WKSET
                END IF

C.................  If current format is NTI, check if current CAS number
C                   is split
                IF( CURFMT == NTIFMT ) THEN
                    NPOLPERCAS = UCASNPOL( UCASPOS )
                ELSE
                    NPOLPERCAS = 1
                    POLFAC = 1.
                END IF
                    
                DO J = 0, NPOLPERCAS - 1

C.....................  If NTI format, set current pollutant
                    IF( CURFMT == NTIFMT ) THEN

C.........................  Make sure current pollutant is kept
                        SCASPOS = UCASIDX( UCASPOS ) + J
                        IF( ITKEEPA( SCASIDX( SCASPOS ) ) ) THEN
                            POLNAM = ITNAMA( SCASIDX( SCASPOS ) ) 
                            POLFAC = ITFACA( SCASIDX( SCASPOS ) )
                        ELSE
                            CYCLE
                        END IF
                    END IF

C.....................  Find code corresponding to current pollutant
                    POLCOD = INDEX1( POLNAM, MXIDAT, INVDNAM )
                    IF( POLCOD == 0 ) THEN
                        WRITE( MESG,94010 ) 'ERROR: Unknown  ' //
     &                      'pollutant ' // TRIM( POLNAM ) // 
     &                      ' at line', IREC
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                
C.....................  Store data in unsorted order
                    SP = SP + 1
                
                    IF( SP <= NRAWBP ) THEN
                        INRECA  ( SP ) = CURSRC    ! map to sorted source number
                        INDEXA  ( SP ) = SP        ! index for sorting POLVLA
                        TPFLGA  ( SP ) = TPF       ! temporal flag
                        IPOSCODA( SP ) = POLCOD    ! pollutant code
                        ICASCODA( SP ) = UCASPOS   ! CAS number (set to 0 for non-toxic sources)
                    
ccs                    POLVLA( SP,NEM ) = INVDCNV( J ) * EANN
                        IF( EANN >= 0. ) THEN
                            POLVLA( SP,NEM ) = EANN * POLFAC
                        ELSE
                            POLVLA( SP,NEM ) = EANN
                        END IF
                        
                        IF( EOZN >= 0. ) THEN
                            POLVLA( SP,NOZ ) = EOZN * POLFAC
                        ELSE
                            POLVLA( SP,NOZ ) = EOZN
                        END IF
                    
                        IF( CATEGORY == 'AREA' .OR. 
     &                      CATEGORY == 'POINT'     ) THEN
                            POLVLA( SP,NEF ) = EMFC
                            POLVLA( SP,NCE ) = CEFF
                            POLVLA( SP,NRE ) = REFF
                            POLVLA( SP,NRP ) = RPEN
                        END IF
                    END IF

                END DO  ! end loop through pols per CAS number

            END DO  ! end loop through pols per line
            
        END DO  ! end loop through records array

C.........  Abort if there was an error
        IF( EFLAG ) THEN
            MESG = 'Error reading data from inventory file ' //
     &              TRIM( FNAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
            
C.........  Sort inventory data by source
        CALL M3MESG( 'Sorting inventory data by source ' //
     &               'and pollutant...' )
        
        CALL SORTI2( SP, INDEXA, INRECA, IPOSCODA )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

        END SUBROUTINE RDINVDATA
