
        SUBROUTINE RDINVSRCS( FDEV, VDEV, SDEV, FNAME, 
     &                        NRAWBP, NRAWSRCS, TFLAG, TOXFLG )

C***********************************************************************
C  subroutine body starts at line 133
C
C  DESCRIPTION:
C      This subroutine controls reading an ASCII inventory file for any source 
C      category from one of many formats.  It determines the format and 
C      calls the appropriate reader subroutines. It controls the looping 
C      through multiple files when a list-formatted file is used as input.
C      This routine only reads the unique source characteristics from the inventories.
C
C  PRECONDITIONS REQUIRED:
C      Input file unit FDEV opened
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
        USE MODSOURC, ONLY: CSOURCA, SRCIDA, 
     &                      NSTRECS, SRCSBYREC, RECIDX

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, CATEGORY, NEMSFILE
        
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: FILFMT, LSTSTR

C.........  This module is for mobile-specific data
        USE MODMOBIL, ONLY: NVTYPE, NRCLAS, IVTIDLST, CVTYPLST, 
     &                      AMSRDCLS, RDWAYTYP, NSCCTBL, SCCTBL,
     &                      SCCRVC

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        INTEGER         GETFLINE
        INTEGER         GETFORMT
        INTEGER         GETINVYR
        INTEGER         GETVMIX
        INTEGER         JUNIT
        INTEGER         FIND1
        INTEGER         FIND1FIRST
        INTEGER         FINDC
        LOGICAL         CHKINT
        INTEGER         STR2INT

        EXTERNAL        CRLF, ENVINT, GETFLINE, GETFORMT, GETINVYR,  
     &                  GETVMIX, JUNIT, FIND1, FIND1FIRST, FINDC,
     &                  CHKINT, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER,          INTENT (IN) :: FDEV         ! unit no. of inv file
        INTEGER,          INTENT (IN) :: VDEV         ! unit no. of vmtmix file
        INTEGER,          INTENT (IN) :: SDEV         ! unit no. of speeds file
        CHARACTER(LEN=*), INTENT (IN) :: FNAME        ! logical name of file
        INTEGER,          INTENT(OUT) :: NRAWBP       ! no. of sources with pols/acts
        INTEGER,          INTENT(OUT) :: NRAWSRCs     ! no. of raw sources
        LOGICAL,          INTENT(OUT) :: TFLAG        ! true: PTREF output
        LOGICAL,          INTENT(OUT) :: TOXFLG       ! true: read toxics inventory

C...........   Local parameters
        INTEGER      , PARAMETER :: MXRECS = 1000000  ! maximum records per iteration
        INTEGER      , PARAMETER :: NSCSEG = 8        ! num. segments in scratch file
        
C...........   Local arrays
        CHARACTER(LEN=SRCLEN3) TMPCSOURC( MXRECS )   ! source information from inventory file(s)
        INTEGER                TCSRCIDX ( MXRECS )   ! index for sorting source info
        INTEGER                FRSNUMS  ( MXRECS,3 ) ! triplets of file, record, and source number
        CHARACTER(LEN=SRCLEN3) SCSEGMENT( NSCSEG )   ! segments from scratch file
        INTEGER, ALLOCATABLE:: CSRCIDX  ( : )        ! index for sorting CSOURCA

C...........   File units and logical/physical names
        INTEGER         EDEV( 5 )   !  up to 5 EMS-95 emissions files
        INTEGER         CDEV        !  scratch file

C...........   Other local variables
        INTEGER         I, J, K, K1, K2, L, S !  counters and indices

        INTEGER         CSRC_LEN     !  length of source characteristics
        INTEGER         CURFMT       !  format of current inventory file
        INTEGER         CURFIL       !  current file from list formatted inventory
        INTEGER         IOS          !  i/o status
        INTEGER         INVFMT       !  inventory format code
        INTEGER         IREC         !  no. of records read
        INTEGER         ISTREC       !  no. of records stored
        INTEGER         IVT          !  vehicle type code
        INTEGER         LDEV         !  device no. for log file
        INTEGER         MXWARN       !  maximum number of warnings
        INTEGER         NINVFILES    !  number of EMS-95 inventory files
        INTEGER         NLINE        !  number of lines in list format file
        INTEGER         NPOLPERLN    !  no. of pollutants per line of inventory file
        INTEGER         NRECPERLN    !  no. of records per line
        INTEGER      :: NWARN = 0    !  current number of warnings
        INTEGER      :: NWRLINE = 0  !  no. of lines in file writting to log
        INTEGER         ROAD         !  road class number
        INTEGER         RWT          !  roadway type
        INTEGER      :: TOTSRCS = 0  !  total number of sources
        INTEGER      :: TOTRECS = 0  !  total number of records
        
        LOGICAL      :: EFLAG   = .FALSE. ! true: error occured
        LOGICAL      :: EMSFLAG = .FALSE. ! true: at least one file is EMS-95 format
        LOGICAL      :: HDRFLAG           ! true: current line is part of header
        LOGICAL      :: LSTTIME = .FALSE. ! true: last time through 

        CHARACTER(LEN=FIPLEN3) CFIP    ! fips code
        CHARACTER(LEN=LNKLEN3) CLNK    ! link ID
        CHARACTER(LEN=VIDLEN3) CIVT    ! vehicle type ID
        CHARACTER(LEN=RWTLEN3) CROAD   ! road class no.
        CHARACTER(LEN=RWTLEN3) CRWT    ! roadway type
        CHARACTER(LEN=VTPLEN3) VTYPE   ! tmp vehicle type        
        CHARACTER(LEN=RWTLEN3+VTPLEN3) CRVC    ! tmp roadway // vehicle type
        
        CHARACTER(LEN=PLTLEN3) FCID    ! facility ID
        CHARACTER(LEN=CHRLEN3) PTID    ! point ID
        CHARACTER(LEN=CHRLEN3) SKID    ! stack ID
        CHARACTER(LEN=CHRLEN3) SGID    ! segment ID
        CHARACTER(LEN=CHRLEN3) DVID    ! device ID
        CHARACTER(LEN=CHRLEN3) PRID    ! process ID
        
        CHARACTER(LEN=SCCLEN3) TSCC    ! scc code
        CHARACTER(LEN=ALLLEN3) TCSOURC ! concatenated src (minus pollutant)

        CHARACTER(LEN=10)      CREC    ! record number
        CHARACTER(LEN=4)       CFIL    ! file number
        CHARACTER(LEN=300)     BUFFER  ! tmp line buffer
        CHARACTER(LEN=300)     OUTLINE ! line to write to scratch file
        CHARACTER(LEN=300)     INFILE  !  input file line buffer
        CHARACTER(LEN=300)     LINE    !  input file line buffer
        CHARACTER(LEN=300)     MESG    !  message buffer
        CHARACTER(LEN=20)      VIDFMT  ! vehicle type ID format
        CHARACTER(LEN=20)      RWTFMT  ! roadway type number format

        CHARACTER(LEN=200)     PATHNM           ! path name for tmp file
        CHARACTER(LEN=300)     TENLINES( 10 )   ! first ten lines of inventory file

        CHARACTER*16 :: PROGNAME =  'RDINVSRCS' ! program name

C***********************************************************************
C   begin body of subroutine RDINVSRCS

C.........  Get log file number for reports
        LDEV = INIT3()

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )
        
C.........  Get temporary directory location
        MESG = 'Path where temporary import file will be written'
        CALL ENVSTR( 'SMK_TMPDIR', MESG, '.', PATHNM, IOS )
        
        IF( IOS /= 0 ) THEN
            IF( NWARN < MXWARN ) THEN
                MESG = 'WARNING: Temporary input file will be ' //
     &                 'placed in executable directory because ' // 
     &                 CRLF() // BLANK10 // 'SMK_TMPDIR environment '//
     &                 'variable is not set properly'
                CALL M3MSG2( MESG )
                NWARN = NWARN + 1
            END IF
        END IF

C.........  Initialize toxics flag to false
        TOXFLG = .FALSE.

C.........  Create formats for mobile data
        IF( CATEGORY == 'MOBILE' ) THEN
            WRITE( VIDFMT, '("(I",I2.2,")")' ) VIDLEN3
            WRITE( RWTFMT, '("(I",I2.2,")")' ) RWTLEN3
        END IF
        
C.........  Determine file format of inventory file
        INVFMT = GETFORMT( FDEV )

C.........  If SMOKE list format, read file and check file for formats.
C           NOTE- LSTFMT defined in EMCNST3.EXT
        IF( INVFMT .EQ. LSTFMT ) THEN

C.............  Generate message for GETFLINE and RDLINES calls
            MESG = TRIM( CATEGORY ) // ' inventory file, ' //
     &             TRIM( FNAME ) // ', in list format'

C.............  Get number of lines of inventory files in list format
            NLINE = GETFLINE( FDEV, MESG )

C.............  Allocate memory for storing contents of list-format'd file
            ALLOCATE( FILFMT( NLINE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FILFMT', PROGNAME )
            ALLOCATE( LSTSTR( NLINE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LSTSTR', PROGNAME )

C.............  Store lines of PTINV file
            CALL RDLINES( FDEV, MESG, NLINE, LSTSTR )

C.............  Check the format of the list-formatted inventory file and
C               return the code for the type of files it contains
            CALL CHKLSTFL( NLINE, FNAME, LSTSTR, FILFMT )

C.............  Close original inventory file (will reuse device number for individual files)
            CLOSE( FDEV )

C.........  If not list format, then set FILFMT to the type of file (IDA,EPS)
        ELSE

            NLINE = 1
            ALLOCATE( FILFMT( NLINE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FILFMT', PROGNAME )
            FILFMT = INVFMT
 
        END IF

C.........  Set default inventory characteristics (declared in MODINFO) used
C           by the IDA and EPS formats, including NPPOL
        DO I = 1, NLINE
            IF( FILFMT( I ) > 0 ) THEN
                CALL INITINFO( FILFMT( I ) )
                EXIT
            END IF
        END DO

C.........  Check if any files are EMS-95 format
        DO I = 1, NLINE
            IF( FILFMT( I ) == EMSFMT ) THEN
                EMSFLAG = .TRUE.
                EXIT
            END IF
        END DO

C.........  If EMS-95 point source inventory, check for correct number of files
        IF( EMSFLAG .AND. CATEGORY == 'POINT' ) THEN
            NINVFILES = 0
            DO I = 1, NLINE
                CURFMT = FILFMT( I )
                
C.................  Skip INVYEAR packet lines
                IF( CURFMT < 0 ) CYCLE
                                
C.................  Make sure all files are EMS format
                IF( CURFMT /= EMSFMT ) THEN
                    MESG = 'Cannot use EMS-95 with other formats ' //
     &                     'for point source inventory'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Make sure file opens properly                
                INFILE = LSTSTR( I )
                OPEN( FDEV, FILE=INFILE, STATUS='OLD', IOSTAT=IOS )
                IF( IOS /= 0 ) THEN
                    WRITE( MESG,94010 ) 'Problem at line ', I, 'of ' //
     &                     TRIM( FNAME ) // '. Could not open file:' //
     &                     CRLF() // BLANK5 // TRIM( INFILE )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                CLOSE( FDEV )
                
                NINVFILES = NINVFILES + 1
            END DO
            
            IF( MOD( NINVFILES, NEMSFILE ) /= 0 ) THEN
                MESG = 'EMS-95 point source files must be in groups ' //
     &                 'of five.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

C.........  Read vehicle mix, if it is available
C.........  The tables are passed through MODMOBIL and MODXREF
        IF( VDEV .GT. 0 ) THEN
            CALL RDVMIX( VDEV )
        ELSE

C.............  Check if VMT mix is required
            IF( CATEGORY == 'MOBILE' .AND. EMSFLAG ) THEN
                MESG = 'Mobile VMT mix data are required ' //
     &                 'for import of EMS-95 mobile format' //
     &                 CRLF() // BLANK10 //
     &                 'Set the IMPORT_VMTMIX_YN environment' //
     &                 ' variable to Y and try again.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

        CURFIL = 1

C.........  If file is list format, check for inventory year packet
        IF( INVFMT == LSTFMT ) THEN
            LINE = LSTSTR( CURFIL )
            
            IF( GETINVYR( LINE ) > 0 ) THEN
                CURFIL = CURFIL + 1  ! move to next file in list
            ELSE
                IF( EMSFLAG ) THEN
                    MESG = 'Must set inventory year using ' //
     &                     'INVYEAR packet for EMS-95 input.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            END IF

C.............  If EMS point inventory, advance to emission file
            IF( EMSFLAG .AND. CATEGORY == 'POINT' ) THEN              
                CURFIL = CURFIL + 1
            END IF

C.............  Make sure there are more files
            IF( CURFIL > NLINE ) THEN
                MESG = 'No individual inventory files in list file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
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
            
            CURFMT = FILFMT( CURFIL )
            
        ELSE

C.............  If not list format, set current format to inventory format
            CURFMT = INVFMT
        
        END IF
        
        IREC = 0   ! current record number

C.........  Open scratch file for writing record numbers
        CDEV = JUNIT()
        OPEN( CDEV, FILE=TRIM( PATHNM ) // '/import_tmp', IOSTAT=IOS )
        
        IF( IOS /= 0 ) THEN
            MESG = 'Could not open temporary import file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
C.........  Loop over files and multiples of MXRECS
        DO

C.............  Reset counters        
            S = 0       ! source number
            ISTREC= 0   ! number of records stored

C.............  Loop through records in current file
            DO
                IF( ISTREC == MXRECS ) EXIT
              
                READ( FDEV, 93000, IOSTAT=IOS ) LINE
            
                IREC = IREC + 1
            
                IF( IOS > 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'I/O error', IOS,
     &                 'reading inventory file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Check if we've reached the end of the file            
                IF( IOS < 0 ) THEN

C.....................  If list format, try to open next file
                    IF( INVFMT == LSTFMT ) THEN

C.........................  Close current file and reset counter
                        CLOSE( FDEV )
                        IREC = 0

C.........................  Advance to next file
                        IF( EMSFLAG .AND. CATEGORY == 'POINT' ) THEN
                            CURFIL = CURFIL + 4
                        ELSE                    
                            CURFIL = CURFIL + 1
                        END IF

C.........................  Check if there are more files to read
                        IF( CURFIL <= NLINE ) THEN 
                            LINE = LSTSTR( CURFIL )

C.............................  Make sure current line is not INVYEAR packet                    
                            IF( GETINVYR( LINE ) > 0 ) THEN
                                CURFIL = CURFIL + 1  ! move to next file in list
                            END IF
                            
C.............................  If EMS point inventory, skip ahead to emission file
                            IF( EMSFLAG .AND. CATEGORY == 'POINT' ) THEN
                                CURFIL = CURFIL + 1
                            END IF

C.............................  Make sure there are still files to read                            
                            IF( CURFIL > NLINE ) THEN
                                LSTTIME = .TRUE.
                                EXIT
                            END IF
                             
                            INFILE = LSTSTR( CURFIL )
                                        
                            OPEN( FDEV, FILE=INFILE, STATUS='OLD', 
     &                            IOSTAT=IOS )
                    
C.............................  Check for errors while opening file
		         	        IF( IOS /= 0 ) THEN
						
					            WRITE( MESG,94010 ) 'Problem at line ', 
     &  			               CURFIL, 'of ' // TRIM( FNAME ) // 
     &  			               '.' // ' Could not open file:' //
     &		                       CRLF() // BLANK5 // TRIM( INFILE ) 
					            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
			
				            ELSE
					            WRITE( MESG,94010 ) 
     &  			              'Successful OPEN for ' //
     &		                      'inventory file(s):' // CRLF() // 
     &                            BLANK5 // TRIM( INFILE )
					            CALL M3MSG2( MESG ) 
			
					        END IF

C.............................  Set default inventory characteristics that depend on file format
				            CALL INITINFO( FILFMT( CURFIL ) )
				            CURFMT = FILFMT( CURFIL )
				            
				            NWRLINE = 0
				  
C.............................  Skip back to the beginning of the loop
                            CYCLE
				  
C.........................  Otherwise, no more files to read, so exit
				        ELSE
				            LSTTIME = .TRUE.
				            EXIT
				        END IF

C.....................  Otherwise, not a list file, so exit
                    ELSE
                        LSTTIME = .TRUE.
                        EXIT
                    END IF
                 
                END IF   ! end check for end of file
            
C.................  Skip blank lines
                IF( LINE == ' ' ) CYCLE            

C.................  Process line depending on file format and source category
                SELECT CASE( CURFMT )
                CASE( IDAFMT )
                    SELECT CASE( CATEGORY )
                    CASE( 'AREA' )
                        CALL RDSRCIDAAR( LINE, CFIP, TSCC, NPOLPERLN,
     &                                   HDRFLAG, EFLAG )
                    CASE( 'MOBILE' )
                        CALL RDSRCIDAMB( LINE, CFIP, CLNK, TSCC, 
     &                                   NPOLPERLN, HDRFLAG, EFLAG )
                    CASE( 'POINT' )
                        CALL RDSRCIDAPT( LINE, CFIP, FCID, PTID, SKID,
     &                                   SGID, TSCC, NPOLPERLN, 
     &                                   HDRFLAG, EFLAG )
                    END SELECT
                CASE( EPSFMT )
                    MESG = 'EPS 2.0 format is not currently supported'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                CASE( EMSFMT )
                    SELECT CASE( CATEGORY )
                    CASE( 'AREA' )
                        CALL RDSRCEMSAR( LINE, CFIP, TSCC, NPOLPERLN,
     &                                   HDRFLAG, EFLAG )
                    CASE( 'MOBILE' )
                        CALL RDSRCEMSMB( LINE, CFIP, CROAD, CLNK,
     &                                   NPOLPERLN, HDRFLAG, EFLAG )
                        TSCC = ' '   ! set fake SCC code
                    CASE( 'POINT' )
                        TFLAG = .TRUE. 
                        CALL RDSRCEMSPT( LINE, CFIP, FCID, SKID, DVID, 
     &                                   PRID, NPOLPERLN, 
     &                                   HDRFLAG, EFLAG )
                        TSCC = ' '   ! set fake SCC code
                    END SELECT
                CASE( NTIFMT )
                    TOXFLG = .TRUE.
                    
                    SELECT CASE( CATEGORY )
                    CASE( 'AREA' )
                        CALL RDSRCNTIAR( LINE, CFIP, TSCC, NPOLPERLN,
     &                                   HDRFLAG, EFLAG )
                    CASE( 'MOBILE' )
                        CALL RDSRCNTIMB( LINE, CFIP, CLNK, TSCC,
     &                                   NPOLPERLN, HDRFLAG, EFLAG )
                    END SELECT
                END SELECT

C.................  Check for header lines
                IF( HDRFLAG ) THEN
                    CYCLE
                END IF

C.................  Write first ten lines of inventory to log file
                IF( NWRLINE < 10 ) THEN
                    NWRLINE = NWRLINE + 1
                    TENLINES( NWRLINE ) = BLANK10 // TRIM( LINE )
                    
                    IF( NWRLINE == 10 ) THEN
                        MESG = BLANK10 // 
     &                      'First 10 lines of current inventory:'
                        WRITE( LDEV, '(A)' ) TRIM( MESG )
                        
                        DO I = 1,NWRLINE
                            WRITE( LDEV, '(A)' ) TRIM( TENLINES( I ) )
                        END DO
                    END IF
                END IF

C.................  Make sure some emissions are kept for this source
                IF( NPOLPERLN == 0 ) THEN
                    IF( NWARN < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: No kept '//
     &                     'pollutants found at line', IREC, '. ' //
     &                     'The source will be dropped.'
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                    CYCLE
                END IF

C.................  Check that source characteristics are correct
            	IF( .NOT. CHKINT( CFIP ) ) THEN
            	    EFLAG = .TRUE.
            	    WRITE( MESG,94010 ) 'ERROR: State and/or county ' //
     &      	           'code is non-integer at line', IREC
            	    CALL M3MESG( MESG )
            	END IF
            	
            	IF( CFIP( 2:3 ) == '00' .OR.
     &      	    CFIP( 4:6 ) == '000'     ) THEN
            	    EFLAG = .TRUE.
            	    WRITE( MESG,94010 ) 'ERROR: State and/or county ' //
     &      	           'code is missing at line', IREC
            	    CALL M3MESG( MESG )
            	END IF
	
c                IF( .NOT. CHKINT( TSCC ) ) THEN
c                    EFLAG = .TRUE.
c                    WRITE( MESG,94010 ) 'ERROR: SCC code is ' //
c     &      	       'non-integer at line', IREC
c                    CALL M3MESG( MESG )
c                END IF

C.................  Check source specific characteristics
                IF( CATEGORY == 'MOBILE' ) THEN

                    IF( CURFMT == IDAFMT .OR. CURFMT == NTIFMT ) THEN

C.........................  Check if SCC has proper length
                        IF( LEN_TRIM( TSCC ) /= SCCLEN3 ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: SCC code not ',
     &                         SCCLEN3, ' characters wide at line', IREC
                            CALL M3MESG( MESG )
                        END IF

C.........................  Set vehicle type and road class
                        IVT = STR2INT( TSCC( 3:6 ) )
                        ROAD = STR2INT( TSCC( 8:10 ) )
                                                
C.........................  Ensure that vehicle type is valid
                        DO J = 1, NVTYPE
                            IF( IVT == IVTIDLST( J ) ) EXIT
                        END DO
                        
                        IF( J > NVTYPE ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 
     &                         'ERROR: Vehicle type "' //
     &                         TSCC( 3:6 ) // '" at line ', IREC,
     &                         ' was not found in list of valid types'
                            CALL M3MESG( MESG )
                        END IF
                    END IF

                    IF( CURFMT == EMSFMT ) THEN

C.........................  Match FIP, road, and link with vehicle mix table
C                           Error messages are written out in GETVMIX
                        K1 = GETVMIX( CFIP, CROAD, CLNK )
                        IF( K1 <= 0 ) THEN
                            EFLAG = .TRUE.
                            CYCLE
                        END IF
                        
C.........................  Check that road class is an integer
                        IF( .NOT. CHKINT( CROAD ) ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Road class ' //
     &                         'is not an integer at line', IREC
                            CALL M3MESG( MESG )
                        END IF
                        
                        ROAD = STR2INT( CROAD )
                    END IF
                    
C.....................  Ensure that road class is valid and convert from road class
                    J = FIND1( ROAD, NRCLAS, AMSRDCLS )
                    
                    IF( J <= 0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Road class "' //
     &                         CROAD // '" at line', IREC,
     &                         ' was not found in list of valid classes'
                        CALL M3MESG( MESG )
                    ELSE
                        RWT = RDWAYTYP( J )
                    END IF

                ELSE IF( CATEGORY == 'POINT' ) THEN
                
                    IF( CURFMT == IDAFMT ) THEN
                    
C.........................  Make sure we have a plant ID
                        IF( PTID == ' ' ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Missing ' //
     &                          'plant ID code at line', IREC
                            CALL M3MESG( MESG )
                        END IF

                    END IF
                    
                    IF( CURFMT == EMSFMT ) THEN

C.........................  Check for facility, stack, device, and process IDs
!                        IF( FCID == ' ' ) THEN
!                            EFLAG = .TRUE.
!                            WRITE( MESG,94010 ) 'ERROR: Missing ' //
!     &                          'facility ID code at line', IREC
!                            CALL M3MESG( MESG )
!                        END IF
!                        
!                        IF( SKID == ' ' ) THEN
!                            EFLAG = .TRUE.
!                            WRITE( MESG,94010 ) 'ERROR: Missing ' //
!     &                          'stack ID code at line', IREC
!                            CALL M3MESG( MESG )
!                        END IF
!                        
!                        IF( DVID == ' ' ) THEN
!                            EFLAG = .TRUE.
!                            WRITE( MESG,94010 ) 'ERROR: Missing ' //
!     &                          'device ID code at line', IREC
!                            CALL M3MESG( MESG )
!                        END IF
!                        
!                        IF( PRID == ' ' ) THEN
!                            EFLAG = .TRUE.
!                            WRITE( MESG,94010 ) 'ERROR: Missing ' //
!     &                          'process ID code at line', IREC
!                            CALL M3MESG( MESG )
!                        END IF
                    
                    END IF
                
                END IF

C.................  Skip rest of loop if an error has occured
                IF( EFLAG ) CYCLE

C.................  If file format is mobile, EMS-95, make sure adding
C                   sources won't go over ISTREC limit
                IF( CATEGORY == 'MOBILE' .AND. CURFMT == EMSFMT ) THEN
                    IF( ISTREC + NVTYPE > MXRECS ) THEN
                        BACKSPACE( FDEV )
                        EXIT
                    END IF
                END IF

C.................  Build concatenated source information
                SELECT CASE( CATEGORY )
                CASE( 'AREA' )
                    CALL PADZERO( CFIP )
                    CALL PADZERO( TSCC )
                
                    CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3, 
     &                            CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                            CHRBLNK3, TCSOURC )
                CASE( 'MOBILE' )
                    CALL FLTRNEG( CLNK )
                    WRITE( CRWT,RWTFMT ) RWT
                    
                    IF( CURFMT == EMSFMT ) THEN

C.........................  Loop through vehicle types
                        DO J = 1, NVTYPE
                            IVT   = IVTIDLST( J )
                            VTYPE = CVTYPLST( J )
                            WRITE( CIVT,VIDFMT ) IVT

C.............................  Search for and assign SCC                            
                            CRVC = CROAD // VTYPE
                            K2 = FINDC( CRVC, NSCCTBL, SCCRVC )
                            
                            IF( K2 <= 0 ) THEN
                                EFLAG = .TRUE.
                                MESG = 'ERROR: Could not find SCC ' //
     &                              'for ' // CRLF() // BLANK10 //
     &                              'Road Type: ' // CRWT // 
     &                              'Vehicle Type: ' // VTYPE
                                CALL M3MESG( MESG )
                                CYCLE
                            ELSE
                                TSCC = SCCTBL( K2 )
                            END IF
                            
                            CALL BLDCSRC( CFIP, CRWT, CLNK, CIVT,
     &                                    TSCC, CHRBLNK3, CHRBLNK3,
     &                                    CHRBLNK3, TCSOURC )
                            CSRC_LEN = LEN_TRIM( TCSOURC )

C.............................  Store source info on first time through
                            IF( S == 0 ) THEN
                                S = S + 1
                                TCSRCIDX ( S ) = S
                                TMPCSOURC( S ) = TCSOURC
                            END IF

C.............................  On subsequent passes, only store source info 
C                               if it does not match previous source
                            IF( TCSOURC /= TMPCSOURC( S ) ) THEN
                                S = S + 1
                                TCSRCIDX ( S ) = S
                                TMPCSOURC( S ) = TCSOURC
                            END IF

C.............................  Store current source number for this record
                            ISTREC = ISTREC + 1
                            FRSNUMS( ISTREC,1 ) = CURFIL
                            FRSNUMS( ISTREC,2 ) = IREC
                            FRSNUMS( ISTREC,3 ) = S

C.............................  Update total number of sources with pollutants
                            NRAWBP = NRAWBP + NPOLPERLN
                        END DO
                        CYCLE
     
                    ELSE
                        WRITE( CIVT,VIDFMT ) IVT
                        CALL BLDCSRC( CFIP, CRWT, CLNK, CIVT, TSCC, 
     &                                CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                                TCSOURC )
                    END IF
                CASE( 'POINT' )
                    CALL PADZERO( TSCC )
                
                    SELECT CASE( CURFMT ) 
                    CASE( IDAFMT )
                        CALL BLDCSRC( CFIP, FCID, PTID, SKID, SGID, 
     &                                TSCC, CHRBLNK3, CHRBLNK3, 
     &                                TCSOURC )
                    CASE( EMSFMT )
                        CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID,
     &                                CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                                TCSOURC )
                    END SELECT
                END SELECT
                
                CSRC_LEN = LEN_TRIM( TCSOURC )

C.................  Store source info on first time through
                IF( S == 0 ) THEN
                    S = S + 1
                    TCSRCIDX ( S ) = S
                    TMPCSOURC( S ) = TCSOURC
                END IF

C.................  On subsequent passes, only store source info 
C                   if it does not match previous source
                IF( TCSOURC /= TMPCSOURC( S ) ) THEN
                    S = S + 1
                    TCSRCIDX ( S ) = S
                    TMPCSOURC( S ) = TCSOURC
                END IF

C.................  Store current source number for this record
                ISTREC = ISTREC + 1
                FRSNUMS( ISTREC,1 ) = CURFIL
                FRSNUMS( ISTREC,2 ) = IREC
                FRSNUMS( ISTREC,3 ) = S

C.................  Update total number of sources with pollutants
                NRAWBP = NRAWBP + NPOLPERLN
        
            END DO  ! loop through MXRECS lines

C.............  Abort if there was a reading error
            IF( EFLAG ) THEN
                MESG = 'Error reading raw inventory file ' // 
     &                  TRIM( FNAME )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Update total number of stored records
            TOTRECS = TOTRECS + ISTREC
        
C.............  Sort source info
            CALL SORTIC( S, TCSRCIDX, TMPCSOURC )

C.............  Write source info and record numbers to file
            DO I = 1, S
                J = TCSRCIDX( I )

C.................  On first time through, set up first line and counter
                IF( I == 1 ) THEN
                    OUTLINE = TMPCSOURC( J )( 1:CSRC_LEN )
                    TCSOURC = TMPCSOURC( J )
                    NRECPERLN = 0
                    TOTSRCS = TOTSRCS + 1
                END IF
                
C.................  If current source does not match previous, write old output line
C                   and start new line
                IF( TMPCSOURC( J ) /= TCSOURC ) THEN
                    WRITE( CDEV, '(A)' ) TRIM( OUTLINE )
                    OUTLINE = TMPCSOURC( J )( 1:CSRC_LEN )
                    TCSOURC = TMPCSOURC( J )   ! store source info for next comparison
                    NRECPERLN = 0              ! reset number of read records
                    TOTSRCS = TOTSRCS + 1      ! increment total number of sources
                END IF

C.................  Find source number in records array               
                K = FIND1FIRST( J, ISTREC, FRSNUMS( :,3 ) )
                
                IF( K <= 0 ) THEN
                    WRITE( MESG,94010 ) 'INTERNAL ERROR: Could not ' //
     &                 'find source number', J, 'in file and record ' //
     &                  'number array'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
                END IF

C.................  Loop over records corresponding to current source
                DO
                
C.....................  Make sure not to go outside the array                
                    IF( K > ISTREC ) EXIT
                    
                    IF( FRSNUMS( K,3 ) == J ) THEN

C.........................  If already read NSCSEG-1 records, write line with continuation
C                           character and start new line
                        IF( NRECPERLN == NSCSEG-1 ) THEN
#ifdef DBLBSLASH
                            OUTLINE = TRIM( OUTLINE ) // ' \\'
#else
                            OUTLINE = TRIM( OUTLINE ) // ' \'
#endif
                            WRITE( CDEV, '(A)' ) TRIM( OUTLINE )
                            OUTLINE = TRIM( TMPCSOURC( J ) )
                            NRECPERLN = 0
                        END IF
                            
                        WRITE( CFIL, '(I4)' ) FRSNUMS( K,1 )
                        WRITE( CREC, '(I10)' ) FRSNUMS( K,2 )

C.........................  If writing the first record, can't use trim otherwise
C                           will lose any blank source characteristics
                        IF( NRECPERLN == 0 ) THEN
                            OUTLINE = OUTLINE( 1:CSRC_LEN ) // ' ' //
     &                                TRIM( ADJUSTL( CFIL ) ) // '/' //
     &                                TRIM( ADJUSTL( CREC ) )
                        ELSE
                            OUTLINE = TRIM( OUTLINE ) // ' ' //
     &                                TRIM( ADJUSTL( CFIL ) ) // '/' // 
     &                                TRIM( ADJUSTL( CREC ) )
                        END IF
                        
                        NRECPERLN = NRECPERLN + 1
                        K = K + 1
                    ELSE
                        EXIT
                    END IF
                END DO    ! loop over records

C.................  If last source, write final line                
                IF( I == S ) THEN
                    WRITE( CDEV, '(A)' ) TRIM( OUTLINE )
                END IF
            END DO   ! loop to write sources to file

C.............  Check if this is last time
            IF( LSTTIME ) EXIT

        END DO

C.........  Rewind scratch file        
        REWIND( CDEV )
        
C.........  Allocate memory to read complete scratch file
        ALLOCATE( CSRCIDX( TOTSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )
        ALLOCATE( CSOURCA( TOTSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSOURCA', PROGNAME )
        ALLOCATE( SRCSBYREC( TOTRECS,3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCSBYREC', PROGNAME )
        ALLOCATE( RECIDX( TOTRECS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RECIDX', PROGNAME )
        
        S = 0
        ISTREC = 0
        TCSOURC = ' '
        
C.........  Read source info from scratch file
        DO
        
            READ( CDEV, 93000, IOSTAT=IOS ) LINE

C.............  Check for I/O errors            
            IF( IOS > 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &                 'reading scratch file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Check for end of file            
            IF( IOS < 0 ) EXIT

C.............  Parse line after source info into segments
            BUFFER = LINE( CSRC_LEN+1:LEN_TRIM( LINE ) )
            CALL PARSLINE( BUFFER, NSCSEG, SCSEGMENT )
            
C.............  Check if this is a continuation of a previous line
            IF( LINE( 1:CSRC_LEN ) /= TCSOURC ) THEN
                S = S + 1
            
C.................  Store information from line
                CSRCIDX( S ) = S
                CSOURCA( S ) = LINE( 1:CSRC_LEN )
                TCSOURC = CSOURCA( S )
            END IF
            
C.............  Loop through segments (file and record numbers)
            DO I = 1,NSCSEG-1

C.................  Exit if segment is blank (reached end of line)
                IF( SCSEGMENT( I ) == ' ' ) EXIT

C.................  Increment record counter and initialize sorting array
                ISTREC = ISTREC + 1
                RECIDX( ISTREC ) = ISTREC
                
C.................  Find location of separator
                K = INDEX( SCSEGMENT( I ), '/' )
                L = LEN_TRIM( SCSEGMENT( I ) )

C.................  Store file number                
                SRCSBYREC( ISTREC,1 ) = 
     &              STR2INT( SCSEGMENT( I )( 1:K-1 ) )
                SRCSBYREC( ISTREC,2 ) = 
     &              STR2INT( SCSEGMENT( I )( K+1:L ) )
                SRCSBYREC( ISTREC,3 ) = S
            END DO
        
        END DO

        NRAWSRCS = S
        NSTRECS  = ISTREC

C.........  Sort sources by record array by file number then record number
        CALL M3MESG( 'Sorting sources by file and line number...' )

C.........  If processing mobile EMS-95 format, add third dimension to sort
C           since there will be multiple sources with the same file and 
C           record number
        IF( CATEGORY == 'MOBILE' .AND. EMSFLAG ) THEN
            CALL SORTI3( NSTRECS, RECIDX, SRCSBYREC( :,1 ), 
     &                   SRCSBYREC( :,2 ), SRCSBYREC( :,3 ) )
        ELSE
            CALL SORTI2( NSTRECS, RECIDX, SRCSBYREC( :,1 ),
     &                   SRCSBYREC( :,2 ) )
        END IF
     
C.........  Sort inventory sources if needed (only if have more than MXRECS values)
        IF( NSTRECS > MXRECS ) THEN
            CALL M3MESG( 'Sorting sources by characteristics...' )
            
            CALL SORTIC( NRAWSRCS, CSRCIDX, CSOURCA )
        END IF

C.........  Allocate memory for source numbers
        ALLOCATE( SRCIDA( NRAWSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCIDA', PROGNAME )

C.........  Loop through sources to determine source IDs
        TCSOURC = ' '
        NSRC = 0

        DO I = 1, NRAWSRCS
            J = CSRCIDX( I )

            IF( CSOURCA( J ) /= TCSOURC ) THEN
                NSRC = NSRC + 1
                TCSOURC = CSOURCA( J )
            END IF
            
            SRCIDA( J ) = NSRC
        END DO

C.........  Deallocate arrays that are no longer needed
        DEALLOCATE( CSRCIDX )
     
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

94070   FORMAT( I3, A1, I8 )

        END SUBROUTINE RDINVSRCS
