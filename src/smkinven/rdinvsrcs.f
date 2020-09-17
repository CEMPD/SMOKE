
        SUBROUTINE RDINVSRCS( FDEV, FNAME, NRAWBP, NRAWSRCS, ORLFLG )

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
C      Updated 2/06 by B. Baek (adding wildfire format)
C
C**************************************************************************
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

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURCA, SRCIDA, 
     &                      NSTRECS, SRCSBYREC, RECIDX

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, CATEGORY
        
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: FILFMT, LSTSTR, FIREFLAG, FF10FLAG

C.........  This module is for mobile-specific data
        USE MODMOBIL, ONLY: NSCCTBL, SCCTBL, SCCRVC, SCCMAPFLAG,
     &                      NSCCMAP, SCCMAPLIST, EXCLSCCFLAG

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)    CRLF
        INTEGER         ENVINT
        INTEGER         GETFLINE
        INTEGER         GETFORMT
        INTEGER         GETINVYR
        INTEGER         JUNIT
        INTEGER         FIND1
        INTEGER         FIND1FIRST
        INTEGER         FINDC
        LOGICAL         CHKINT
        INTEGER         STR2INT
        INTEGER         INDEX1
        LOGICAL         BLKORCMT
        LOGICAL         SETENVVAR
        INTEGER*4       GETPID   
        LOGICAL         USEEXPGEO
        LOGICAL         ENVYN

        EXTERNAL        CRLF, ENVINT, GETFLINE, GETFORMT, GETINVYR, 
     &                  JUNIT, FIND1, FIND1FIRST, FINDC, ENVYN,
     &                  CHKINT, STR2INT, INDEX1, BLKORCMT, SETENVVAR,
     &                  USEEXPGEO

C...........   SUBROUTINE ARGUMENTS
        INTEGER,      INTENT (IN) :: FDEV         ! unit no. of inv file
        CHARACTER(*), INTENT (IN) :: FNAME        ! logical name of file
        INTEGER,      INTENT(OUT) :: NRAWBP       ! no. of sources with pols/acts
        INTEGER,      INTENT(OUT) :: NRAWSRCs     ! no. of raw sources
        LOGICAL,      INTENT(OUT) :: ORLFLG       ! true: read ORL inventory

C...........   Local parameters
        INTEGER      , PARAMETER :: MXRECS = 3000000  ! maximum records per iteration
        INTEGER      , PARAMETER :: NSCSEG = 8        ! num. segments in scratch file
        INTEGER      , PARAMETER :: NSEG   = 70       ! maximum no of segments
        
C...........   Local arrays
        CHARACTER(ALLLEN3+OBRLEN3),ALLOCATABLE :: TMPCSOURC( : )   ! source information from inventory file(s)
        INTEGER,           ALLOCATABLE :: TCSRCIDX ( : )   ! index for sorting source info
        INTEGER,           ALLOCATABLE :: FRSNUMS  ( :,: ) ! triplets of file, record, and source number
        CHARACTER(SRCLEN3) SCSEGMENT( NSCSEG )   ! segments from scratch file
        CHARACTER( 40 )    SEGMENT  ( NSEG )     ! segments of line

        INTEGER,            ALLOCATABLE:: CSRCIDX  ( : )    ! index for sorting CSOURCA

C...........   File units and logical/physical names
        INTEGER         CDEV        !  scratch file

C...........   Other local variables
        INTEGER         I, J, JJ, K, KK, K1, K2, L, NP, S !  counters and indices
        INTEGER         L0, L1, L2, L3, L4, L5, L6, L7, L8, L9

        INTEGER         CSRC_LEN     !  length of source characteristics
        INTEGER         CURFMT       !  format of current inventory file
        INTEGER         CURFIL       !  current file from list formatted inventory
        INTEGER         IOS          !  i/o status
        INTEGER         INVFMT       !  inventory format code
        INTEGER         IREC         !  no. of records read
        INTEGER         ISTREC       !  no. of records stored
        INTEGER         IVT          !  vehicle type code
        INTEGER         LDEV         !  device no. for log file
        INTEGER         NSCC         !  tmp no of reference SCCs
        INTEGER         MXWARN       !  maximum number of warnings
        INTEGER         NLINE        !  number of lines in list format file
        INTEGER         NPOLPERLN    !  no. of pollutants per line of inventory file
        INTEGER         NRECPERLN    !  no. of records per line
        INTEGER      :: NWARN0= 0    !  current number of warnings
        INTEGER      :: NWARN1= 0    !  current number of warnings 1
        INTEGER      :: NWRLINE = 0  !  no. of lines in file writting to log
        INTEGER*4       PID          !  UNIX process ID at runtime
        INTEGER         ROAD         !  road class number
        INTEGER         RWT          !  roadway type
        INTEGER      :: TOTSRCS = 0  !  total number of sources
        INTEGER      :: TOTRECS = 0  !  total number of records
        INTEGER      :: NORSID = 0   !  no of Oris IDs under same EIS unit
        
        LOGICAL      :: ORSFLAG = .FALSE. ! true: process multiple ORIS units under the same EIS unit
        LOGICAL      :: EFLAG   = .FALSE. ! true: error occured
        LOGICAL      :: HDRFLAG           ! true: current line is part of header
        LOGICAL      :: LSTTIME = .FALSE. ! true: last time through 

        CHARACTER(FIPLEN3) CFIP    ! fips code
        CHARACTER(LNKLEN3) CLNK    ! link ID
        CHARACTER(VIDLEN3) CIVT    ! vehicle type ID
        CHARACTER(RWTLEN3) CROAD   ! road class no.
        CHARACTER(RWTLEN3) CRWT    ! roadway type
        CHARACTER(VTPLEN3) VTYPE   ! tmp vehicle type        
        CHARACTER(RWTLEN3+VTPLEN3) CRVC    ! tmp roadway // vehicle type
        
        CHARACTER(PLTLEN3) FCID    ! facility/plant ID
        CHARACTER(CHRLEN3) PTID, NPTID    ! point ID
        CHARACTER(CHRLEN3) SKID    ! stack ID
        CHARACTER(CHRLEN3) SGID    ! segment ID
        CHARACTER(CHRLEN3) DVID    ! device ID
        CHARACTER(CHRLEN3) PRID    ! process ID

        CHARACTER(ORSLEN3) :: CORS = ' '  ! DOE plant ID
        CHARACTER(BLRLEN3) :: BLID = ' '  ! boiler ID

        CHARACTER(SCCLEN3) TSCC    ! scc code
        CHARACTER(ALLLEN3) TCSOURC ! concatenated src (minus pollutant)
        CHARACTER(ALLLEN3+OBRLEN3) FCSOURC ! concatenated src (minus pollutant)

        CHARACTER(10)      CREC    ! record number
        CHARACTER(4)       CFIL    ! file number
        CHARACTER(300)     BUFFER  ! tmp line buffer
        CHARACTER(300)     OUTLINE ! line to write to scratch file
        CHARACTER(300)     INFILE  ! input file line buffer
        CHARACTER(500)     LINE    ! input file line buffer
        CHARACTER(300)     MESG    ! message buffer
        CHARACTER(20)      VIDFMT  ! vehicle type ID format
        CHARACTER(20)      RWTFMT  ! roadway type number format
        CHARACTER(1024)    TMPFILNAM  ! File name of tmp file

        CHARACTER(512)     PATHNM           ! path name for tmp file
        CHARACTER(300)     TENLINES( 10 )   ! first ten lines of inventory file

        CHARACTER(16) :: PROGNAME =  'RDINVSRCS' ! program name

C***********************************************************************
C   begin body of subroutine RDINVSRCS

C.........  Get log file number for reports
        LDEV = INIT3()

C.........  Allocate tmp arrays for storing source information
        ALLOCATE( TMPCSOURC( MXRECS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPCSOURC', PROGNAME )
        ALLOCATE( TCSRCIDX( MXRECS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TCSRCIDX', PROGNAME )
        ALLOCATE( FRSNUMS( MXRECS,3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FRSNUMS', PROGNAME )

        TMPCSOURC = ' ' ! array
        TCSRCIDX  = 0   ! array
        FRSNUMS   = 0   ! array

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

C.........  Handle multiple ORIS units under same EIS Unit
        MESG = 'Process multiple ORIS/Boiler units under the same source unit'
        ORSFLAG = ENVYN( 'PROCESS_MULT_ORIS_UNITS_YN', MESG, .FALSE., IOS )

C.........  Get temporary directory location
        MESG = 'Path where temporary import file will be written'
        CALL ENVSTR( 'SMK_TMPDIR', MESG, '.', PATHNM, IOS )

        IF( IOS /= 0 ) THEN
            IF( NWARN0 < MXWARN ) THEN
                MESG = 'WARNING: Temporary input file will be ' //
     &                 'placed in executable directory because ' // 
     &                 CRLF() // BLANK10 // 'SMK_TMPDIR environment '//
     &                 'variable is not set properly'
                CALL M3MSG2( MESG )
                NWARN0 = NWARN0 + 1
            END IF
        END IF

C.........  Get process ID for using in tmp file name
        PID = GETPID()

C.........  Build tmp file name and set environment variable to its value,
C           so calling script can clean up tmp file if desired.
        WRITE( TMPFILNAM, '(A,I8)') TRIM( PATHNM )// '/import_tmp_', PID
        IF ( .NOT. SETENVVAR( 'SMKINVEN_TMPFILE', TMPFILNAM )) THEN
            MESG = 'Could not set environment variable for Smkinven '//
     &            'temporary file name:'// CRLF()// BLANK10// TMPFILNAM
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
C.........  Initialize ORL (fire) flag to false
        ORLFLG = .FALSE.
        FIREFLAG = .FALSE.

C.........  Create formats for mobile data
        IF( CATEGORY == 'MOBILE' ) THEN
            WRITE( VIDFMT, '("(I",I2.2,")")' ) VIDLEN3
            WRITE( RWTFMT, '("(I",I2.2,")")' ) RWTLEN3
        END IF
        
C.........  Determine file format of inventory file
        INVFMT = GETFORMT( FDEV, -1 )

C.........  If SMOKE list format, read file and check file for formats.
C           NOTE- LSTFMT defined in EMCNST3.EXT
        IF( INVFMT == LSTFMT ) THEN

C.............  Generate message for GETFLINE and RDLINES calls
            MESG = TRIM( CATEGORY ) // ' inventory file, ' //
     &             TRIM( FNAME ) // ', in list format'

C.............  Get number of lines of inventory files in list format
            NLINE = GETFLINE( FDEV, MESG )

C.............  Allocate memory for storing contents of list-formatted file
            ALLOCATE( FILFMT( NLINE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FILFMT', PROGNAME )
            ALLOCATE( LSTSTR( NLINE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LSTSTR', PROGNAME )

            FILFMT = -1  ! array
            LSTSTR = ' ' ! array

C.............  Store lines of PTINV file
            CALL RDLINES( FDEV, MESG, NLINE, LSTSTR )

C.............  Reset number of lines to remove blanks
C               (RDLINES does not store blank lines)
            DO I = 1, NLINE
                IF( LSTSTR( I ) == ' ' ) THEN
                    NLINE = I - 1
                    EXIT
                END IF    
            END DO

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

        CURFIL = 1

C.........  If file is list format, open first file
        IF( INVFMT == LSTFMT ) THEN

            LINE = LSTSTR( CURFIL )
                
C.............  Skip #LIST lines (must be first)
            IF( INDEX( LINE, 'LIST' ) > 0 ) THEN
                CURFIL = CURFIL + 1
                LINE = LSTSTR( CURFIL )
            END IF

C.............  Check for inventory year packet
            IF( GETINVYR( LINE ) > 0 ) THEN
                CURFIL = CURFIL + 1  ! move to next file in list
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

            CURFMT = FILFMT( CURFIL )

C.............  Set default inventory characteristics that depend on file format
            CALL INITINFO( CURFMT )

        ELSE

C.............  If not list format, set current format to inventory format
            CURFMT = INVFMT
        
        END IF
        
        IREC = 0   ! current record number

C.........  Open scratch file for writing record numbers
        CDEV = JUNIT()
        OPEN( CDEV, FILE=TMPFILNAM, IOSTAT=IOS, STATUS='NEW' )
        
        IF( IOS /= 0 ) THEN
            MESG = 'Could not open temporary import file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
C.........  Loop over files and multiples of MXRECS
        DO

C...........  Reset counters        
          S = 0       ! source number
          ISTREC= 0   ! number of records stored

C...........  Loop through records in current file
          DO
              IF( ISTREC == MXRECS ) EXIT
              
              READ( FDEV, 93000, IOSTAT=IOS ) LINE
            
              IREC = IREC + 1          
              IF( IOS > 0 ) THEN
                  EFLAG = .TRUE.
                  WRITE( MESG, 94010 ) 'I/O error', IOS,
     &               'reading inventory file at line', IREC
                  CALL M3MESG( MESG )
                  CYCLE
              END IF

C...............  Check if we've reached the end of the file            
              IF( IOS < 0 ) THEN

C...................  If list format, try to open next file
                  IF( INVFMT == LSTFMT ) THEN

C.......................  Close current file and reset counter
                       CLOSE( FDEV )
                      IREC = 0

C.......................  Advance to next file                 
                      CURFIL = CURFIL + 1

C.......................  Check if there are more files to read
                      IF( CURFIL <= NLINE ) THEN 
                          LINE = LSTSTR( CURFIL )

C...........................  Check for #LIST line
                          IF( INDEX( LINE, 'LIST' ) > 0 ) THEN
                              CURFIL = CURFIL + 1  ! move to next file in list
                          END IF

C...........................  Make sure current line is not INVYEAR packet                    
                          IF( GETINVYR( LINE ) > 0 ) THEN
                              CURFIL = CURFIL + 1  ! move to next file in list
                          END IF

C...........................  Make sure there are still files to read                            
                          IF( CURFIL > NLINE ) THEN
                              LSTTIME = .TRUE.
                              EXIT
                          END IF
                             
                          INFILE = LSTSTR( CURFIL )
                                        
                          OPEN( FDEV, FILE=INFILE, STATUS='OLD', 
     &                          IOSTAT=IOS )

C...........................  Check for errors while opening file
                          IF( IOS /= 0 ) THEN
                        
                              WRITE( MESG,94010 ) 'Problem at line ', 
     &                           CURFIL, 'of ' // TRIM( FNAME ) // 
     &                           '.' // ' Could not open file:' //
     &                           CRLF() // BLANK5 // TRIM( INFILE ) 
                              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
          
                          ELSE
                              WRITE( MESG,94010 ) 
     &                          'Successful OPEN for ' //
     &                          'inventory file(s):' // CRLF() // 
     &                          BLANK5 // TRIM( INFILE )
                              CALL M3MSG2( MESG ) 
            
                          END IF

C...........................  Set default inventory characteristics that depend on file format
                          CALL INITINFO( FILFMT( CURFIL ) )
                          CURFMT = FILFMT( CURFIL )
                          NWRLINE = 0
                  
C...........................  Skip back to the beginning of the loop
                          CYCLE
                  
C.......................  Otherwise, no more files to read, so exit
                      ELSE
                          LSTTIME = .TRUE.
                          EXIT
                      END IF

C...................  Otherwise, not a list file, so exit
                  ELSE
                      LSTTIME = .TRUE.
                      EXIT
                  END IF
                 
              END IF   ! end check for end of file

C...............  Skip blank lines
              IF( LINE == ' ' ) CYCLE

C...............  Check if format works with expanded geographic codes
              IF( USEEXPGEO() ) THEN
              IF( CURFMT == ORLFMT .OR. CURFMT == ORLNPFMT .OR. CURFMT == ORLFIREFMT ) THEN
                  MESG = 'ERROR: Expanded geographic codes are only ' //
     &                   'supported for inventories in FF10 format.'
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              END IF
              END IF

C...............  Process line depending on file format and source category
              SELECT CASE( CURFMT )

                CASE( MEDSFMT )
                    SELECT CASE( CATEGORY )
                    CASE( 'POINT' )   ! used for pregridded MEDS format inv
                        CALL RDSRCMEDSPT( LINE, CFIP, FCID, PTID, SKID,
     &                                   SGID, TSCC, NPOLPERLN, HDRFLAG, 
     &                                   EFLAG )
                    END SELECT

                CASE( FF10FMT )
                    ORLFLG = .TRUE.
                    FF10FLAG = .TRUE.

                    SELECT CASE( CATEGORY )
                    CASE( 'AREA' )   ! used for nonroad only
                        CALL RDSRCFF10AR( LINE, CFIP, TSCC, NPOLPERLN,
     &                                   HDRFLAG, EFLAG )
                    CASE( 'MOBILE' )
                        CALL RDSRCFF10MB( LINE, CFIP, CLNK, TSCC,
     &                                   NPOLPERLN, HDRFLAG, EFLAG )
                    CASE( 'POINT' )
                        CALL RDSRCFF10PT( LINE, CFIP, FCID, PTID, SKID,
     &                                   SGID, TSCC, CORS, BLID, NPOLPERLN,
     &                                   HDRFLAG, EFLAG )
                    END SELECT

                CASE( ORLFMT )
                    ORLFLG = .TRUE.

                    SELECT CASE( CATEGORY )
                    CASE( 'AREA' )   ! used for nonroad only
                        CALL RDSRCORLAR( LINE, CFIP, TSCC, NPOLPERLN,
     &                                   HDRFLAG, EFLAG )
                    CASE( 'MOBILE' )
                        CALL RDSRCORLMB( LINE, CFIP, CLNK, TSCC,
     &                                   NPOLPERLN, HDRFLAG, EFLAG )
                    CASE( 'POINT' )
                        CALL RDSRCORLPT( LINE, CFIP, FCID, PTID, SKID,
     &                                   SGID, TSCC, CORS, BLID, NPOLPERLN,
     &                                   HDRFLAG, EFLAG )
                    END SELECT

              CASE( ORLNPFMT )
                  ORLFLG = .TRUE.
                 
                  CALL RDSRCORLNP( LINE, CFIP, TSCC, NPOLPERLN,
     &                             HDRFLAG, EFLAG )

              CASE( ORLFIREFMT )
                  ORLFLG   = .TRUE.
                  FIREFLAG = .TRUE.

                  CALL RDSRCORLFR( LINE, CFIP, FCID, PTID, SKID,
     &                             SGID, TSCC, NPOLPERLN, HDRFLAG, 
     &                             EFLAG )

              CASE DEFAULT
                  WRITE( MESG, 94010 ) 'Routine rdinvsrc.f not '//
     &                   'expecting to read file of format code', 
     &                    CURFMT
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

              END SELECT

C...............  Check for header lines
              IF( HDRFLAG ) THEN
                  CYCLE
              END IF

C...............  SCC mapping loop : Mobile activity inventory use only.
              KK = 0
              NSCC = 0
              IF( CATEGORY == 'MOBILE' ) THEN
              IF( SCCMAPFLAG ) THEN
                  CALL PADZERO( TSCC )
                  KK   = INDEX1( TSCC, NSCCMAP, SCCMAPLIST( :,1 ) )
                  IF( KK > 0 ) THEN
                      NSCC = STR2INT( SCCMAPLIST( KK,3 ) )
                  ELSE
                      IF( EXCLSCCFLAG ) THEN
                          MESG = 'WARNING: Dropping SCC "' //
     &                         TRIM( TSCC ) //
     &                        '" not listed in SCCXREF file'
                          CALL M3MESG( MESG )
                          CYCLE  ! skip SCC not found in SCCXREF file
                      END IF
                  END IF
              END IF
              END IF

C.................  loop over mapped SCC
              DO JJ = 0, NSCC

                IF( JJ > 0 .AND. KK > 0 ) IREC = IREC + 1     ! increment no of records by reference SCCs
                IF( SCCMAPFLAG .AND. KK > 0 ) TSCC = SCCMAPLIST( KK+JJ,2 )

C.................  Write first ten lines of inventory to log file
                IF( NWRLINE < 10 .AND. .NOT. FIREFLAG ) THEN
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

C.................  Check that source characteristics are correct          
C.................  Make sure some emissions are kept for this source
                IF( NPOLPERLN == 0 ) THEN
                    IF( NWARN0 < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: No kept '//
     &                         'pollutants found at line', IREC, '. ' //
     &                         'The source will be dropped.'
                        CALL M3MESG( MESG )
                        NWARN0 = NWARN0 + 1
                    END IF
                    CYCLE
                END IF

                IF( .NOT. USEEXPGEO() .AND. .NOT. CHKINT( CFIP ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: State and/or ' //
     &                     'county code is non-integer at line', IREC
                    IF( NWARN1 < MXWARN ) CALL M3MESG( MESG )
                    NWARN1 = NWARN1 + 1
                END IF

                IF( .NOT. USEEXPGEO() .AND.
     &              ( CFIP( FIPEXPLEN3+2:FIPEXPLEN3+3 ) == '00' .OR.
     &                CFIP( FIPEXPLEN3+4:FIPEXPLEN3+6 ) == '000'     ) ) THEN
                    WRITE( MESG,94010 ) 'WARNING: State and/or ' //
     &                     'county code is zero (missing) at line', IREC
                    IF( NWARN1 < MXWARN ) CALL M3MESG( MESG )
                    NWARN1 = NWARN1 + 1
                END IF

C.................  Check source specific characteristics
                IF( CATEGORY == 'AREA' ) THEN
                
C.....................  Make sure SCC is at least 8 characters long
                    IF( LEN_TRIM( TSCC ) < 8 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: SCC code must ' //
     &                      'be at least 8 characters long at line', 
     &                      IREC
                        CALL M3MESG( MESG )
                    END IF
                
                END IF
                
                IF( CATEGORY == 'MOBILE' ) THEN

C.........................  Set vehicle type and road class
                    CALL PADZERO( TSCC )

                    IVT = STR2INT( TSCC( 13:14 ) )
                    RWT = STR2INT( TSCC( 15:16 ) )

                ELSE IF( CATEGORY == 'POINT' ) THEN
                
                    IF( CURFMT == ORLFMT .OR. CURFMT == MEDSFMT .OR.
     &                  CURFMT == ORLFIREFMT .OR. CURFMT == FF10FMT ) THEN

C.........................  Make sure SCC is at least 8 characters long
                        IF( LEN_TRIM( TSCC ) < 8 ) THEN
                            IF( NWARN0 < MXWARN ) THEN
                                WRITE( MESG,94010 ) 'WARNING: SCC ' //
     &                             'code is less than 8 characters ' //
     &                             'long at line', IREC, '. Adding ' //
     &                             'trailing zeros.'
                                CALL M3MESG( MESG )
                                NWARN0 = NWARN0 + 1
                            END IF
                            
                            DO I = LEN_TRIM( TSCC )+1, 8
                                TSCC( I:I ) = '0'
                            END DO
                        END IF

C.........................  Make sure we have a facility/plant ID
                        IF( FCID == ' ' ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Missing ' //
     &                          'plant/facility ID code at line', IREC
                            CALL M3MESG( MESG )
                        ELSE IF( LEN_TRIM( FCID ) > PLTLEN3 ) THEN
                            IF( NWARN0 < MXWARN ) THEN
                                WRITE( MESG,94010 ) 'WARNING: Facility ' //
     &                             'ID is longer than maximum 20 characters ' //
     &                             'long at line', IREC
                                CALL M3MESG( MESG )
                                NWARN0 = NWARN0 + 1
                            END IF
                        END IF

                    END IF
                
                END IF

C.................  Skip rest of loop if an error has occured
                IF( EFLAG ) CYCLE

C.................  Build concatenated source information
                SELECT CASE( CATEGORY )
                CASE( 'AREA' )
                    CALL PADZERO( TSCC )
                
                    CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3, 
     &                            CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                            CHRBLNK3, TCSOURC )
                CASE( 'MOBILE' )
                    CALL FLTRNEG( CLNK )
                    WRITE( CRWT,RWTFMT ) RWT
                    WRITE( CIVT,VIDFMT ) IVT

                    CALL BLDCSRC( CFIP, CRWT, CLNK, CIVT, TSCC, 
     &                            CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                            TCSOURC )

                CASE( 'POINT' )
                    CALL PADZERO( TSCC )
                
                    CALL BLDCSRC( CFIP, FCID, PTID, SKID, SGID, 
     &                            TSCC, CHRBLNK3, CHRBLNK3, 
     &                            TCSOURC )
                END SELECT
                
                CSRC_LEN = LEN_TRIM( TCSOURC )

C.................  Append DOE plant ORIS and Boiler IDs to define src
                FCSOURC = TCSOURC // CORS // BLID

C.................  Store source info on first time through
                IF( S == 0 ) THEN
                    S = S + 1
                    TCSRCIDX ( S ) = S
                    TMPCSOURC( S ) = FCSOURC
                END IF

C.................  On subsequent passes, only store source info 
C                   if it does not match previous source
                IF( FCSOURC /= TMPCSOURC( S ) ) THEN
                    S = S + 1
                    TCSRCIDX ( S ) = S
                    TMPCSOURC( S ) = FCSOURC
                END IF

C.................  Store current source number for this record
                ISTREC = ISTREC + 1
                FRSNUMS( ISTREC,1 ) = CURFIL
                FRSNUMS( ISTREC,2 ) = IREC
                FRSNUMS( ISTREC,3 ) = S

C.................  Update total number of sources with pollutants
                NRAWBP = NRAWBP + NPOLPERLN

              END DO  ! loop through SCCMAPLIST (NSCC)

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
                    TCSOURC = TMPCSOURC( J )( 1:CSRC_LEN )
                    FCSOURC = TMPCSOURC( J )
                    NORSID = 0
                    NRECPERLN = 0
                    TOTSRCS = TOTSRCS + 1
                END IF

C.................  If current source does not match previous, write old output line
C                   and start new line
                IF( TMPCSOURC( J )( 1:CSRC_LEN ) /= FCSOURC( 1:CSRC_LEN) ) THEN
                    WRITE( CDEV, '(A)' ) TRIM( OUTLINE )
                    OUTLINE = TMPCSOURC( J )( 1:CSRC_LEN )
                    TCSOURC = TMPCSOURC( J )( 1:CSRC_LEN )
                    FCSOURC = TMPCSOURC( J )
                    NORSID = 0                 ! reset number of oris/boiler IDs
                    NRECPERLN = 0              ! reset number of read records
                    TOTSRCS = TOTSRCS + 1      ! increment total number of sources

                ELSE
C.....................  Added new source when there are more than one
C                       ORIS and Boilers under samee plant ID
                    IF( ORSFLAG ) THEN
                    IF( TMPCSOURC( J ) /= FCSOURC ) THEN 
                        NORSID = NORSID + 1
                        CFIP = TMPCSOURC( J )( PTBEGL3( 1 ):PTENDL3( 1 ) )
                        FCID = TMPCSOURC( J )( PTBEGL3( 2 ):PTENDL3( 2 ) )
                        PTID = TMPCSOURC( J )( PTBEGL3( 3 ):PTENDL3( 3 ) )
                        SKID = TMPCSOURC( J )( PTBEGL3( 4 ):PTENDL3( 4 ) )
                        SGID = TMPCSOURC( J )( PTBEGL3( 5 ):PTENDL3( 5 ) )
                        TSCC = TMPCSOURC( J )( PTBEGL3( 6 ):PTENDL3( 6 ) )

C.........................  Update Unit ID with ## when multiple oris IDs
                        PTID = ADJUSTL( PTID )
                        L1 = LEN_TRIM( PTID )
                        IF( L1 > CHRLEN3-3 ) L1 = CHRLEN3 - 3
                        WRITE( NPTID, '( A,A,I2.2)' ) PTID( 1:L1 ),'_',NORSID

                        MESG = 'WARNING: Multiple ORIS IDs under same '//
     &                         'Unit ID: "' // TRIM(PTID) //'"'// CRLF() //BLANK10 //
     &                         'Renamed original Unit ID "'//TRIM(PTID)//'" to "'// TRIM(NPTID) //'"'
                        CALL M3MSG2( MESG )

                        CALL BLDCSRC( CFIP, FCID, NPTID, SKID, SGID,
     &                                TSCC, CHRBLNK3, CHRBLNK3,
     &                                TCSOURC )

                        WRITE( CDEV, '(A)' ) TRIM( OUTLINE )
                        OUTLINE = TCSOURC( 1:CSRC_LEN )
                        FCSOURC = TMPCSOURC( J )   ! preserve FCSOURC before it gest updated/modified
                        TMPCSOURC( J ) = TCSOURC
                        NRECPERLN = 0              ! reset number of read records
                        TOTSRCS = TOTSRCS + 1      ! increment total number of sources
                    END IF
                    END IF

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
                            OUTLINE = TRIM( OUTLINE ) // ' \'
                            WRITE( CDEV, '(A)' ) TRIM( OUTLINE )
                            OUTLINE = TRIM( TMPCSOURC( J )( 1:CSRC_LEN ) )
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
            DO I = 1, NSCSEG-1

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

        CALL SORTI2( NSTRECS, RECIDX, SRCSBYREC( :,1 ),
     &               SRCSBYREC( :,2 ) )
     
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
