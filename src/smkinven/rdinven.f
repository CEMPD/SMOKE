
        SUBROUTINE RDPTINV( FDEV, MXIPOL, INVPCOD, INVPNAM, INVFMT, 
     &                      NPSRC, TFLAG, EFLAG, INVSTAT )

C***********************************************************************
C  subroutine body starts at line 154
C
C  DESCRIPTION:
C      This subroutine controls reading a point source inventory file from
C      one of many formats.  It dertermines the format and call the appropriate
C      reader subroutines. It sorts the raw data and store it in sorted 
C      format.
C
C  PRECONDITIONS REQUIRED:
C      Input file unit FDEV opened
C      Inventory pollutant list created: MXIPOL, INVPCOD, and INVPNAM
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, BLDENAMS, CHECKMEM, FMTCSRC, RDEMSPT, 
C                   RDEPSPT, RDIDAPT, RDLINES
C      Functions: I/O API functions, GETFLINE, GETFMTPT, GETIDASZ, GETINVYR,
C         GETISIZE
C
C  REVISION  HISTORY:
C      Created 10/98 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
        USE PTMODULE           !  NOTE: includes EMCNST3.EXT

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         GETFLINE
        INTEGER         GETFMTPT
        INTEGER         GETIDASZ
        INTEGER         GETINVYR
        INTEGER         GETISIZE
        INTEGER         JUNIT
        INTEGER         STR2INT
        INTEGER         TRIMLEN

        EXTERNAL        CRLF, ENVYN, GETFLINE, GETFMTPT, GETIDASZ, 
     &                  GETINVYR, GETISIZE, JUNIT, STR2INT, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        INTEGER       FDEV              ! unit number of inv input file (in)
        INTEGER       MXIPOL            ! max no of inv pollutants (in)
        INTEGER       INVPCOD( MXIPOL ) ! 5-digit pollutant codes (in)
        CHARACTER*(*) INVPNAM( MXIPOL ) ! name of pollutants (in)
        INTEGER       INVFMT            ! inventory format code (out)
        INTEGER       NPSRC             ! actual source count (out)
        LOGICAL       TFLAG             ! output, true iff PTREF output (out)
        LOGICAL       EFLAG             ! output, error flag (out)
        INTEGER       INVSTAT( MXIPOL ) ! status (0=not in inventory) (out)

C...........   Contents of PTFILE 
        CHARACTER*256,ALLOCATABLE:: PNLSTSTR( : )! Char strings in PTINV file

C...........   Dropped emissions
        INTEGER         NDROP             !  number of records dropped
        REAL            EDROP  ( MXIPOL ) !  total emis dropped for each pol

C...........   Variables for reading dummy names of emission output
        INTEGER :: NC1 = 0 !  pos in 2nd dim of POLVLA of primary control code
        INTEGER :: NC2 = 0 !  pos in 2nd dim of POLVLA of secondary control code
        INTEGER :: NCE = 0 !  pos in 2nd dim of POLVLA of control efficiency
        INTEGER :: NEF = 0 !  pos in 2nd dim of POLVLA of emission factors
        INTEGER :: NEM = 0 !  pos in 2nd dim of POLVLA of annual emissions
        INTEGER :: NOZ = 0 !  pos in 2nd dim of POLVLA of ozone season emis
        INTEGER :: NRE = 0 !  pos in 2nd dim of POLVLA of rule effectivenss

        INTEGER                IDUMARR( 1,NPTPPOL3 ) !  integer dummy array
        CHARACTER(LEN=IOVLEN3) ENAMES ( 1,NPTPPOL3 ) !  dummy names
        CHARACTER(LEN=IODLEN3) CDUMARR( 1,NPTPPOL3 ) !  character dummy array
    
C...........   File units and logical/physical names
        INTEGER         EDEV( 5 )   !  5 EMS-95 emissions files
        INTEGER         TDEV        !  emissions file in list format file

C...........   Other local variables
        INTEGER         I, J, K, L, L1, L2, LK, LS, S  !  counters and indices

        INTEGER         ERRIOS      !  error i/o stat from sub call(s)
        INTEGER         ERRREC      !  record number for error msgs
        INTEGER         INY         !  tmp inventory year
        INTEGER         IOS         !  i/o status
        INTEGER         FILFMT      !  file format code
        INTEGER         NEDIM1      !  1st dimension for sparse emis arrays
        INTEGER         NEDIM2      !  2nd dimension for sparse emis arrays
        INTEGER         NLINE       !  number of lines
        INTEGER         NRAWBP      !  actual total raw records by pollutants
        INTEGER         NRAWIN      !  total raw record-count (estimate)
        INTEGER         NRAWOUT     !  no. of valid entries in emis file(s)
        INTEGER         PIPCOD      !  IPOSCOD of previous iteration of loop
        INTEGER         PREVFMT     !  file format code of previous iteration

        REAL            EMISI       !  inverse emissions value
        REAL            EMISN       !  new emissions value
        REAL            EMISO       !  old emissions value

        LOGICAL         DFLAG       !  input verification:  TRUE iff ERROR
        LOGICAL         FIRSTITER   !  true when on first iteration of a loop

        CHARACTER(LEN=SRCLEN3)  LSRCCHR     !  previous CSOURC
        CHARACTER(LEN=SRCLEN3)  TSRCCHR     !  tmporary CSOURC

        CHARACTER*5     TPOLPOS     !  Temporary pollutant position
        CHARACTER*16    ERFILDSC    !  desc of file creating an error from sub
        CHARACTER*300   BUFFER      !  input file line buffer
        CHARACTER*300   INFILE      !  input file line buffer
        CHARACTER*300   LINE        !  input file line buffer
        CHARACTER*300   MESG        !  message buffer

        CHARACTER*16 :: PROGNAME =  'RDPTINV' ! program name

C***********************************************************************
C   begin body of subroutine RDPTINV

C.........  Get settings from the environment
        DFLAG = ENVYN( 'RAW_DUP_CHECK',
     &                 'Check for duplicate species-records',
     &                 .FALSE., IOS )

C.........  Determine file format of PTINV file
        INVFMT = GETFMTPT( FDEV )

C.........   Initialize variables for keeping track of dropped emissions
        NDROP = 0
        EDROP = 0.  ! array

C.........  If SMOKE list format, read file and check file for formats.
C           NOTE: LSTFMT defined in EMCNST3.EXT
        IF( INVFMT .EQ. LSTFMT ) THEN

C.............  Get number of lines of PTINV file in list format
            NLINE = GETFLINE( FDEV, 'PTINV file in list format' )

C.............  Allocate memory for storing contents of list-format'd PTINV file
            ALLOCATE( PNLSTSTR( NLINE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PNLSTSTR', PROGNAME )

C.............  Store lines of PTINV file
            MESG = 'POINT inventory file in list format'
            CALL RDLINES( FDEV, MESG, NLINE, PNLSTSTR )

C.............  Loop through lines of PTINV file to check the formats
            FIRSTITER = .TRUE.
            INY = IMISS3
            DO J = 1, NLINE

                LINE = PNLSTSTR( J )

C.................  Skip INVYEAR packet 
                I = GETINVYR( LINE )
                IF( I .GT. 0 ) THEN
                    INY = I
                    CYCLE
                ENDIF

                IF( INY .LT. 0 ) THEN  ! Final check to ensure it's set !
                    MESG = 'Must set inventory year using ' //
     &                     'INVYEAR packet for list-style input.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
     
                INFILE = LINE

C.................  Open INFILE
                TDEV = JUNIT()
                OPEN( TDEV, ERR=1006, FILE=INFILE, STATUS='OLD' )

C.................  Determine format of INFILE
                FILFMT = GETFMTPT( TDEV )

C.................  If first iteration, save format, if not, make sure 
C                   that different formats are not used in same PTINV list
                IF( FIRSTITER ) THEN
                    FIRSTITER = .FALSE.
                    PREVFMT = FILFMT

                ELSEIF( FILFMT .NE. PREVFMT ) THEN
                
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 
     &                     'ERROR: In SMOKE list input, previous ' //
     &                     'file was format', PREVFMT, 
     &                     'but file at line', J, 'was format', FILFMT
                    CALL M3MESG( MESG )

                ENDIF

                CLOSE( TDEV )

            ENDDO     ! End of loop through list-formatted PTINV file

C.............  Exit if files in list-formatted file were of inconsistent type
            IF( EFLAG ) THEN
                MESG = 'Problem reading SMOKE list format'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

C.........  If not list format, then set FILFMT to the type of file (IDA,EPS)
        ELSE

            FILFMT = INVFMT
 
        ENDIF

C.........  Get the total number of records (Srcs x Non-missing pollutants)
C.........  Also, make sure file format is known

        NEDIM2 = NPTPPOL3 ! for all posible pollutant-specific fields

        IF( FILFMT .EQ. IDAFMT ) THEN

            NRAWIN = GETIDASZ( FDEV, 'POINT', 1 ) ! No. actual records
            NEDIM1 = GETIDASZ( FDEV, 'POINT', 2 ) ! No. actual records x pols

        ELSEIF( FILFMT .EQ. EPSFMT .OR. 
     &          FILFMT .EQ. EMSFMT      ) THEN

            NRAWIN = GETISIZE( FDEV, 'POINT', INVFMT ) ! Estimate
            NEDIM1 = NRAWIN

        ELSE
            WRITE( MESG,94010 ) 'File format ', FILFMT, 
     &             'not known by program "' // PROGNAME // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF 

C.........  Allocate memory for (unsorted) input arrays using dimensions set
C           based on type of inventory being input
        ALLOCATE( IFIPA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IFIPA', PROGNAME )
        ALLOCATE( ISCCA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISCCA', PROGNAME )
        ALLOCATE( ISICA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISICA', PROGNAME )
        ALLOCATE( IORISA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IORISA', PROGNAME )
        ALLOCATE( TPFLGA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TPFLGA', PROGNAME )
        ALLOCATE( INVYRA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVYRA', PROGNAME )
        ALLOCATE( IDIUA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDIUA', PROGNAME )
        ALLOCATE( IWEKA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IWEKA', PROGNAME )
        ALLOCATE( NPCNTA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NPCNTA', PROGNAME )
        ALLOCATE( SRCIDA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCIDA', PROGNAME )
        ALLOCATE( XLOCAA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'XLOCAA', PROGNAME )
        ALLOCATE( YLOCAA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'YLOCAA', PROGNAME )
        ALLOCATE( STKHTA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STKHTA', PROGNAME )
        ALLOCATE( STKDMA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STKDMA', PROGNAME )
        ALLOCATE( STKTKA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STKTKA', PROGNAME )
        ALLOCATE( STKVEA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STKVEA', PROGNAME )
        ALLOCATE( CBLRIDA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CBLRIDA', PROGNAME )
        ALLOCATE( CPDESCA( NRAWIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CPDESCA', PROGNAME )

        ALLOCATE( INDEXA( NEDIM1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )
        ALLOCATE( IPOSCOD( NEDIM1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPOSCOD', PROGNAME )
        ALLOCATE( CSOURCA( NEDIM1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSOURCA', PROGNAME )
        ALLOCATE( POLVLA( NEDIM1, NEDIM2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLVLA', PROGNAME )

C.........   Initialize sorting index
        DO I = 1, NEDIM1
            INDEXA( I ) = I
        ENDDO

C.........  Initialize pollutant-specific values as missing
        POLVLA = AMISS3  ! array

C.........  Get dummy output pollutant variable names
        CALL BLDENAMS( 'POINT', 1, NPTPPOL3, 'DUM', ENAMES, CDUMARR,
     &                 IDUMARR, CDUMARR )

C.........  Based on the order of the oiutput names, find the positions in the
C           second dimension of POLVLA and POLVAL for the pollutant-specific
C           data.
        DO I = 1, NPTPPOL3
            IF( ENAMES( 1,I )(1:IOVLEN3)  .EQ. 'DUM' ) NEM = I
            IF( ENAMES( 1,I )(1:CPRTLEN3) .EQ. OZNSEART ) NOZ = I
            IF( ENAMES( 1,I )(1:CPRTLEN3) .EQ. CTLEFFRT ) NCE = I
            IF( ENAMES( 1,I )(1:CPRTLEN3) .EQ. RULEFFRT ) NRE = I
            IF( ENAMES( 1,I )(1:CPRTLEN3) .EQ. EMISFCRT ) NEF = I
            IF( ENAMES( 1,I )(1:CPRTLEN3) .EQ. CECOD1RT ) NC1 = I
            IF( ENAMES( 1,I )(1:CPRTLEN3) .EQ. CECOD2RT ) NC2 = I
        ENDDO

        IF( NEM .EQ. 0 .OR. NOZ .EQ. 0 .OR. NCE .EQ. 0 .OR.
     &      NRE .EQ. 0 .OR. NEF .EQ. 0 .OR. NC1 .EQ. 0 .OR.
     &      NC2 .EQ. 0 ) THEN
            MESG = 'INTERNAL ERROR: Could not find position of ' //
     &             'one or more of the pollutant-specific data.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        ENDIF

C.........  Read emissions from raw file(s) depending on input format...

C.........  IDA format (single file)
        IF( INVFMT .EQ. IDAFMT ) THEN

C            CALL RDIDAPT( FDEV, INY, NRAWIN, NEDIM1, NEM, NOZ, NCE, 
C     &                    NRE, NEF, NC1, NC2, MXIPOL, INVPCOD, INVPNAM,
C     &                    EFLAG, NDROP, EDROP )

            NRAWBP = NEDIM1  ! it was exact to start with b/c from RDIDASIZ

C.........  EPS format (single file)
        ELSEIF( INVFMT .EQ. EPSFMT ) THEN

C            CALL RDEPSPT(  )
            CBLRIDA = BLRBLNK3
            IORISA  = IMISS3

C.........  SMOKE list format requires a loop for multiple files
C.........  Includes EMS-95 format
        ELSEIF( INVFMT .EQ. LSTFMT ) THEN  

            INY = IMISS3
            J   = 0
            DO             ! Loop through lines of the list-formatted file

                J = J + 1  ! Can't use standard loop because J used also below
                IF( J .GT. NLINE ) EXIT

                LINE = PNLSTSTR( J )

                I = GETINVYR( LINE )

                IF( I .GT. 0 ) THEN
                    INY = I
                    CYCLE
                ENDIF

                IF( INY .LT. 0 ) THEN  ! Final check to ensure it's set !
                    MESG = 'Must set inventory year using ' //
     &                     'INVYEAR packet for list-style input.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
     
C.................  Store path of file name (if no INVYEAR packet on this line)
                INFILE = LINE

C.................  Open INFILE
                TDEV = JUNIT()
                OPEN( TDEV, ERR=1006, FILE=INFILE, STATUS='OLD' )

                WRITE( MESG,94010 ) 'Successful OPEN for ' //
     &                 'inventory file(s):' // CRLF() // BLANK5 //
     &                 INFILE( 1:TRIMLEN( INFILE ) )
                CALL M3MSG2( MESG ) 

C.................  Read file based on format set above
                IF( FILFMT .EQ. IDAFMT ) THEN

C                   CALL RDIDAPT  ! Will add later

                ELSEIF( FILFMT .EQ. EPSFMT ) THEN

C                   CALL RDEPSPT  ! Will add later

                ELSEIF( FILFMT .EQ. EMSFMT ) THEN

                    TFLAG = .TRUE.
                    EDEV( 1 ) = TDEV  ! Store first file unit number
 
C.....................  Make sure that next 4 files in list are also EMSFMT
C.....................  Increment line, scan for INVYEAR, open file, check file,
C                       write message, and store unit number.
                    DO I = 1, 4

                        J = J + 1
                        LINE = PNLSTSTR( J )
                        INFILE = LINE( 1:TRIMLEN( LINE ) )
                        IF( INDEX( INFILE,'INVYEAR' ) .GT. 0 ) GOTO 1007 ! Error
                        TDEV = JUNIT()
                        OPEN( TDEV, ERR=1006, FILE=INFILE, STATUS='OLD')
                        FILFMT = GETFMTPT( TDEV )
                        IF( FILFMT .NE. EMSFMT ) GO TO 1008  ! Error
                        CALL M3MSG2( INFILE( 1:TRIMLEN( INFILE ) ) )
                        EDEV( I+1 ) = TDEV

                    ENDDO

C.....................  Call EMS-95 reader for current 5 files
C.....................  These calls populate the unsorted inventory 
C                       variables in the Module MDL_PTINV
                    CALL RDEMSPT( EDEV, INY, NRAWIN, NEM, NCE, NRE, 
     &                            MXIPOL, INVPCOD, INVPNAM, NRAWOUT, 
     &                            ERRIOS, ERRREC, ERFILDSC, EFLAG, 
     &                            NDROP, EDROP )

                    IF( ERRIOS .GT. 0 ) THEN

                        L2 = TRIMLEN( ERFILDSC )
                        WRITE( MESG, 94010 ) 
     &                         'Error ', ERRIOS,  'reading ' // 
     &                         ERFILDSC( 1:L2 ) // ' file at line', 
     &                         ERRREC
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                    ELSEIF( EFLAG ) THEN
                        MESG = 'ERROR reading EMS-95 inventory files.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                    ENDIF

                ELSE  ! File format not recognized	

                    MESG = 'File format is not recognized for file. ' //
     &                     CRLF() // BLANK10 // 
     &                     INFILE( 1:TRIMLEN( INFILE ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF

            ENDDO     ! End of loop through list-formatted PTINV file

C.............  Set exact pollutant times records
            NRAWBP = NRAWOUT

        ENDIF

C.........  Report how many records were dropped and the emissions involved
        IF( NDROP .GT. 0 ) THEN

            WRITE( MESG,94010 ) 'WARNING:', NDROP, 
     &             'input emissions records dropped.  This has' //
     &             CRLF() // BLANK5 //
     &             '        resulted in the following amounts of ' //
     &             'lost emissions in tons/year:'
            CALL M3MSG2( MESG )

            DO I = 1, MXIPOL

                IF( EDROP( I ) .GT. 0. ) THEN
                    WRITE( MESG,94060 ) 
     &                     BLANK16 // INVPNAM( I ) // ': ', EDROP( I )
                    CALL M3MSG2( MESG )
                ENDIF

            ENDDO

        END IF          !  if ndrop > 0

        CALL M3MSG2( 'Sorting raw inventory data...' )

C.........  Sort inventory and pollutants (sources x pollutants). Note that
C           sources are sorted based on character string definition of the 
C           source so that source definition can be consistent with that of
C           the input format.
        CALL SORTIC( NRAWBP, INDEXA, CSOURCA )

C.........  Loop through sources X pollutants to determine source IDs and check
C           for duplicates. Also keep a count of the total unique key
C           combinations (CSOURCA without the pollutant position)
C.........  NOTE: The last part of the CSOURCA string is the integer position 
C           of the pollutant for that record in the INVPNAM pollutant array 
        LSRCCHR = EMCMISS3
        LK = IMISS3
        S = 0
        DO I = 1, NRAWBP
            
            J  = INDEXA( I )
            L1 = POLPOS3 - 1

            TSRCCHR = CSOURCA( J )(    1:L1 )      ! Source characteristics
            TPOLPOS = CSOURCA( J )( L1+1:SRCLEN3 ) ! Pos of pollutant (ASCII)

C.............  Update pointer for list of actual pollutants
            K = STR2INT( TPOLPOS )  ! Convert pollutant code to integer
            INVSTAT( K ) = 1
            IPOSCOD( I ) = K
           
C.............  Increment source count by comparing this iteration to previous
            IF( TSRCCHR .NE. LSRCCHR ) THEN
                S = S + 1
                LSRCCHR = TSRCCHR

C.............  Give message of duplicates are not permitted in inventory
C.............  This IF also implies TSRCCHR = LSRCCHR
            ELSEIF( K .EQ. LK ) THEN

                CALL FMTCSRC( TSRCCHR, 8, BUFFER, L2 )

                IF ( DFLAG ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Duplicate emissions found for' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                ELSE
                    MESG = 'WARNING: Duplicate emissions found for' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                ENDIF

                CALL M3MESG( MESG )

            ENDIF

            LK = K  ! Store pollutant index for comparison in next iteration

C.............  Assign source ID (to use as an index) for all inv X pol
            SRCIDA( I ) = S

        ENDDO  ! On sources x pollutants

        NPSRC = S

        IF( EFLAG ) THEN
           MESG = 'Error in raw inventory file(s)'         
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Allocate memory for SMOKE inventory arrays (NOT sources X pollutnts)
        ALLOCATE( IFIP( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IFIP', PROGNAME )
        ALLOCATE( ISCC( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISCC', PROGNAME )
        ALLOCATE( ISIC( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISIC', PROGNAME )
        ALLOCATE( IORIS( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IORIS', PROGNAME )
        ALLOCATE( TPFLAG( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TPFLAG', PROGNAME )
        ALLOCATE( INVYR( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVYR', PROGNAME )
        ALLOCATE( IDIU( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDIU', PROGNAME )
        ALLOCATE( IWEK( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IWEK', PROGNAME )
        ALLOCATE( NPCNT( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NPCNT', PROGNAME )
        ALLOCATE( XLOCA( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'XLOCA', PROGNAME )
        ALLOCATE( YLOCA( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'YLOCA', PROGNAME )
        ALLOCATE( STKHT( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STKHT', PROGNAME )
        ALLOCATE( STKDM( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STKDM', PROGNAME )
        ALLOCATE( STKTK( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STKTK', PROGNAME )
        ALLOCATE( STKVE( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STKVE', PROGNAME )
        ALLOCATE( CBLRID( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CBLRID', PROGNAME )  
        ALLOCATE( CPDESC( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CPDESC', PROGNAME )  

        ALLOCATE( CSOURC( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSOURC', PROGNAME )

C.........  Loop through sources x pollutants to store sorted arrays for output
C           to I/O API file.
        LS = IMISS3
        L  = POLPOS3 - 1
        DO I = 1, NRAWBP

            J = INDEXA( I )
            S = SRCIDA( I )

            IF( S .NE. LS ) THEN
                LS  = S
                IFIP  ( S )  = IFIPA  ( J )
                ISCC  ( S )  = ISCCA  ( J )
                ISIC  ( S )  = ISICA  ( J )
                IORIS ( S )  = IORISA ( J )
                IDIU  ( S )  = IDIUA  ( J )
                IWEK  ( S )  = IWEKA  ( J )
                TPFLAG( S )  = TPFLGA ( J )
                INVYR ( S )  = INVYRA ( J )
                XLOCA ( S )  = XLOCAA ( J )
                YLOCA ( S )  = YLOCAA ( J )
                STKHT ( S )  = STKHTA ( J )
                STKDM ( S )  = STKDMA ( J )
                STKTK ( S )  = STKTKA ( J )
                STKVE ( S )  = STKVEA ( J )
                CBLRID( S )  = CBLRIDA( J )
                CPDESC( S )  = CPDESCA( J )
                CSOURC( S )  = CSOURCA( J )( 1:L )
            ENDIF

        ENDDO

C.........  Deallocate local memory for per-source temporal x-ref      

        DEALLOCATE( IFIPA, ISCCA, ISICA, IORISA, IDIUA, IWEKA, 
     &              NPCNTA, TPFLGA, INVYRA, XLOCAA, YLOCAA, STKHTA,
     &              STKDMA, STKTKA, STKVEA, CBLRIDA, CPDESCA, 
     &              CSOURCA, PNLSTSTR )

C.........  Allocate memory for aggregating any duplicate pol-specific data
C.........  Note that the POLVAL array that contains the pollutant-specific
C           data is dimensioned to output only one pollutant at a time. This is 
C           because we may need to be able to handle many pollutants in the 
C           future, and the memory requirements would be prohibitive if all of
C           the memory were allocated at the same time.           
        ALLOCATE( POLVAL( NRAWBP, NEDIM2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLVAL', PROGNAME )

        POLVAL = AMISS3

C.........  Initialize pollutant count for EPS and EMS-95 data structure
        IF( FILFMT .EQ. EPSFMT .OR. FILFMT .EQ. EMSFMT ) THEN
            NPCNT = 0  ! array
        ENDIF

C.........  Store pollutant-specific data in sorted order.  For EPS and EMS-95
C           formats, ensure that any duplicates are aggregated.
C.........  Aggregate duplicate pollutant-specific data (not possible 
C           for IDA format)
C.........  NOTE: We have already checked to ensure that if there are duplicate
C           emissions, they are allowed
C.........  NOTE: Pollutants are stored in output order because they've been
C           previously sorted in part based on their position in the master
C           array of output pollutants
        IF( FILFMT .EQ. IDAFMT ) THEN
            DO I = 1, NRAWBP

                J = INDEXA( I )

                DO K = 1, NEDIM2
                    POLVAL( I, K ) = POLVLA( J, K )
                ENDDO
            ENDDO

        ELSEIF( FILFMT .EQ. EPSFMT .OR. FILFMT .EQ. EMSFMT ) THEN
        
            K = 0
            PIPCOD = IMISS3  ! Previous iteration IPOSCOD 
            LS = IMISS3      ! Previous iteration S
            DO I = 1, NRAWBP

                J = INDEXA( I )
                S = SRCIDA( I )

                IF( S .NE. LS .OR. IPOSCOD( I ) .NE. PIPCOD ) THEN

C.....................  Sum up the number of pollutants by source, but do this
C                       here only, because this part of the IF statement is for
C                       new pollutants
                    NPCNT( S ) = NPCNT( S ) + 1

                    K = K + 1

                    POLVAL( K, NEM ) = POLVLA( J, NEM )
                    POLVAL( K, NCE ) = POLVLA( J, NCE )
                    POLVAL( K, NRE ) = POLVLA( J, NRE )

                    PIPCOD = IPOSCOD( I ) 

C.................  If the existing value is defined, sum with new emissions
C                   and use weighted average for control factors
C.................  No need to change NPCNT because it is already 1 for all
                ELSE

                    EMISN = POLVLA( J, NEM )
                    EMISO = POLVAL( K, NEM )
                    POLVAL( K, NEM ) = EMISO + EMISN

                    EMISI = 1. / POLVAL( K, NEM ) ! Compute inverse only once
                    POLVAL( K,NCE ) = ( POLVAL( K,NCE )*EMISO + 
     &                                  POLVLA( J,NCE )*EMISN  ) * EMISI
                    POLVAL( K,NRE ) = ( POLVAL( K,NRE )*EMISO + 
     &                                  POLVLA( J,NRE )*EMISN  ) * EMISI

                ENDIF

                LS = S

            ENDDO

        ENDIF

        DEALLOCATE( INDEXA, SRCIDA, POLVLA )
                
        RETURN

C******************  ERROR MESSAGES WITH EXIT **************************
 
C.........  Error because improper grouping of raw input files
1005    MESG = 'EMS95 input file list ended without complete set'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Error opening raw input file
1006    WRITE( MESG,94010 ) 'ERROR at line ', J, 'of PTINV. ' // 
     &         'Could not open file:' //
     &         CRLF() // BLANK5 // INFILE( 1:TRIMLEN( INFILE ) )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Error with INVYEAR packet read
1007    WRITE( MESG,94010 ) 'ERROR at line ', J, 'of PTINV.' // 
     &         CRLF() // BLANK5 // 'INVYEAR packet can be ' //
     &         'used only once for each group of five EMS-95 files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Error with number of EMS-95 files
1008    WRITE( MESG,94010 ) 'ERROR at line ', J, 'of PTINV.' //
     &         CRLF() // BLANK5 // 
     &        'EMS-95 files must be in groups of five.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

        END
