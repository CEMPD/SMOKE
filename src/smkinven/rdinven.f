
        SUBROUTINE RDINVEN( FDEV, FNAME, MXIPOL, INVPCOD, INVPNAM, 
     &                      FILFMT, NRAWBP, PRATIO, TFLAG )

C***********************************************************************
C  subroutine body starts at line 157
C
C  DESCRIPTION:
C      This subroutine controls reading an inventory file for any source 
C      category from one of many formats.  It determines the format and 
C      call the appropriate reader subroutines. It controls the looping 
C      through multiple files when a list-formatted file is used as input.
C
C  PRECONDITIONS REQUIRED:
C      Input file unit FDEV opened
C      Inventory pollutant list created: MXIPOL, INVPCOD, and INVPNAM
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, BLDENAMS, CHECKMEM, FMTCSRC, RDEMSPT, 
C                   RDEPSPT, RDIDAPT, RDLINES
C      Functions: I/O API functions, GETFLINE, GETFORMT, GETIDASZ, GETINVYR,
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
C...........   This module is the inventory arrays
        USE MODSOURC 

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         GETFLINE
        INTEGER         GETFORMT
        INTEGER         GETIDASZ
        INTEGER         GETINVYR
        INTEGER         GETISIZE
        INTEGER         JUNIT

        EXTERNAL        CRLF, ENVYN, GETFLINE, GETFORMT, GETIDASZ, 
     &                  GETINVYR, GETISIZE, JUNIT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV              ! unit no. of inv file
        CHARACTER(*), INTENT (IN) :: FNAME             ! logical name of file
        INTEGER     , INTENT (IN) :: MXIPOL            ! max no. inv pols
        INTEGER     , INTENT (IN) :: INVPCOD( MXIPOL ) ! 5-dig pol codes
        CHARACTER(*), INTENT (IN) :: INVPNAM( MXIPOL ) ! name of pols
        INTEGER     , INTENT(OUT) :: FILFMT            ! file format code
        INTEGER     , INTENT(OUT) :: NRAWBP            ! no.raw records x pol
        REAL        , INTENT(OUT) :: PRATIO            ! position ratio
        LOGICAL     , INTENT(OUT) :: TFLAG             ! true: PTREF output

C...........   Contents of PTFILE 
        CHARACTER*300,ALLOCATABLE:: NLSTSTR( : )! Char strings in list-fmt file

C...........   Dropped emissions
        INTEGER         NDROP             !  number of records dropped
        REAL            EDROP  ( MXIPOL ) !  total emis dropped for each pol

C...........   File units and logical/physical names
        INTEGER         EDEV( 5 )   !  up to 5 EMS-95 emissions files
        INTEGER         TDEV        !  emissions file in list format file

C...........   Other local variables
        INTEGER         I, J, L, L1, L2 !  counters and indices

        INTEGER         ERRIOS      !  error i/o stat from sub call(s)
        INTEGER         ERRREC      !  record number for error msgs
        INTEGER         INY         !  tmp inventory year
        INTEGER         IOS         !  i/o status
        INTEGER         INVFMT      !  inventory format code
        INTEGER         FLEN        !  length of FNAME string
        INTEGER         NEDIM1      !  1st dimension for sparse emis arrays
        INTEGER         NLINE       !  number of lines
        INTEGER         NRAWBP      !  actual total raw records by pollutants
        INTEGER         NRAWIN      !  total raw record-count (estimate)
        INTEGER         NRAWOUT     !  no. of valid entries in emis file(s)

        LOGICAL      :: EFLAG  = .FALSE.  ! true: error occured

        CHARACTER*16    ERFILDSC    !  desc of file creating an error from sub
        CHARACTER*300   INFILE      !  input file line buffer
        CHARACTER*300   LINE        !  input file line buffer
        CHARACTER*300   MESG        !  message buffer

        CHARACTER*16 :: PROGNAME =  'RDINVEN' ! program name

C***********************************************************************
C   begin body of subroutine RDINVEN

        FLEN   = LEN_TRIM( FNAME )

C.........  Determine file format of PTINV file
        INVFMT = GETFORMT( FDEV )

C.........   Initialize variables for keeping track of dropped emissions
        NDROP = 0
        EDROP = 0.  ! array

C.........  If SMOKE list format, read file and check file for formats.
C           NOTE: LSTFMT defined in EMCNST3.EXT
        IF( INVFMT .EQ. LSTFMT ) THEN

C.............  Generate message for GETFLINE and RDLINES calls
            MESG = CATEGORY( 1:CATLEN ) // ' inventory file, ' //
     &             FNAME( 1:FLEN ) // ', in list format'

C.............  Get number of lines of PTINV file in list format
            NLINE = GETFLINE( FDEV, MESG )

C.............  Allocate memory for storing contents of list-format'd PTINV file
            ALLOCATE( NLSTSTR( NLINE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLSTSTR', PROGNAME )

C.............  Store lines of PTINV file
            CALL RDLINES( FDEV, MESG, NLINE, NLSTSTR )

C.............  Check the format of the list-formatted inventory file and
C               return the code for the type of files it contains
            CALL CHKLSTFL( NLINE, FNAME, NLSTSTR, FILFMT )

C.........  If not list format, then set FILFMT to the type of file (IDA,EPS)
        ELSE

            FILFMT = INVFMT
 
        ENDIF

C.........  Set default inventory characteristics (declared in MODINFO) used
C           by the IDA and EPS formats, including NPPOL
        CALL INITINFO( FILFMT )

C.........  Get the total number of records (Srcs x Non-missing pollutants)
C.........  Also, make sure file format is known

        IF( FILFMT .EQ. IDAFMT ) THEN

            NRAWIN = GETIDASZ( FDEV, CATEGORY, 1 ) ! No. actual records
            NEDIM1 = GETIDASZ( FDEV, CATEGORY, 2 ) ! No. actual records x pols

        ELSEIF( FILFMT .EQ. EPSFMT .OR. 
     &          FILFMT .EQ. EMSFMT      ) THEN

            NRAWIN = GETISIZE( FDEV, CATEGORY, INVFMT ) ! Estimate
            NEDIM1 = NRAWIN

        ELSE
            WRITE( MESG,94010 ) 'INTERNAL ERROR: File format ', FILFMT, 
     &             'not known by program "' // PROGNAME // '"'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        ENDIF 

C.........  Set ratio of NRAWIN and NEDIM1 to use in saving file source arrays
C           which are in different data structures for EMS-95 ibput than IDA
C           input.
        PRATIO = REAL( NRAWIN ) / MAX( 1,NEDIM1 )

C.........  Allocate memory for (unsorted) input arrays using dimensions set
C           based on the source category and type of inventory being input
        CALL SRCMEM( CATEGORY, 'UNSORTED', .TRUE., .FALSE., NRAWIN, 
     &               NEDIM1, NPPOL )

        CALL SRCMEM( CATEGORY, 'UNSORTED', .TRUE., .TRUE., NRAWIN, 
     &               NEDIM1, NPPOL )

C.........   Initialize sorting index
        DO I = 1, NEDIM1
            INDEXA( I ) = I
        ENDDO

C.........  Initialize pollutant-specific values as missing
        POLVLA = AMISS3  ! array

C.........  Read emissions from raw file(s) depending on input format...

C.........  IDA format (single file)
        IF( INVFMT .EQ. IDAFMT ) THEN

            SELECT CASE( CATEGORY )
            CASE( 'AREA' )
c                CALL RDIDAAR( FDEV, NRAWIN, NEDIM1, MXIPOL, INVPNAM,
c     &                        NRAWOUT, EFLAG, NDROP, EDROP )

            CASE( 'MOBILE' )
                
                CALL RDIDAMB( FDEV, NRAWIN, MXIPOL, INVPNAM, NRAWOUT, 
     &                        EFLAG, NDROP, EDROP )

            CASE( 'POINT' )
                CALL RDIDAPT( FDEV, NRAWIN, NEDIM1, MXIPOL, INVPNAM,
     &                        NRAWOUT, EFLAG, NDROP, EDROP )

            END SELECT

            NRAWBP = NRAWOUT 

C.........  EPS format (single file)
        ELSEIF( INVFMT .EQ. EPSFMT ) THEN

            SELECT CASE( CATEGORY )
            CASE( 'AREA' )
C                CALL RDEPSAR(  )

            CASE( 'MOBILE' )
c                CALL RDEPSMV(  )

            CASE( 'POINT' )
C                CALL RDEPSPT(  )

c                CBLRIDA = BLRBLNK3! Internal to rdepspt!
c                IORISA  = IMISS3

            END SELECT

C.........  SMOKE list format requires a loop for multiple files
C.........  Includes EMS-95 format
        ELSEIF( INVFMT .EQ. LSTFMT ) THEN  

            INY = IMISS3
            J   = 0
            DO             ! Loop through lines of the list-formatted file

                J = J + 1  ! Can't use standard loop because J used also below
                IF( J .GT. NLINE ) EXIT

                LINE = NLSTSTR( J )

                I = GETINVYR( LINE )

                IF( I .GT. 0 ) THEN
                    INY = I
                    CYCLE
                ENDIF

C.................  Final check to ensure the inventory year is set when needed
                IF( INY .LT. 0 .AND. FILFMT .EQ. EMSFMT ) THEN  
                    MESG = 'Must set inventory year using ' //
     &                     'INVYEAR packet for EMS-95 input.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
     
C.................  Store path of file name (if no INVYEAR packet on this line)
                INFILE = LINE

C.................  Open INFILE
                TDEV = JUNIT()
                OPEN( TDEV, ERR=1006, FILE=INFILE, STATUS='OLD' )

                WRITE( MESG,94010 ) 'Successful OPEN for ' //
     &                 'inventory file(s):' // CRLF() // BLANK5 //
     &                 INFILE( 1:LEN_TRIM( INFILE ) )
                CALL M3MSG2( MESG ) 

c NOTE: If list-format contains IDA or EPS files, there is currently no way
c       that the memory allocation is correct!

C.................  Read file based on format set above
                IF( FILFMT .EQ. IDAFMT ) THEN

                    SELECT CASE( CATEGORY )
                    CASE( 'AREA' )
c                        CALL RDIDAAR( FDEV, NRAWIN, NEDIM1, MXIPOL, 
c     &                                INVPNAM, NRAWOUT, EFLAG, 
c     &                                NDROP, EDROP )

                    CASE( 'MOBILE' )
                        CALL RDIDAMB( FDEV, NRAWIN, MXIPOL, INVPNAM, 
     &                                NRAWOUT, EFLAG, NDROP, EDROP  )


                    CASE( 'POINT' )
                        CALL RDIDAPT( FDEV, NRAWIN, NEDIM1, MXIPOL, 
     &                                INVPNAM, NRAWOUT, EFLAG, 
     &                                NDROP, EDROP )

                    END SELECT

                ELSEIF( FILFMT .EQ. EPSFMT ) THEN

                    SELECT CASE( CATEGORY )
                    CASE( 'AREA' )
C                        CALL RDEPSAR(  )

                    CASE( 'MOBILE' )
c                        CALL RDEPSMV(  )

                    CASE( 'POINT' )
C                        CALL RDEPSPT(  )

C                        CBLRIDA = BLRBLNK3  ! Internal to rdepspt!
C                        IORISA  = IMISS3

                    END SELECT

                ELSEIF( FILFMT .EQ. EMSFMT ) THEN

                    TFLAG = .TRUE.
                    EDEV( 1 ) = TDEV  ! Store first file unit number
 
C.....................  Make sure that next 4 files in list are also EMSFMT
C.....................  Increment line, scan for INVYEAR, open file, check file,
C                       write message, and store unit number.
                    DO I = 2, NEMSFILE

                        J = J + 1
                        LINE = NLSTSTR( J )
                        INFILE = LINE( 1:LEN_TRIM( LINE ) )
                        IF( INDEX( INFILE,'INVYEAR' ) .GT. 0 ) GOTO 1007 !Error
                        TDEV = JUNIT()
                        OPEN( TDEV, ERR=1006, FILE=INFILE, STATUS='OLD')
                        FILFMT = GETFORMT( TDEV )
                        IF( FILFMT .NE. EMSFMT ) GO TO 1008  ! Error
                        CALL M3MSG2( INFILE( 1:LEN_TRIM( INFILE ) ) )
                        EDEV( I ) = TDEV

                    END DO

C.....................  Call EMS-95 reader for current 5 files
C.....................  These calls populate the unsorted inventory 
C                       variables in the module MODSOURC
                    SELECT CASE( CATEGORY )
                    CASE( 'AREA' )
C                        CALL RDEMSAR(  )

                    CASE( 'MOBILE' )
c                        CALL RDEMSMV(  )

                    CASE( 'POINT' )
                        CALL RDEMSPT( EDEV, INY, NRAWIN, MXIPOL, 
     &                                INVPCOD, INVPNAM, NRAWOUT, 
     &                                ERRIOS, ERRREC, ERFILDSC, EFLAG, 
     &                                NDROP, EDROP )
 
                    END SELECT

                    IF( ERRIOS .GT. 0 ) THEN

                        L2 = LEN_TRIM( ERFILDSC )
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
     &                     INFILE( 1:LEN_TRIM( INFILE ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF

            ENDDO     ! End of loop through list-formatted PTINV file

C.............  Set exact pollutant times records
            NRAWBP = NRAWOUT

        END IF

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

C.........  Abort if there was a reading error
        IF( EFLAG ) THEN
           MESG = 'Error reading raw inventory file ' // FNAME( 1:FLEN )
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  ERROR MESSAGES WITH EXIT **************************

C.........  Error opening raw input file
1006    WRITE( MESG,94010 ) 'Problem at line ', J, 'of ' //
     &         FNAME( 1:FLEN ) // '.' // 'Could not open file:' //
     &         CRLF() // BLANK5 // INFILE( 1:LEN_TRIM( INFILE ) )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Error with INVYEAR packet read
1007    WRITE( MESG,94010 ) 'Problem at line ', J, 'of ' // 
     &         FNAME( 1:FLEN ) // '.' // CRLF() // BLANK10 // 
     &         'INVYEAR packet can be used only once for each ' //
     &         'group of five EMS-95 files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Error with number of EMS-95 files
1008    WRITE( MESG,94010 ) 'Problem at line ', J, 'of ' //
     &         FNAME( 1:FLEN ) // '.' // CRLF() // BLANK10 // 
     &        'EMS-95 files must be in groups of five.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

        END SUBROUTINE RDINVEN
