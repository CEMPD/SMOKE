
        SUBROUTINE GETISIZE( FDEV, CATEGORY, INVFMT, NREC, NRECDAT )

C**************************************************************************
C  subroutine body starts at line 92
C
C  DESCRIPTION:
C      This subroutine returns an approximate dimension needed for reading in
C      inventory files to be used to allocate memory for the unsorted 
C      inventory records.
C
C  PRECONDITIONS REQUIRED:
C      File opened on unit FDEV
C      Source category of interest CATEGORY specified
C      Inventory format INVFMT specified
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, CHECKMEM, RDLINES
C      Functions: I/O API functions, GETFLINE, GETFORMT
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 12/98
C
C**************************************************************************
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        INTEGER         GETFLINE
        INTEGER         GETFORMT
        INTEGER         JUNIT

        EXTERNAL        CRLF, GETFLINE, GETFORMT, JUNIT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN)  :: FDEV     !  unit number of input file
        CHARACTER(*), INTENT (IN)  :: CATEGORY !  desc of src category
        INTEGER     , INTENT (IN)  :: INVFMT   !  inven format code
        INTEGER     , INTENT (OUT) :: NREC     !  no. input recs
        INTEGER     , INTENT (OUT) :: NRECDAT  ! no. recs times data vars

C...........   Contents of file 
        CHARACTER*300,ALLOCATABLE:: LSTSTR( : )! Char strings in PTINV file

C...........   File units and logical/physical names
        INTEGER         TDEV        !  tmp emissions file for multiple input

C...........   Other local variables
        INTEGER         I, J, L1, L2, N1, N2   !  counters and indices

        INTEGER         CATLEN      !  string length of category
        INTEGER         IOS         !  i/o/ status
        INTEGER         FILFMT      !  file format code
        INTEGER         NLINE       !  number of lines

        CHARACTER*300   INFILE      !  input file line buffer
        CHARACTER*300   LINE        !  input file line buffer
        CHARACTER*300   MESG        !  message buffer

        CHARACTER*16 :: PROGNAME = 'GETISIZE' ! program name

C***********************************************************************
C   begin body of subroutine GETISIZE

        CATLEN = LEN_TRIM( CATEGORY )

        IF( INVFMT .EQ. LSTFMT ) THEN

C.............  Get number of lines of inventory file in list format
            MESG = CATEGORY( 1:CATLEN ) // 
     &             ' inventory file in list format'

            NLINE = GETFLINE( FDEV, MESG )

C.............  Allocate memory for storing contents of list-formatted inv file
            ALLOCATE( LSTSTR( NLINE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LSTSTR', PROGNAME )

C.............  Store lines of inventory file in list format
            CALL RDLINES( FDEV, MESG, NLINE, LSTSTR )

C.............   Initialize counters for input records
            NREC    = 0
            NRECDAT = 0

C.............  Loop through lines of list file.  Must use this loop structure
C               to be able to change J NREC for EMS95 format
            J = 0
            DO
                J = J + 1

                IF( J .GT. NLINE ) EXIT   ! End the loop 

                LINE = LSTSTR( J )
                L1 = INDEX( LINE, 'INVYEAR' )
                L2 = LEN_TRIM( LINE )

C.................  Store path of file name (if no INVYEAR packet on this line)
                IF( L1 .LE. 0 ) THEN
                    INFILE = LINE( 1:L2 )
                ELSE
                    CYCLE  ! To head of loop
                END IF

C.................  Open INFILE
                TDEV = JUNIT()
                OPEN( TDEV, ERR=1006, FILE=INFILE, STATUS='OLD' )

C.................  Make sure read pointer is at the beginning of the file
                REWIND( TDEV )

C.................  Determine format of INFILE
                FILFMT = GETFORMT( TDEV )

C.................  For EPS2 format files inside a SMOKE list file
                IF( FILFMT .EQ. EPSFMT ) THEN

                    MESG = CATEGORY( 1:CATLEN ) // 
     &                     ' inventory file "' // LINE( 1:L2 ) // '"'
                    NREC = NREC + GETFLINE( TDEV, MESG )
                    NRECDAT = NREC

C.................  For IDA format files inside a SMOKE list file
                ELSE IF( FILFMT .EQ. IDAFMT ) THEN
                    CALL GETIDASZ( TDEV, CATEGORY, N1, N2 )
                    NREC = NREC + N1
                    NRECDAT = NRECDAT + N2

C.................  For EMS-95 format files inside a SMOKE list file
                ELSE IF( FILFMT .EQ. EMSFMT ) THEN

C.....................  Make sure that next 4 files in list are also EMSFMT
C.....................  Increment line, scan for INVYEAR, open file, check file,
C                       write message, and store unit number.
                    IF( CATEGORY .EQ. 'POINT' ) THEN

                        J = J + 1
                        LINE = LSTSTR( J )
                        INFILE = LINE( 1:LEN_TRIM( LINE ) )
                        IF( INDEX( INFILE,'INVYEAR' ) .GT. 0 ) GOTO 1007 ! Error

C......................... Close previous file, get unit, and open this file
                        CLOSE( TDEV )  
                        TDEV = JUNIT()
                        OPEN( TDEV, ERR=1006, FILE=INFILE, STATUS='OLD')
                        FILFMT = GETFORMT( TDEV )
                        IF( FILFMT .NE. EMSFMT ) GO TO 1008  ! Error
                        J = J + 3 ! Skip past other three files

                    END IF

C.....................  Treat mobile sources like the IDA format because it
C                       can have multiple data values on each line
                    IF( CATEGORY .EQ. 'MOBILE' ) THEN
                        CALL GETIDASZ( TDEV, CATEGORY, N1, N2 )
                        NREC = NREC + N1
                        NRECDAT = NRECDAT + N2

C.....................  Sum number of lines in EMS-95 emission files
                    ELSE
                        NREC = NREC + GETFLINE( TDEV, MESG )
                        NRECDAT = NREC
                    END IF

                    CLOSE( TDEV )

                ELSE  ! File format not recognized	

                    MESG = 'File format is not recognized for file. ' //
     &                     CRLF() // BLANK10 // 
     &                     INFILE( 1:LEN_TRIM( INFILE ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

            END DO   ! End of loop through list-formatted PTINV file

            DEALLOCATE( LSTSTR )

C.........  For a EPS2 format file or
C.........  For an EMS-95 format file for non-point sources
        ELSE IF( INVFMT .EQ. EPSFMT                  .OR.
     &         ( INVFMT .EQ. EMSFMT           .AND.
     &         ( CATEGORY .EQ. 'AREA'   .OR.
     &           CATEGORY .EQ. 'MOBILE'     )       )    ) THEN

            MESG = CATEGORY( 1:CATLEN ) // ' inventory file'
            NREC = GETFLINE( FDEV, MESG )
            NRECDAT = NREC

C.........  For an IDA format file
        ELSE IF( INVFMT .EQ. IDAFMT ) THEN
            CALL GETIDASZ( FDEV, CATEGORY, NREC, NRECDAT )

        ELSE
            WRITE( MESG,94010 ) 'INTERNAL ERROR: Illegal call to ' //
     &             'function ' // PROGNAME // ' with format ', INVFMT
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

        IF( NREC .EQ. 0 ) THEN
            MESG = FMTNAMES( INVFMT ) // '-formatted inventory ' //
     &             'file has no valid lines of inventory data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN


C******************  ERROR MESSAGES WITH EXIT **************************
 
C.........  Error opening raw input file
1006    WRITE( MESG,94010 ) 'ERROR at line ', J, 'of PTINV. ' // 
     &         'Could not open file:' //
     &         CRLF() // BLANK5 // INFILE( 1:LEN_TRIM( INFILE ) )
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

        END
