
        INTEGER FUNCTION GETISIZE( FDEV, CATEGORY, INVFMT )

C***********************************************************************
C  function body starts at line 92
C
C  DESCRIPTION:
C      This function returns an approximate dimension needed for reading in
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
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
        INTEGER         TRIMLEN

        EXTERNAL        CRLF, GETFLINE, GETFORMT, JUNIT, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        INTEGER       FDEV        !  unit number of input file
        CHARACTER*(*) CATEGORY    !  description of source category
        INTEGER       INVFMT      !  inventory format code (from EMDIMS3.EXT)

C...........   Contents of file 
        CHARACTER*300,ALLOCATABLE:: LSTSTR( : )! Char strings in PTINV file

C...........   File units and logical/physical names
        INTEGER         TDEV        !  emissions file

C...........   Other local variables
        INTEGER         I, J, L1, L2   !  CNTRers and indices

        INTEGER         CNTR       !  CNTR of records numbers
        INTEGER         IOS         !  i/o/ status
        INTEGER         FILFMT      !  file format code
        INTEGER         NLINE       !  number of lines

        CHARACTER*300   INFILE      !  input file line buffer
        CHARACTER*300   LINE        !  input file line buffer
        CHARACTER*300   MESG        !  message buffer

        CHARACTER*16 :: PROGNAME = 'GETISIZE' ! program name

C***********************************************************************
C   begin body of function GETISIZE

c note: must be updated for all source categories.  note PTINV file name below

C.........   Initialize CNTR of lines to read in
        CNTR = 0

        IF( INVFMT .EQ. LSTFMT ) THEN

C.............  Get number of lines of inventory file in list format
            MESG = CATEGORY( 1:TRIMLEN( CATEGORY ) ) // 
     &             ' inventory file in list format'

            NLINE = GETFLINE( FDEV, MESG )

C.............  Allocate memory for storing contents of list-formatted PTINV file
            ALLOCATE( LSTSTR( NLINE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LSTSTR', PROGNAME )

C.............  Store lines of inventory file in list format
            CALL RDLINES( FDEV, MESG, NLINE, LSTSTR )

C.............   Initialize CNTR of lines to read in
            CNTR = 0

C.............  Loop through lines of PTINV file.  Must use this loop structure
C               to be able to change J CNTR for EMS95 format
            J = 0
            DO
                J = J + 1

                IF( J .GT. NLINE ) EXIT   ! End the loop 

                LINE = LSTSTR( J )
                L1 = INDEX( LINE, 'INVYEAR' )
                L2 = TRIMLEN( LINE )

C.................  Store path of file name (if no INVYEAR packet on this line)
                IF( L1 .LE. 0 ) THEN
                    INFILE = LINE( 1:L2 )
                ELSE
                    CYCLE  ! To head of loop
                ENDIF

C.................  Open INFILE
                TDEV = JUNIT()
                OPEN( TDEV, ERR=1006, FILE=INFILE, STATUS='OLD' )

C.................  Determine format of INFILE
                FILFMT = GETFORMT( TDEV )

                IF( FILFMT .EQ. EPSFMT .OR. FILFMT .EQ. IDAFMT ) THEN

                    MESG = CATEGORY( 1:TRIMLEN( CATEGORY ) ) // 
     &                     ' inventory file "' // LINE( 1:L2 ) // '"'
                    CNTR = CNTR + GETFLINE( TDEV, MESG )

                ELSEIF( FILFMT .EQ. EMSFMT ) THEN

C.....................  Make sure that next 4 files in list are also EMSFMT
C.....................  Increment line, scan for INVYEAR, open file, check file,
C                       write message, and store unit number.
                    IF( CATEGORY .EQ. 'POINT' ) THEN

                        J = J + 1
                        LINE = LSTSTR( J )
                        INFILE = LINE( 1:TRIMLEN( LINE ) )
                        IF( INDEX( INFILE,'INVYEAR' ) .GT. 0 ) GOTO 1007 ! Error

C......................... Close previous file, get unit, and open this file
                        CLOSE( TDEV )  
                        TDEV = JUNIT()
                        OPEN( TDEV, ERR=1006, FILE=INFILE, STATUS='OLD')
                        FILFMT = GETFORMT( TDEV )
                        IF( FILFMT .NE. EMSFMT ) GO TO 1008  ! Error
                        J = J + 3 ! Skip past other three files

                    ENDIF

C.....................  Sum number of lines in EMS-95 emission files
                    CNTR = CNTR + GETFLINE( TDEV, MESG )

                    CLOSE( TDEV )

                ELSE  ! File format not recognized	

                    MESG = 'File format is not recognized for file. ' //
     &                     CRLF() // BLANK10 // 
     &                     INFILE( 1:TRIMLEN( INFILE ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF

            ENDDO   ! End of loop through list-formatted PTINV file

            DEALLOCATE( LSTSTR )

        ELSEIF( INVFMT .EQ. EPSFMT .OR. INVFMT .EQ. IDAFMT ) THEN

            MESG = CATEGORY( 1:TRIMLEN( CATEGORY ) ) // 
     &             ' inventory file'
            CNTR = GETFLINE( TDEV, MESG )

            CLOSE( TDEV )

        ELSE
            WRITE( MESG,94010 ) 'INTERNAL ERROR: Illegal call to ' //
     &             'function ' // PROGNAME // ' with format ', INVFMT
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        ENDIF

        IF( CNTR .EQ. 0 ) THEN
            MESG = FMTNAMES( INVFMT ) // '-formatted inventory ' //
     &             'file has no valid lines of inventory data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        GETISIZE = CNTR

        RETURN

C******************  ERROR MESSAGES WITH EXIT **************************
 
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

        END
