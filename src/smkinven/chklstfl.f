
        SUBROUTINE CHKLSTFL( NLINE, FNAME, NLSTSTR, FILFMT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine checks the format of the list-formatted inventory
C      file and returns the code for the type of files it contains.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 10/98 by M. Houyoux
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

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2  CRLF
        INTEGER      GETFORMT
        INTEGER      GETINVYR
        INTEGER      JUNIT

        EXTERNAL     CRLF, GETFORMT, GETINVYR, JUNIT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NLINE            ! number of lines in file
        CHARACTER(*), INTENT (IN) :: FNAME            ! logical name of file
        CHARACTER(*), INTENT (IN) :: NLSTSTR( NLINE ) ! contents of file by line
        INTEGER     , INTENT(OUT) :: FILFMT           ! file format code

C...........   File units and logical/physical names
        INTEGER      TDEV        !  emissions file in list format file

C...........   Other local variables
        INTEGER      I, J

        INTEGER      FLEN        !  length of string FNAME
        INTEGER      PREVFMT     !  file format code of previous iteration

        LOGICAL   :: EFLAG     = .FALSE. !  true: error found
        LOGICAL   :: FIRSTITER = .TRUE.  !  true: first iteration of a loop

        CHARACTER*300   INFILE      !  input file line buffer
        CHARACTER*300   MESG        !  message buffer

        CHARACTER*16 :: PROGNAME =  'CHKLSTFL' ! program name

C***********************************************************************
C   begin body of subroutine CHKLSTFL

        FLEN = LEN_TRIM( FNAME )

        FIRSTITER = .TRUE.
C.........  Loop through lines of list-formatted file to check the formats
        DO J = 1, NLINE

C.............  Store the current line's file name  
            INFILE = NLSTSTR( J )

C.............  Skip INVYEAR packet 
            I = GETINVYR( INFILE )
            IF( I .GT. 0 ) CYCLE

C.............  Skip the date range packet
            I = INDEX( INFILE, 'DATERANGE' )
            IF( I .GT. 0 ) CYCLE

C.............  Open INFILE
            TDEV = JUNIT()
            OPEN( TDEV, ERR=1006, FILE=INFILE, STATUS='OLD' )

C.............  Determine format of INFILE
            FILFMT = GETFORMT( TDEV )

C.............  Make sure that file format was found
            IF( FILFMT .LT. 0 ) THEN
                
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 
     &                 'ERROR: In SMOKE list-formatted inventory file, '
     &                 // FNAME( 1:LEN_TRIM( FNAME ) )// ', could '//
     &                 CRLF() // BLANK10 // 
     &                 'not determine format of file listed at line', J
                CALL M3MESG( MESG )

C.............  If first iteration, save format, if not, make sure 
C               that different formats are not used in same PTINV list
            ELSE IF( FIRSTITER ) THEN
                FIRSTITER = .FALSE.
                PREVFMT = FILFMT

            ELSEIF( FILFMT .NE. PREVFMT ) THEN
                
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 
     &                 'ERROR: In SMOKE list-formatted inventory file, '
     &                 // FNAME( 1:LEN_TRIM( FNAME ) )// ', previous '//
     &                 CRLF() // BLANK10 // 
     &                 'file was ' // FMTNAMES( PREVFMT ) // ' format, '
     &                 // 'but file at line', J, 'is ' // 
     &                 FMTNAMES( FILFMT ) // ' format.'
                CALL M3MESG( MESG )

            ENDIF

            CLOSE( TDEV )

        END DO     ! End of loop through list-formatted PTINV file

C.........  Exit if files in list-formatted file were of inconsistent type
        IF( EFLAG ) THEN
            MESG = 'Problem reading SMOKE list format'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF
 
        RETURN

C******************  ERROR MESSAGES WITH EXIT **************************
 
C.........  Error opening raw input file
1006    WRITE( MESG,94010 ) 'Problem at line ', J, 'of ' //
     &         FNAME( 1:FLEN ) // '.' // ' Could not open file:' //
     &         CRLF() // BLANK5 // INFILE( 1:LEN_TRIM( INFILE ) )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE CHKLSTFL

