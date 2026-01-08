
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
C      09/2025 by HT UNC-IE:  Use M3UTILIO
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
        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
c       INCLUDE 'PARMS3.EXT'    !  I/O API parameters
c       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
c       INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

c       CHARACTER(2) CRLF
c       INTEGER      GETFORMT
c       INTEGER      GETINVYR
c       INTEGER      JUNIT

c       EXTERNAL     CRLF, GETFORMT, GETINVYR, JUNIT
        INTEGER, EXTERNAL :: GETFORMT
        INTEGER, EXTERNAL :: GETINVYR

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NLINE            ! number of lines in file
        CHARACTER(*), INTENT (IN) :: FNAME            ! logical name of file
        CHARACTER(*), INTENT (IN) :: NLSTSTR( NLINE ) ! contents of file by line
        INTEGER     , INTENT(OUT) :: FILFMT ( NLINE ) ! file format code

C...........   File units and logical/physical names
        INTEGER      TDEV        !  emissions file in list format file

C...........   Other local variables
        INTEGER      I, J
        INTEGER      EXTFORMAT           !  used when a #LIST entry has format info
        INTEGER      IOS                 !  I/O status

        LOGICAL   :: EMSFLAG   = .FALSE. !  true: at least one file is EMS format
        LOGICAL   :: ORLFLAG  = .FALSE.  !  true: at least one file is ORL format

        CHARACTER(300)  INFILE      !  input file line buffer
        CHARACTER(500)  MESG        !  message buffer

        CHARACTER(16) :: PROGNAME =  'CHKLSTFL' ! program name

C***********************************************************************
C   begin body of subroutine CHKLSTFL

        EMSFLAG  = .FALSE.   ! Need to reset for each each subroutine call
        ORLFLAG  = .FALSE.
        EXTFORMAT = -1
        
C.........  Loop through lines of list-formatted file to check the formats
        DO J = 1, NLINE

C.............  Skip blank lines
            IF( NLSTSTR( J ) == ' ' ) CYCLE

C.............  Store the current line's file name  
            INFILE = NLSTSTR( J )

C.............  Skip INVYEAR packet 
            I = GETINVYR( INFILE )
            IF( I .GT. 0 ) THEN
                FILFMT( J ) = -1
                CYCLE
            END IF

C.............  Skip the date range packet
            I = INDEX( INFILE, 'DATERANGE' )
            IF( I .GT. 0 ) THEN
                FILFMT( J ) = -1
                CYCLE
            END IF

C.............  Check for #LIST entry
            I = INDEX( INFILE, '#LIST' )
            IF( I .GT. 0 ) THEN
                FILFMT( J ) = -1
                    
                IF( INDEX( INFILE, 'EMS-95' ) > 0 ) THEN
                    EXTFORMAT = EMSFMT
                    
                ELSE IF( INDEX( INFILE, 'CEM' ) > 0 ) THEN
                    EXTFORMAT = CEMFMT
                    
                ELSE IF( INDEX( INFILE, 'FF10' ) > 0 ) THEN
                    EXTFORMAT = FF10FMT 

                ELSE IF( INDEX( INFILE, 'MEDS' ) > 0 ) THEN
                    EXTFORMAT = MEDSFMT

                ELSE IF( INDEX( INFILE, 'GRID' ) > 0 ) THEN
                    EXTFORMAT = NCDFMT

                ELSE IF( INDEX( INFILE, 'ORL' ) > 0 ) THEN

                    IF( INDEX( INFILE, 'NONPOINT' ) > 0 ) THEN
                        EXTFORMAT = ORLNPFMT

                    ELSE IF( INDEX( INFILE, 'FIRE' ) > 0 ) THEN
                        EXTFORMAT = ORLFIREFMT

                    ELSE IF( INDEX( INFILE, 'FIREEMIS' ) > 0 ) THEN
                        EXTFORMAT = ORLDYFRFMT

                    ELSE
                        EXTFORMAT = ORLFMT
                        WRITE( MESG,94010 ) 'WARNING: No format '//
     &                    'included with #LIST ORL entry in '//
     &                    'inventory input file.'// CRLF() // BLANK5 // 
     &                    'Assuming ORL point, nonroad, or on-road.'

                        CALL M3MSG2( MESG )

                    END IF

                END IF
                
                CYCLE
            END IF

C.............  Open INFILE
            TDEV = JUNIT()

            OPEN( TDEV, FILE=INFILE, STATUS='OLD', IOSTAT=IOS )

C.............  Check for problems opening raw input file
            IF( IOS /= 0 ) THEN
                WRITE( MESG,94010 ) 'Problem at line ', J, 'of ' //
     &             TRIM( FNAME ) // '.' // ' Could not open file:' //
     &             CRLF() // BLANK5 // TRIM( INFILE )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        
C.............  Determine format of INFILE
            FILFMT( J ) = GETFORMT( TDEV, EXTFORMAT )
            CLOSE( TDEV )

C.............  Set flag based on format
            IF( FILFMT( J ) == EMSFMT ) EMSFLAG = .TRUE.
            IF( FILFMT( J ) == FF10FMT    .OR.
     &          FILFMT( J ) == FF10DYFMT  .OR.
     &          FILFMT( J ) == FF10HRFMT  .OR.
     &          FILFMT( J ) == MEDSFMT    .OR.
     &          FILFMT( J ) == ORLFMT     .OR.
     &          FILFMT( J ) == ORLNPFMT   .OR.
     &          FILFMT( J ) == ORLFIREFMT .OR.
     &          FILFMT( J ) == ORLDYFRFMT     ) ORLFLAG = .TRUE.

C.............  Check that file formats are consistent
            IF( EMSFLAG .AND. FILFMT( J ) /= EMSFMT ) THEN
                WRITE( MESG,94010 ) 
     &                 'ERROR: In SMOKE list-formatted inventory file, '
     &                 // TRIM( FNAME ) // ', at least one file is ' //
     &                 CRLF() // BLANK10 // 'EMS-95 format ' //
     &                 'while another is not. When using an EMS-95 ' //
     &                 'formatted' // CRLF() // BLANK10 // 
     &                 'inventory, all other inventory files must ' //
     &                 'also be in EMS-95 format.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( ORLFLAG                   .AND. 
     &          FILFMT( J ) /= FF10FMT    .AND.
     &          FILFMT( J ) /= FF10DYFMT  .AND.
     &          FILFMT( J ) /= FF10HRFMT  .AND.
     &          FILFMT( J ) /= MEDSFMT    .AND.
     &          FILFMT( J ) /= ORLFMT     .AND.
     &          FILFMT( J ) /= ORLNPFMT   .AND.
     &          FILFMT( J ) /= ORLFIREFMT .AND.
     &          FILFMT( J ) /= ORLDYFRFMT      ) THEN
                WRITE( MESG,94010 )
     &                 'ERROR: In SMOKE list-formatted inventory file, '
     &                 // TRIM( FNAME ) // ', at least one file is ' //
     &                 CRLF() // BLANK10 // 'ORL ' //
     &                 'format while another is not ' //
     &                 'ORL format.' // CRLF() // BLANK10 //
     &                 'When using ORL formatted ' //
     &                 'inventories, all other inventories ' // 
     &                 CRLF() // BLANK10 // 'must also be ORL '
     &                 // 'format.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            EXTFORMAT = FILFMT( J )

        END DO     ! End of loop through list-formatted file

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE CHKLSTFL

