
        SUBROUTINE RDCODNAM( PDEV, VDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the master inventory pollutants list or activity list, concatenates
C     the names if they are longer than 16 characters, checks for duplicates,
C     renames if there are duplicates, and converts all names to uppercase. The 
C     subroutine returns the names in their original, unsorted order. If the
C     arrays are not empty on input, the subroutine appends the data from the
C     file to the input arrays.
C
C  PRECONDITIONS REQUIRED:
C     File unit FDEV already is opened
C     Memory allocated for CODES and NAMES
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: Models-3 subroutines
C     Functions: Models-3 functions
C
C  REVISION  HISTORY:
C     Created 10/98 by M. Houyoux
C
C****************************************************************************
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

C.........  MODULES for public variables
C.........  This module contains the lists of unique inventory information
        USE MODLISTS

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         GETFLINE
        INTEGER         STR2INT

        EXTERNAL CRLF, GETFLINE, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (IN) :: PDEV   ! input pollutants file unit no.
        INTEGER , INTENT (IN) :: VDEV   ! input activities file unit no.

C...........   Unsorted pollutant records
        INTEGER, ALLOCATABLE :: INDX1A( : )
        INTEGER, ALLOCATABLE :: INDX2A( : )
        INTEGER, ALLOCATABLE :: CODESA( : )
        INTEGER, ALLOCATABLE :: ISTATA( : )
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: NAMESA( : ) !  pollutant names
        CHARACTER(LEN=IOULEN3), ALLOCATABLE :: UNITSA( : ) !  units

C...........   Other local variables
        INTEGER         I, J, L, L1, L2     !  counters and indices

        INTEGER         COD     !  tmp for pollutant code
        INTEGER         IOS     !  i/o status
        INTEGER         LCOD    !  previous pollutant code
        INTEGER      :: MXADAT = 0  !  maximum number of valid activities
        INTEGER      :: MXPDAT = 0  !  maximum number of valid pollutants
        INTEGER         NDAT    !  number of data values
        INTEGER         STAT    !  tmp variable status (1=pol;-1=act)

        LOGICAL, SAVE :: FIRSTIME = .TRUE.

        CHARACTER*300   LNAM    !  previous pollutant name
        CHARACTER*300   PNAM    !  tmp for pollutant name (NOT correct length)
        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDCODNAM' ! program name

C***********************************************************************
C   begin body of subroutine RDCODNAM

C.........  Get file sizes and allocate memory...

C.........  Get no. lines in pollutant codes & activities files for allocating
C           memory     
        IF( PDEV .GT. 0 ) THEN
            MXPDAT = GETFLINE( PDEV, 'Pollutant codes and names file' )
        END IF
        IF( VDEV .GT. 0 ) THEN
            MXADAT = GETFLINE( VDEV, 'Activity names file' )
        END IF

        MXIDAT = MXPDAT + MXADAT

C.........  Allocate memory for storing contents of pollutants & activities
C           files, units, and local conversion factors
        ALLOCATE( INVDCOD( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDCOD', PROGNAME )
        ALLOCATE( INVDNAM( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDNAM', PROGNAME )
        ALLOCATE( INVSTAT( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSTAT', PROGNAME )
        ALLOCATE( INVDUNT( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDUNT', PROGNAME )
        ALLOCATE( INVDCNV( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDCNV', PROGNAME )

        INVDCOD = 0
        INVDNAM = ' '
        INVSTAT = 1             ! array (expected by PROCINVEN)
        INVDUNT = EMCMISS3
        INVDCNV = 1.            ! until otherwise found by header scan

C.........  Allocate local memory
        ALLOCATE( INDX1A( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDX1A', PROGNAME )
        ALLOCATE( INDX2A( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDX2A', PROGNAME )
        ALLOCATE( CODESA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CODESA', PROGNAME )
        ALLOCATE( ISTATA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISTATA', PROGNAME )
        ALLOCATE( NAMESA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NAMESA', PROGNAME )
        ALLOCATE( UNITSA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UNITSA', PROGNAME )

C.........  Read pollutants and activities files
        NDAT = 0
        CALL READ_POL_OR_ACT( PDEV,  1, NDAT )
        CALL READ_POL_OR_ACT( VDEV, -1, NDAT )

C.........  Check dimensions
        IF( NDAT .GT. MXIDAT ) THEN
            WRITE( MESG,94010 ) 
     &             'ERROR: Number of pollutant records :', NDAT, 
     &             CRLF() // BLANK10 // 'Memory allocated :', MXIDAT
            CALL M3MSG2( MESG )

            MESG = 'Insufficient memory allocated for ' //
     &             'codes/names file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE IF( NDAT .EQ. 0 ) THEN
            MESG ='No entries in any codes/names file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Sort by number for checking for duplicates
        CALL SORTI1( NDAT, INDX1A, CODESA )

C.........  Check for duplicate numbers - delete duplicates
        LCOD = -9
        DO I = 1, NDAT

            J = INDX1A( I )
            COD = CODESA( J )

            IF( COD .EQ. LCOD ) THEN
                WRITE( MESG,94010 )
     &                 'WARNING: Duplicate code "', COD,
     &                 '" at line', J, 'in codes/names file.' //
     &                 CRLF() // BLANK5 // 'Skipping record.'
                CALL M3MESG( MESG )

                CODESA( J ) = 0   ! flag for skipping

            END IF

            LCOD = COD

        END DO

C.........  Sort by name for checking for duplicates
        CALL SORTIC( NDAT, INDX2A, NAMESA )

C.........  Check for duplicate names - skip second duplicate
        LNAM = EMCMISS3
        DO I = 1, NDAT

            J = INDX2A( I )
            PNAM = NAMESA( J )

            IF( PNAM .EQ. LNAM .AND. CODESA( J ) .NE. 0 ) THEN

                L2 = LEN_TRIM( PNAM )
                WRITE( MESG,94010 )
     &                 'WARNING: Duplicate name "' // PNAM( 1:L2 ) //
     &                 '" at line', J, 'in codes/names file.' //
     &                 CRLF() // BLANK5 // 'Skipping record.'
                CALL M3MESG( MESG )

                CODESA( J ) = 0   ! flag for skipping

            END IF

            LNAM = PNAM

        END DO

C.........  Store valid entries in original (unsorted) order and convert to 
C           uppercase
C.........  Recompute maximum pollutants and activities
        J = 0
        DO I = 1, NDAT

            IF( CODESA( I ) .GT. 0 ) THEN
                J = J + 1

                STAT = ISTATA( I ) 

                INVDCOD( J ) = CODESA( I )
                INVDNAM( J ) = NAMESA( I )
                INVDUNT( J ) = UNITSA( I )
                INVSTAT( J ) = STAT

            END IF

        END DO

        MXIDAT = J

C.........  Rewind files
        REWIND( PDEV )
        REWIND( VDEV )

C.........  Deallocate local memory
        DEALLOCATE( INDX1A, INDX2A, CODESA, ISTATA, NAMESA )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram is for reading the pollutant or 
C               activity list with needed checks and possible error messages.
            SUBROUTINE READ_POL_OR_ACT( FDEV, STATVAL, CNT )

C.............  Subroutine arguments 
            INTEGER, INTENT (IN)     :: FDEV    ! unit number
            INTEGER, INTENT (IN)     :: STATVAL ! 1=pol; -1=act
            INTEGER, INTENT (IN OUT) :: CNT     ! ongoing count of all entries

C.............  Local arrays
            CHARACTER*50   SEGMENT( 3 )

C.............  Other local variables
            INTEGER         J       !  counters and indices
            INTEGER         IOS     !  i/o status
            INTEGER         IREC    !  record counter   
            INTEGER         L, L2   !  length indices

            CHARACTER*300   PNAM    !  tmp for pollutant name (NOT correct length)
            CHARACTER*300   UNIT    !  tmp for units (NOT correct length)
            CHARACTER*300   LINE    !  line buffer

C----------------------------------------------------------------------

C.............  Do not read if file is not opened
            IF( FDEV .LE. 0 ) RETURN

C.............  Loop through input file...
            IREC = 0
            DO

                READ( FDEV, 93000, END=22, IOSTAT=IOS ) LINE
                IREC = IREC + 1

                IF( IOS .GT. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                 'Error', IOS,  'reading names ' // 
     &                 'and codes file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Parse the line into its sections
                CALL PARSLINE( LINE, 3, SEGMENT )

                COD = STR2INT( SEGMENT( 1 ) )

                PNAM = SEGMENT( 2 )
                L1 = LEN_TRIM( PNAM )

                UNIT = SEGMENT( 3 )
                L2 = LEN_TRIM( UNIT )

C.................  Check code number
                IF( COD .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: ' //
     &                 'Blank, alphabetic, or 0 data code at line',
     &                 IREC, CRLF() // BLANK10 //
     &                 'in code/names file. Skipping record.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Truncate name to IOVLEN3 characters
                IF( L1 .GT. IOVLEN3 ) THEN

                    PNAM = PNAM( 1:IOVLEN3 )

                    WRITE( MESG,94010 )
     &                 'WARNING: Name too long for ' //
     &                 'I/O API variable names at line', IREC,
     &                 CRLF() // BLANK5 //
     &                 'in pollutant file. Truncating to "' //
     &                 PNAM( 1:IOVLEN3 ) // '"'
                    CALL M3MESG( MESG )

                END IF

C.................  Truncate units to IOULEN3 characters
                IF( L2 .GT. IOULEN3 ) THEN

                    UNIT = UNIT( 1:IOULEN3 )

                    WRITE( MESG,94010 )
     &                 'WARNING: Units too long for ' //
     &                 'I/O API units at line', IREC,
     &                 CRLF() // BLANK5 //
     &                 'in pollutant file. Truncating to "' //
     &                 UNIT( 1:IOULEN3 ) // '"'
                    CALL M3MESG( MESG )

                END IF

C.................  Store unsorted variables
                CNT = CNT + 1
                IF( CNT .LE. MXIDAT ) THEN
                    INDX1A( CNT ) = CNT
                    INDX2A( CNT ) = CNT
                    CODESA( CNT ) = COD
                    ISTATA( CNT ) = STATVAL
                    NAMESA( CNT ) = PNAM( 1:IOVLEN3 )
                    UNITSA( CNT ) = UNIT( 1:IOULEN3 )
                END IF

            END DO          !  end read loop

22          CONTINUE        !  exit from loop reading FDEV

            RETURN

C-------------------  FORMAT  STATEMENTS   ----------------------------

C...............   Formatted file I/O formats............ 93xxx

93000       FORMAT( A )

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE READ_POL_OR_ACT

        END SUBROUTINE RDCODNAM
