
      PROGRAM INVSPLIT

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C       This program splits an IDA and toxics inventory file into 
C       multiple output files
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       PROMPTFFILE, ENVINT, STR2INT, FIND1, INDEX1
C
C  REVISION  HISTORY:
C       Created  12/02 by M Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: $Id$ 
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C*************************************************************************

C...........   MODULES for public variables
C.........  This module contains the information about the source category
        USE M3UTILIO

        USE MODINFO, ONLY: CATEGORY, CRL

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'     ! emissions constant parameters
C        INCLUDE 'PARMS3.EXT'      ! I/O API constants
C        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
C        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations

C...........   PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: SCCSW = '%W%'

C.........  EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)    CRLF
C       INTEGER         ENVINT
C       INTEGER         FIND1
        INTEGER         GETFLINE
        INTEGER         GETFORMT
C       INTEGER         INDEX1
C       INTEGER         PROMPTFFILE
C       INTEGER         STR2INT
 
C        EXTERNAL    CRLF, ENVINT, FIND1, GETFLINE, GETFORMT, INDEX1, 
C     &              PROMPTFFILE, STR2INT
        EXTERNAL     GETFLINE, GETFORMT

C...........   LOCAL PARAMETERS
        INTEGER, PARAMETER ::   NSEG = 80

C.........  Allocatable arrays...

C.........  Splits file arrays
        INTEGER, ALLOCATABLE :: SINDEX  ( : )      ! sorting index
        INTEGER, ALLOCATABLE :: STATLIST( : )      ! state codes
        INTEGER, ALLOCATABLE :: FNUMLIST( : )      ! file numbers
        INTEGER, ALLOCATABLE :: OUTNUM  ( : )      ! output file codes

C.........  Logical file names and unit numbers
        INTEGER, ALLOCATABLE :: ODEV( : )  ! output unit numbers
        INTEGER   IDEV         ! input inventory file
        INTEGER   LDEV         ! log file unit number
        INTEGER   SDEV         ! log file unit number

C.........  Local arrays
        CHARACTER(40) SEGMENT( NSEG )

C.........  Local variables

        INTEGER   I, CNT, F, J, PF         ! indices and counters

        INTEGER   IFMT         ! inventory format code
        INTEGER   IOS          ! i/o status
        INTEGER   IREC         ! record counter
        INTEGER   NLIST        ! number of entries in split file
        INTEGER   NOUT         ! number of output files
        INTEGER   STA          ! tmp state ID

        LOGICAL :: EFLAG = .FALSE.   ! true: error found

        CHARACTER(256)  MESG    ! temporary message array
        CHARACTER(2560) LINE    ! line buffer

        CHARACTER(IOVLEN3) INNAME   ! input file name
        CHARACTER(IOVLEN3) OUTNAME  ! output file names

        CHARACTER(16) :: PROGNAME = 'INVSPLIT'   !  program name

C***********************************************************************
C   begin body of program INVSPLIT

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Find name name of raw inventory file
        J = INDEX1( CATEGORY, NCAT, CATLIST )
        IF( J .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: Do not know about category ' //
     &             TRIM( CATEGORY ) // ' in program ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE 
            INNAME = ANAMLIST( J )

        END IF

C.........  Prompt for name of input inventory file
        MESG = 'Enter logical name of the RAW ' // 
     &         TRIM( CATEGORY ) // ' AVERAGE INVENTORY ' // 'file'

        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INNAME, PROGNAME )

C.........  Prompt for name of input splits definitions file
        MESG = 'Enter logical name of the SPLITS DEFINITIONS file' 
        SDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SPLITS', PROGNAME )

C.........   Get format of input file
        IFMT = GETFORMT( IDEV, -1 )

C.........   Get size for splits file
        NLIST = GETFLINE( SDEV, 'Splits definitions' )

C.........  Allocate memory for arrays
        ALLOCATE( STATLIST( NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STATLIST', PROGNAME )
        ALLOCATE( FNUMLIST( NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FNUMLIST', PROGNAME )
        ALLOCATE( SINDEX( NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SINDEX', PROGNAME )

C.........  Read split definitions and create list of output file
C           numbers
        IREC = 0
        CNT = 0
        DO I = 1, NLIST

            READ( SDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010)
     &              'I/O error', IOS, 'reading splits definitions '//
     &              'file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            IF( LINE( 1:1 ) .EQ. CINVHDR ) CYCLE

            CALL PARSLINE( LINE, NSEG, SEGMENT )

            CNT = CNT + 1
            IF( CNT .LE. NLIST ) THEN
                SINDEX  ( CNT ) = CNT
                STATLIST( CNT ) = STR2INT( SEGMENT( 1 ) )
                FNUMLIST( CNT ) = STR2INT( SEGMENT( 2 ) )
            END IF

        END DO
        NLIST = CNT

C.........  Exit if problem with splits input file
        IF( EFLAG ) THEN
            MESG = 'Problem reading splits definitions file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Create list of output file numbers and count the number
C           of outputs...
C.........  Sort the number of output files
        CALL SORTI1( NLIST, SINDEX, FNUMLIST )

C.........  Count number of outputs
        PF = IMISS3
        NOUT = 0
        DO I = 1, NLIST
            J = SINDEX( I )
            F = FNUMLIST( J )
            IF( F .NE. PF ) NOUT = NOUT + 1
            PF = F
        END DO

C.........  Allocate memory for output numbers
        ALLOCATE( OUTNUM( NOUT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTNUM', PROGNAME )
        ALLOCATE( ODEV( NOUT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ODEV', PROGNAME )

C.........  Store the output numbers
        PF = IMISS3
        NOUT = 0
        DO I = 1, NLIST
            J = SINDEX( I )
            F = FNUMLIST( J )
            IF( F .NE. PF ) THEN
                NOUT = NOUT + 1
                OUTNUM( NOUT ) = F
            END IF
            PF = F
        END DO

C.........  Prompt for output file names
        DO I = 1, NOUT
            WRITE( OUTNAME, '(A,I2.2)' ) 'OUTFILE', OUTNUM( I )
            WRITE( MESG, '(A,1X,I2.2)' ) 'Enter logical name ' //
     &             'for output file', OUTNUM( I ) 
            ODEV( I ) = PROMPTFFILE( MESG, .FALSE., .TRUE., OUTNAME, 
     &                               PROGNAME )
        END DO

C.........  Loop through input file and write output files
        DO

            READ( IDEV, 93000, END=199, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading inventory file at line', IREC
                CALL M3MESG( MESG )
                EXIT
            END IF

C.............  Skip header lines and write them to the output file
            IF( LINE( 1:1 ) .EQ. CINVHDR ) THEN
                DO I = 1, NOUT
                    WRITE( ODEV( I ), '(A)' ) TRIM( LINE )
                END DO
                CYCLE
            END IF

C.............  Convert state code to integer
            SELECT CASE ( IFMT )

            CASE ( ORLFMT, ORLNPFMT )
                CALL PARSLINE( LINE, NSEG, SEGMENT )
                STA = INT( STR2INT( SEGMENT(1) ) / 1000 )

            CASE ( FF10FMT, FF10DYFMT, FF10HRFMT )
                CALL PARSLINE( LINE, NSEG, SEGMENT )
                STA = INT( STR2INT( SEGMENT(2) ) / 1000 )

            CASE DEFAULT
                WRITE( MESG,94010 ) 'Cannot split file with ' // 
     &                              TRIM( FMTNAMES( IFMT ) ) // 'format'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END SELECT

C.............  Find state in list
            J = FIND1( STA, NLIST, STATLIST )

C.............  Error if state is not found
            IF( J .LE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: State', STA, 'is not ' //
     &                 'found in splits file.' // CRLF() // BLANK10 //
     &                 'Ensure that the state is in the list and  '// 
     &                 'the states are ' // CRLF() // BLANK10 //
     &                 'listed in sorted order.'
                CALL M3MSG2( MESG )
            END IF     

C.............  No point in writing output files if there is an error
            IF( EFLAG ) CYCLE

C.............  Check if output file number
            F = FNUMLIST( J )

C.............  Write line to correct output file
            WRITE( ODEV( F ), '(A)' ) TRIM( LINE )

        END DO

199     CONTINUE   ! exit from read loop

C.........  Check if error found
        IF( IOS .GT. 0 ) THEN
            MESG = 'Problem reading input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Normal completion of program
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )


      END PROGRAM INVSPLIT

