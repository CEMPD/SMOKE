
        PROGRAM EXTRACTIDA

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION: 
C      Program to read a list of FIPS and PLANTS, make sure that the 
C      list is unique, then output only those FIPS and plans to a 
C      new IDA file
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Written memorial day weekend (argh), 5/2002 by M. Houyoux
C
C***********************************************************************
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
C*************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(2)    CRLF
        INTEGER         FINDC
        INTEGER         GETFLINE
        INTEGER         PROMPTFFILE

        EXTERNAL   CRLF, FINDC, GETFLINE, PROMPTFFILE

C...........   LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.8.1_Jan2021' ! CVS release tag

C............   Allocatable arrays
        INTEGER,               ALLOCATABLE :: INDEXA  ( : ) ! sorting index
        CHARACTER(PTENDL3(2)), ALLOCATABLE :: CFPLISTA( : ) ! unsorted dups
        CHARACTER(PTENDL3(2)), ALLOCATABLE :: CFPLIST ( : ) ! sorted no dups

C...........   Logical file names and unit numbers

        INTEGER         LDEV    !  log-device
        INTEGER         IDEV    !  for  input PTINV file
        INTEGER         ODEV    !  for output PTINV file
        INTEGER         XDEV    !  for extraction list

C...........   Other local variables
        INTEGER          I, J, L, N             !  counters, subscripts
        INTEGER          IOS                    !  tmp i/o status
        INTEGER          IREC                   !  tmp record number
        INTEGER          NRECS                  !  number of records in XDEV
        INTEGER          NFPLIST                !  no. entries in CFPLIST
        INTEGER          NLINES                 !  number of lines in IDEV

        LOGICAL       :: EFLAG = .FALSE.        !  true: error found

        CHARACTER(FIPLEN3)    CFIP          !  just FIPS code
        CHARACTER(PTENDL3(2)) CFIPPLANT     !  FIPS // PLANT
        CHARACTER(PTENDL3(2)) PREVFP        !  previous iteration CFIPPLANT
        CHARACTER(PLTLEN3)    PLANT         !  just plant ID
        CHARACTER(256)        MESG          !  message buffer
        CHARACTER(512)        LINE          !  tmp PTINV input line

        CHARACTER(16) :: PROGNAME = 'EXTRACTIDA'  !  program name

C***********************************************************************
C   begin body of program EXTRACTIDA

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Open input FIPS/plant list file
        XDEV = PROMPTFFILE( 
     &         'Enter logical name for input FIPS/plants list',
     &         .TRUE., .TRUE., 'INRECS', PROGNAME )

        IDEV = PROMPTFFILE( 
     &         'Enter logical name for input PTINV file',
     &         .TRUE., .TRUE., 'INFILE', PROGNAME )
               
        ODEV = PROMPTFFILE( 
     &         'Enter logical name for output PTINV file',
     &         .FALSE., .TRUE., 'OUTFILE', PROGNAME )

        NRECS  = GETFLINE( XDEV, 'FIPS/plants list' )
        NLINES = GETFLINE( IDEV, 'Input PTINV file' )

C........  Allocat memory for FIP/plant lists
        ALLOCATE( CFPLISTA( NRECS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFPLISTA', PROGNAME )
        ALLOCATE( CFPLIST( NRECS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFPLIST', PROGNAME )
        ALLOCATE( INDEXA( NRECS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )

C........  Read FIP/plant lists
        IREC = 0
        I = 0
        DO N = 1, NRECS

            READ( XDEV, 93000, END = 9991, IOSTAT = IOS ) CFIPPLANT
            IREC = IREC + 1
                
            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'I/O error', IOS,
     &              'reading FIPS/plant list file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C...........  Check for header lines and skip them
            IF( CFIPPLANT( 1:1 ) .EQ. CINVHDR ) CYCLE

            CFIP  = ADJUSTR( CFIPPLANT( PTBEGL3(1):PTENDL3(1) ) )
            PLANT = ADJUSTL( CFIPPLANT( PTBEGL3(2):PTENDL3(2) ) )

            I = I + 1
            CFPLISTA( I ) = CFIP // PLANT
            INDEXA  ( I ) = N

        END DO
        NRECS = I

C........  Sort FIP/plant list
        CALL SORTIC( NRECS, INDEXA, CFPLISTA )

C........  Restore FIP/plant list without duplicates
        I = 0
        PREVFP = ' '
        DO N = 1, NRECS

            J = INDEXA( N )
            IF( CFPLISTA( J ) .NE. PREVFP ) THEN
                I = I + 1
                CFPLIST( I ) = CFPLISTA( J )
            END IF

            PREVFP = CFPLISTA( J )

        END DO
        NFPLIST = I

C........  Read input file and output records that are in FIP/plant list
        IREC = 0
        DO N = 1, NLINES

            READ( IDEV, 93000, END = 9992, IOSTAT = IOS ) LINE
            IREC = IREC + 1
                
            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'I/O error', IOS,
     &              'reading input PTINV file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C...........  CHeck for header line and write it out
            IF( LINE( 1:1 ) .EQ. CINVHDR ) THEN
                L = LEN_TRIM( LINE )
                WRITE( ODEV, '(A)' ) LINE( 1:L )
                CYCLE
            END IF

C...........  Search for FIP / plant combo in input file 
            CFIP  = ADJUSTR( LINE( PTBEGL3(1):PTENDL3(1) ) )
            PLANT = ADJUSTL( LINE( PTBEGL3(2):PTENDL3(2) ) )

            CFIPPLANT = CFIP // PLANT

            I = FINDC( CFIPPLANT, NFPLIST, CFPLIST )

C...........  Ouptput line if record matches input list 
            IF( I .GT. 0 ) THEN
                L = LEN_TRIM( LINE )
                WRITE( ODEV, '(A)' ) LINE( 1:L )
            END IF

        END DO

        IF( EFLAG ) THEN
            MESG = 'Problem reading one or more input files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C..........   Normal completion
        MESG = ' ' 
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 0 )

9991    WRITE( MESG,94010 ) 'End of FIPS/plant list reached ' // 
     &         'unexpectedly at line', IREC
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

9992    WRITE( MESG,94010 ) 'End of PTINV input file reached ' // 
     &         'unexpectedly at line', IREC
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )

        END PROGRAM EXTRACTIDA
