
        SUBROUTINE RDEPROC( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     This subroutine reads the emission processes file, which contains columns
C     for the activity, associated process, and associated pollutants.  If
C     there is more than one process per activity, then these are listed on
C     separate lines in the file. If there is more than one pollutant per
C     activity and process, then these are listed in additional columns.
C     During the read, the column number is set dynamically.  Only processes
C     for activities that are in the inventory are read.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 10/99 by M. Houyoux
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

C.........  MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     GETFLINE
        INTEGER     GETNLIST
        INTEGER     INDEX1

        EXTERNAL    GETFLINE, GETNLIST, INDEX1

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! cross-reference file unit no.

C...........   Local allocatable arrays
        INTEGER               , ALLOCATABLE :: INDX   ( : )  ! POLA sorting indx
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: SEGMENT( : )  ! line segments
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: POLNAM ( : )  ! pol names
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: POLA   ( : )  ! all unsorted pols

C...........   Other local variables
        INTEGER         I, J, K, L, L2, V    !  counters and indices

        INTEGER         IOS     !  i/o status
        INTEGER         L1, L2  !  tmp string lengths
        INTEGER         LJ      !  string length for emis type joiner
        INTEGER         MXCOLS  !  maximum number of columns in the file
        INTEGER         MXPOL   !  maximum number of pollutants per process
        INTEGER         NCOLS   !  no. columns in a row
        INTEGER         NLINES  !  number of lines in file
        INTEGER         NPOL    !  no. pollutants in a row
        INTEGER         NPUNSRT !  no. pols 

        LOGICAL      :: EFLAG = .FALSE.   !  true: error found

        CHARACTER*300          LINE     !  line buffer
        CHARACTER*300          MESG     !  message buffer
        CHARACTER(LEN=IOVLEN3) ACT      !  tmp activity name
        CHARACTER(LEN=IOVLEN3) CPOL     !  tmp pollutant name
        CHARACTER(LEN=IOVLEN3) LPOL     !  tmp pollutant name from previous iter
        CHARACTER(LEN=IOVLEN3) PRC      !  tmp process name

        CHARACTER*16 :: PROGNAME = 'RDEPROC' ! program name

C***********************************************************************
C   begin body of subroutine RDEPROC

C.........  Get the number of lines for the file and allocate array so that
C           the type of the line can be stored
        NLINES = GETFLINE( FDEV, 'Emission processes file' )

C.........  Read through file to determine the maximum no. of columns
        MXCOLS = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading emission processes file at line', I
                CALL M3MESG( MESG )
                CYCLE
            END IF

            L1 = LEN_TRIM( LINE )
            J  = GETNLIST( L1, LINE )

            IF( J .GT. MXCOLS ) MXCOLS = J

        END DO

        MXPOL = MXCOLS - 2

C.........  Allocate memory for parsing line segements and storing pollutants
        ALLOCATE( SEGMENT( MXCOLS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )
        ALLOCATE( POLNAM( MXPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLNAM', PROGNAME )

C.........  Allocate memory for emission types and count for each per activity
        ALLOCATE( EMTNAM( NLINES*MXPOL, NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTNAM', PROGNAME )
        ALLOCATE( EMTIDX( NLINES*MXPOL, NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTIDX', PROGNAME )
        ALLOCATE( NETYPE( NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NETYPE', PROGNAME )

        EMTNAM = ' '  ! array
        EMTIDX = 0    ! array
        NETYPE = 0    ! array

C.........  Rewind file
        REWIND( FDEV )

C.........  Store contents of emissions processes file in output order
        J = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

C.............  Skip blank line
            IF( LINE .EQ. ' ' ) CYCLE

C.............  Separate line into segments
            L1 = LEN_TRIM( LINE )
            NCOLS = GETNLIST( L1, LINE )
            CALL PARSLINE( LINE, NCOLS, SEGMENT )

            ACT = SEGMENT( 1 )
            PRC = SEGMENT( 2 )

            NPOL = NCOLS - 2
            POLNAM( 1:NPOL ) = SEGMENT( 3:NCOLS )

C.............  Make sure activity is in the inventory
            K = INDEX1( ACT, NIACT, ACTVTY )

C.............  Store emission processes and associated pollutants
            IF( K .GT. 0 ) THEN

                DO V = 1, NPOL
                    J = J + 1
                    L1 = LEN_TRIM( PRC )
                    L2 = LEN_TRIM( POLNAM( V ) )
                    EMTNAM( J,K ) = PRC( 1:L1 ) // ETJOIN // 
     &                              POLNAM( V )( 1:L2 )
                END DO

                NETYPE( K ) = NETYPE( K ) + NPOL

            END IF
    
        END DO      ! End loop over file

C.........  Set the maximum number of emission types
        MXETYPE = MAXVAL( NETYPE )

C.........  Create a list of pollutants associated with the emission types...

C.........  Allocate memory for unsorted pollutants list
        ALLOCATE( POLA( MXETYPE * NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLA', PROGNAME )
        ALLOCATE( INDX( MXETYPE * NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDX', PROGNAME )

C.........  Create unsorted pollutants list
        J = 0
        LJ = LEN_TRIM( ETJOIN )
        DO I = 1, NIACT

            DO K = 1, NETYPE( I )
                J = J + 1

                L = INDEX( EMTNAM( K,I ), ETJOIN )
                L2 = LEN_TRIM( EMTNAM( K,I ) )

                POLA( J ) = EMTNAM( K,I )( L+LJ:L2 )
                INDX( J ) = J

            END DO

        END DO
        NPUNSRT = J

C.........  Sort pollutants
        CALL SORTIC( NPUNSRT, INDX, POLA )

C.........  Determine number of actual pollutants
        LPOL = '-9'
        K = 0
        DO I = 1, NPUNSRT

            J = INDX( I )
            CPOL = POLA( J )

            IF( CPOL .NE. LPOL ) THEN
                K = K + 1
                LPOL = CPOL
            END IF

        END DO

        NEPOL = K

C.........  Allocate memory for sorted pollutants
        ALLOCATE( EMTPOL( NEPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTPOL', PROGNAME )

C.........  Store sorted pollutants in a unique list
        LPOL = '-9'
        K = 0
        DO I = 1, NPUNSRT

            J = INDX( I )
            CPOL = POLA( J )

            IF( CPOL .NE. LPOL ) THEN
                K = K + 1
                EMTPOL( K ) = CPOL
                LPOL = CPOL
            END IF

        END DO

C.........  Create an index that references between the emission types and the
C           list of unique pollutants
        DO I = 1, NIACT

            DO K = 1, NETYPE( I )

                L = INDEX( EMTNAM( K,I ), ETJOIN )
                L2 = LEN_TRIM( EMTNAM( K,I ) )

                CPOL = EMTNAM( K,I )( L+LJ:L2 )
                J = INDEX1( CPOL, NEPOL, EMTPOL )
                EMTIDX( K,I ) = J

            END DO

        END DO

C.........  Rewind file
        REWIND( FDEV )

C.........  Deallocate memory for local arrays
        DEALLOCATE( SEGMENT, POLNAM, POLA, INDX )

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of emission processes file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDEPROC
