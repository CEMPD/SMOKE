
        SUBROUTINE RDEPROC( FDEV )

C***********************************************************************
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
C     Modified  3/14 by B.H. Baek
C
C****************************************************************************/
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

C.........  MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC, ONLY: EMTNAM, EMTIDX, NETYPE, MXETYPE, NEPOL, EMTPOL

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIACT, ACTVTY

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     GETFLINE
        INTEGER     GETNLIST
        INTEGER     INDEX1
        LOGICAL     BLKORCMT

        EXTERNAL    BLKORCMT, GETFLINE, GETNLIST, INDEX1

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! cross-reference file unit no.

C...........   Local allocatable arrays
        INTEGER           , ALLOCATABLE :: INDX   ( : )  ! POLA sorting indx
        CHARACTER(IOVLEN3), ALLOCATABLE :: SEGMENT( : )  ! line segments
        CHARACTER(IOVLEN3), ALLOCATABLE :: POLNAM ( : )  ! pol names
        CHARACTER(IOVLEN3), ALLOCATABLE :: POLA   ( : )  ! all unsorted pols

C...........   Local parametes
        CHARACTER, PARAMETER :: CONTCHAR = '\'  ! line continuation character

C...........   Other local variables
        INTEGER         I, J, K, L, M, V    !  counters and indices

        INTEGER         IOS     !  i/o status
        INTEGER         L1, L2  !  tmp string lengths
        INTEGER         LJ      !  string length for emis type joiner
        INTEGER         MXCOLS  !  maximum number of columns in the file
        INTEGER         MXPOL   !  maximum number of pollutants per process
        INTEGER         NCOLS   !  no. columns in a row
        INTEGER         NLINES  !  number of lines in file
        INTEGER         NPOL    !  no. pollutants in a row
        INTEGER         NPUNSRT !  no. pols 
        INTEGER         NTLINES !  no. of lines taking into account continuation lines 

        LOGICAL      :: FNDPOL   = .FALSE.   ! true: found pollutant on master list
        LOGICAL      :: NEWLINE  = .TRUE.    ! true: current line is new (not continued)
        LOGICAL      :: EFLAG    = .FALSE.   ! true: error found

        CHARACTER(3000)    LINE     !  line buffer
        CHARACTER(300)     MESG     !  message buffer
        CHARACTER(IOVLEN3) ACT      !  tmp activity name
        CHARACTER(IOVLEN3) CPOL     !  tmp pollutant name
        CHARACTER(IOVLEN3) LPOL     !  tmp pollutant name from previous iter

        CHARACTER(16) :: PROGNAME = 'RDEPROC' ! program name

C***********************************************************************
C   begin body of subroutine RDEPROC

C.........  Get the number of lines for the file and allocate array so that
C           the type of the line can be stored
        NLINES = GETFLINE( FDEV, 'Emission processes file' )
        NTLINES = NLINES

C.........  Read through file to determine the maximum no. of columns and
C           check packet information
        MXCOLS   = 0
        
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

C.............  Skip any blank and comment lines
            IF( BLKORCMT( LINE ) ) THEN
                NTLINES = NTLINES - 1
                CYCLE
            END IF

            LINE = ADJUSTL( LINE )
            L1 = LEN_TRIM( LINE )
            
C.............  If previous line was continued, add this lines columns to number
C               from previous line; otherwise, count columns in just this line
            IF( .NOT. NEWLINE ) THEN
                NCOLS = NCOLS + GETNLIST( L1, LINE )
            ELSE
                NCOLS = GETNLIST( L1, LINE )
            END IF

C.............  Check for continuation character in current line
            IF( LINE( L1:L1 ) == CONTCHAR ) THEN
                NEWLINE = .FALSE.
                NTLINES = NTLINES - 1
                NCOLS = NCOLS - 1
            ELSE
                NEWLINE = .TRUE.
            END IF
            
            IF( NCOLS > MXCOLS ) THEN
                MXCOLS = NCOLS
            END IF

        END DO

C.........  Check for errors so far
        IF( EFLAG ) THEN
            MESG = 'Problem reading emission processes file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MXPOL = MXCOLS - 1

C.........  Allocate memory for parsing line segements and storing pollutants
        ALLOCATE( SEGMENT( MXCOLS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )
        ALLOCATE( POLNAM( MXPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLNAM', PROGNAME )

C.........  Allocate memory for emission types and count for each per activity
        ALLOCATE( EMTNAM( NTLINES*MXPOL, NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTNAM', PROGNAME )
        ALLOCATE( EMTIDX( NTLINES*MXPOL, NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTIDX', PROGNAME )
        ALLOCATE( NETYPE( NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NETYPE', PROGNAME )

        EMTNAM = ' '  ! array
        EMTIDX = 0    ! array
        NETYPE = 0    ! array

C.........  Rewind file
        REWIND( FDEV )

        NEWLINE  = .TRUE.

C.........  Store contents of emissions processes file in output order
        J = 0
        L = 0
        M = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

            LINE = ADJUSTL( LINE )
            L1 = LEN_TRIM( LINE )
           
C.............  Separate line into segments
            NCOLS = GETNLIST( L1, LINE )
            CALL PARSLINE( LINE, NCOLS, SEGMENT )

            IF( NEWLINE ) THEN
                ACT = SEGMENT( 1 )

                IF( SEGMENT( NCOLS ) == CONTCHAR ) THEN
                    NPOL = NCOLS - 2
                    POLNAM( 1:NPOL ) = SEGMENT( 2:NCOLS-1 )
                ELSE
                    NPOL = NCOLS - 1
                    POLNAM( 1:NPOL ) = SEGMENT( 2:NCOLS )
                END IF
            ELSE
                IF( SEGMENT( NCOLS ) == CONTCHAR ) THEN
                    NPOL = NCOLS - 1
                    POLNAM( 1:NPOL ) = SEGMENT( 1:NCOLS-1 )
                ELSE
                    NPOL = NCOLS
                    POLNAM( 1:NPOL ) = SEGMENT( 1:NCOLS )
                END IF
            END IF

C.............  Make sure activity is in the inventory
            K = INDEX1( ACT, NIACT, ACTVTY )

C.............  Store emission processes and associated pollutants
            IF( K .GT. 0 ) THEN

                DO V = 1, NPOL
                    J = J + 1
                    EMTNAM( J,K ) = POLNAM( V )
                END DO

                NETYPE( K ) = NETYPE( K ) + NPOL

            END IF

C.............  Check if current line is continued
            IF( SEGMENT( NCOLS ) == CONTCHAR ) THEN
                NEWLINE = .FALSE.
            ELSE
                NEWLINE = .TRUE.
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
        DO I = 1, NIACT

            DO K = 1, NETYPE( I )
                J = J + 1
                POLA( J ) = EMTNAM( K,I )
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

                CPOL = EMTNAM( K,I )
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
