
        SUBROUTINE RDMVINFO( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 5/99 by M. Houyoux
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
C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         GETFLINE
        INTEGER         STR2INT

        EXTERNAL  CRLF, GETFLINE, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! cross-reference file unit no.
 
C...........   Local arrays for reading mobile information file
        INTEGER, ALLOCATABLE :: LTYPE   ( : ) !  type of line infile
        INTEGER, ALLOCATABLE :: RDCLSA  ( : ) !  unsorted AIRS AMS road classes
        INTEGER, ALLOCATABLE :: RDWAYA  ( : ) !  unsorted roadway types
        INTEGER, ALLOCATABLE :: RIDXA   ( : ) !  index for sorting road classes
        INTEGER, ALLOCATABLE :: IVTIDA  ( : ) !  vehicle type code
        INTEGER, ALLOCATABLE :: VIDXA   ( : ) !  index for sorting vehicle types

        CHARACTER(LEN=VTPLEN3), ALLOCATABLE :: CVTYPA( : ) ! unsorted veh types

C...........   Other local variables
        INTEGER         I, J, K1, K2, N1, N2    !  counters and indices

        INTEGER         IOS     !  i/o status
        INTEGER         NLINES  !  number of lines

        LOGICAL      :: EFLAG = .FALSE.   !  true: error found
        LOGICAL      :: RFLAG = .FALSE.   !  true: in road class section
        LOGICAL      :: VFLAG = .FALSE.   !  true: in vehicle type section

        CHARACTER*300          LINE     !  line buffer
        CHARACTER*300          MESG     !  message buffer
        CHARACTER(LEN=VTPLEN3) FIELDS( 2 ) !  fields from vehicle type section

        CHARACTER*16 :: PROGNAME = 'RDMVINFO' ! program name

C***********************************************************************
C   begin body of subroutine RDMVINFO

C.........  Get the number of lines for the file and allocate array so that
C           the type of the line can be stored
        NLINES = GETFLINE( FDEV, 'Mobile codes file' )

        ALLOCATE( LTYPE( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LTYPE', PROGNAME )

        LTYPE = 0  ! array

C.........  Get the number of lines for each section
        N1 = 0
        N2 = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading MOBILE INFORMATION file at line', I
                CALL M3MESG( MESG )
                CYCLE
            END IF

            IF( LINE .EQ. ' ' ) CYCLE  ! Skip blank lines

            K1 = INDEX( LINE, '/VEHICLE TYPES/' )
            K2 = INDEX( LINE, '/ROAD CLASSES/'  )

            IF( .NOT. VFLAG .AND. K1 .GT. 0 ) THEN
                VFLAG = .TRUE.
                RFLAG = .FALSE.

            ELSE IF( .NOT. RFLAG .AND. K2 .GT. 0 ) THEN
                VFLAG = .FALSE.
                RFLAG = .TRUE.

            ELSE IF( VFLAG .AND. K1 .GT. 0 .OR.
     &               RFLAG .AND. K2 .GT. 0     ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad format of MOBILE ' //
     &                 'INFORMATION file at line', I
                CALL M3MESG( MESG )

            ELSE IF( VFLAG ) THEN
                N1 = N1 + 1
                LTYPE( I ) = 1
            
            ELSE IF( RFLAG ) THEN
                N2 = N2 + 1
                LTYPE( I ) = 2

            END IF

        END DO

        REWIND( FDEV )
        NVTYPE = N1
        NRCLAS = N2

        IF( EFLAG ) THEN

            MESG = 'Problem reading mobile information file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Allocate memory for storing unsorted vehicle types and road classes
        ALLOCATE( IVTIDA( NVTYPE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IVTIDA', PROGNAME )
        ALLOCATE( CVTYPA( NVTYPE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CVTYPA', PROGNAME )
        ALLOCATE( VIDXA( NVTYPE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VIDXA', PROGNAME )

        ALLOCATE( RDCLSA( NRCLAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RDCLSA', PROGNAME )
        ALLOCATE( RDWAYA( NRCLAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RDWAYA', PROGNAME )
        ALLOCATE( RIDXA( NRCLAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RIDXA', PROGNAME )

C.........  Store contents of this file in unsorted arrays
        N1 = 0
        N2 = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading MOBILE INFORMATION file at line', I
                CALL M3MESG( MESG )
                CYCLE
            END IF

            IF( LTYPE( I ) .EQ. 1 ) THEN
                N1 = N1 + 1
                VIDXA ( N1 ) = N1
                CALL PARSLINE( LINE, 2, FIELDS )
                CVTYPA( N1 ) = FIELDS( 1 )
                IVTIDA( N1 ) = STR2INT( FIELDS( 2 ) )

            ELSE IF ( LTYPE( I ) .EQ. 2 ) THEN
                N2 = N2 + 1
                RIDXA ( N2 ) = N2
                READ( LINE, * ) RDCLSA( N2 ), RDWAYA( N2 )

            END IF

        END DO 

C.........  Allocate memory for sorted arrays
        ALLOCATE( CVTYPLST( NVTYPE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CVTYPLST', PROGNAME )
        ALLOCATE( IVTIDLST( NVTYPE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IVTIDLST', PROGNAME )

        ALLOCATE( AMSRDCLS( NRCLAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AMSRDCLS', PROGNAME )
        ALLOCATE( RDWAYTYP( NRCLAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RDWAYTYP', PROGNAME )

C.........  Sort vehicle types by name for IDA input checking
        CALL SORTIC( NVTYPE, VIDXA, CVTYPA )

C.........  Sort road classes by number for EPS input checking
        CALL SORTI1( NRCLAS, RIDXA, RDCLSA )

C.........  Store sorted vehicle types and numbers
        DO I = 1, NVTYPE
            J = VIDXA( I )
            CVTYPLST( I ) = CVTYPA( J )
            IVTIDLST( I ) = IVTIDA( J )
        END DO

C.........  Store sorted road classes and roadway types
        DO I = 1, NRCLAS
            J = RIDXA( I )
            AMSRDCLS( I ) = RDCLSA( J )
            RDWAYTYP( I ) = RDWAYA( J )
        END DO

C.........  Deallocate temporary unsorted arrays
        DEALLOCATE( CVTYPA, VIDXA, RDCLSA, RDWAYA, RIDXA )

C.........  Rewind file
        REWIND( FDEV )

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of mobile information file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDMVINFO
