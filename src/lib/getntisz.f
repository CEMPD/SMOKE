
        SUBROUTINE GETNTISZ( FDEV, CATEGORY, NLINES )
        
C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine returns an exact number of records for a raw toxics
C      inventory input file opened on unit FDEV.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 1/99
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

C...........   MODULES for public variables
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NUNIQCAS, UNIQCAS, UCASNKEP

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER, EXTERNAL :: FINDC

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN ) :: FDEV     !  unit number of input file
        CHARACTER(*), INTENT (IN ) :: CATEGORY !  description of source category
        INTEGER     , INTENT (OUT) :: NLINES   !  total number of inventory records

C...........   Local allocatable arrays
        CHARACTER(50), ALLOCATABLE :: SEGMENT( : ) ! pieces of input line

C...........   Other local variables
        INTEGER         I            ! counter 
        INTEGER         IOS          ! i/o status
        INTEGER      :: IREC    = 0  ! input line counter
        INTEGER         CASPOS       ! position of CAS number in line
        INTEGER         NSEG         ! number of segments in line

        CHARACTER(300)     LINE         ! input file line buffer
        CHARACTER(CASLEN3) CASNUM       ! CAS number from input line  
        CHARACTER(300)     MESG         ! message buffer

        CHARACTER(16) ::   PROGNAME = 'GETNTISZ' ! program name

C***********************************************************************
C   begin body of subroutine GETNTISZ

C.........  Write message to screen
        CALL M3MSG2( 'Determining size of toxics inventory...' )

C.........  Initialize counters
        NLINES = 0
        IREC = 0

C.........  Allocate memory for line segments
        SELECT CASE( CATEGORY )
        CASE( 'AREA' )
            NSEG = 9
            CASPOS = 4
        
        CASE( 'MOBILE' )
            NSEG = 6
            CASPOS = 4
            
        END SELECT
        
        IF( .NOT. ALLOCATED( SEGMENT ) ) THEN
            ALLOCATE( SEGMENT( NSEG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )
        END IF
        SEGMENT = ' '  ! array
        
C.........  Loop through lines in file
        DO
 
C.............  Read line of file
            READ( FDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1
 
C.............  Check I/O error status
            IF( IOS > 0 ) THEN
                WRITE( MESG, 94010 )
     &                     'Error', IOS,  'reading inventory file ' // 
     &                     'as character strings at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSEIF( IOS < 0 ) THEN  ! reached end of file
                EXIT

            END IF

C.............  Skip blank lines
            IF( LINE == ' ' ) CYCLE

C.............  Skip header lines
            IF( LINE( 1:1 ) == CINVHDR ) CYCLE        

C.............  Extract CAS number from line
            CALL PARSLINE( LINE, NSEG, SEGMENT )
            CASNUM = ADJUSTL( SEGMENT( CASPOS ) )

C.............  Find CAS number in unique sorted list
            I = FINDC( CASNUM, NUNIQCAS, UNIQCAS )
            IF( I < 1 ) CYCLE  ! skip entries with invalid CAS numbers
            
C.............  Add number of records for this CAS
            NLINES = NLINES + UCASNKEP( I )

        END DO

        IF( NLINES .EQ. 0 ) THEN
            MESG = 'Inventory file has no valid ' //
     &             'lines of inventory data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        REWIND( FDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GETNTISZ
