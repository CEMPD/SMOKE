
        SUBROUTINE ALOCUTBL( NUOVAR, ICSIZE )
        
C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine allocates memory for the uncertainty cross-reference
C      table.
C      
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  
C      Subrountines:  CHECKMEM
C
C  REVISION  HISTORY:
C     Created 9/2001 by A. Holland
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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

C...........   This module is for cross reference tables
        USE MODXREF

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: NUOVAR        ! number of pollutants + activities
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups


C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER*16 :: PROGNAME = 'ALOCUTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCUTBL

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( ICTL02 ) ) THEN

            DEALLOCATE( ICTL02, ICTL03, ICTL04, ICTL05, ICTL06 )
            DEALLOCATE( ICTL07, ICTL08, ICTL09, ICTL10, ICTL11 )
            DEALLOCATE( ICTL12, ICTL13, ICTL14, ICTL15, ICTL16 )

        END IF

        ALLOCATE( ICTL01( NUOVAR ), STAT=IOS )         ! SCC=0, FIP=0
        CALL CHECKMEM( IOS, 'ICTL01', PROGNAME )
        ICTL01 = IMISS3

        J = ICSIZE( 2 )                               ! SCC=left, FIP=0
        ALLOCATE( ICTL02( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL02', PROGNAME )
        ICTL02 = IMISS3

        J = ICSIZE( 3 )                               ! SCC=all, FIP=0
        ALLOCATE( ICTL03( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL03', PROGNAME )
        ICTL03 = IMISS3
  
        J = ICSIZE( 4 )                               ! SCC=0, FIP=state
        ALLOCATE( ICTL04( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL04', PROGNAME )
        ICTL04 = IMISS3

        J = ICSIZE( 5 )                               ! SCC=left, FIP=state
        ALLOCATE( ICTL05( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL05', PROGNAME )
        ICTL05 = IMISS3
            
        J = ICSIZE( 6 )                               ! SCC=all, FIP=state
        ALLOCATE( ICTL06( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL06', PROGNAME )
        ICTL06 = IMISS3
                        
        J = ICSIZE( 7 )                               ! SCC=0, FIP=all
        ALLOCATE( ICTL07( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL07', PROGNAME )
        ICTL07 = IMISS3
            
        J = ICSIZE( 8 )                               ! SCC=left, FIP=all
        ALLOCATE( ICTL08( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL08', PROGNAME )
        ICTL08 = IMISS3
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( ICTL09( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL09', PROGNAME )
        ICTL09 = IMISS3
            
        J = ICSIZE( 10 )                              ! PLANT=non-blank, SCC=0
        ALLOCATE( ICTL10( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL10', PROGNAME )
        ICTL10 = IMISS3
            
        J = ICSIZE( 11 )                              ! PLANT=non-blank, SCC=all     
        ALLOCATE( ICTL11( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL11', PROGNAME )
        ICTL11 = IMISS3
            
        J = ICSIZE( 12 )                              ! CHAR1=non-blank, SCC=all     
        ALLOCATE( ICTL12( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL12', PROGNAME )
        ICTL12 = IMISS3
            
        J = ICSIZE( 13 )                              ! CHAR2=non-blank, SCC=all
        ALLOCATE( ICTL13( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL13', PROGNAME )
        ICTL13 = IMISS3
            
        J = ICSIZE( 14 )                              ! CHAR3=non-blank, SCC=all
        ALLOCATE( ICTL14( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL14', PROGNAME )
        ICTL14 = IMISS3
          
        J = ICSIZE( 15 )                              ! CHAR4=non-blank, SCC=all
        ALLOCATE( ICTL15( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL15', PROGNAME )
        ICTL15 = IMISS3
            
        J = ICSIZE( 16 )                              ! CHAR5=non-blank, SCC=all
        ALLOCATE( ICTL16( J,NUOVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL16', PROGNAME )
        ICTL16 = IMISS3
            
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCUTBL
