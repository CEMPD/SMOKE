
        SUBROUTINE ALOCCHRT( ICSIZE )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the grouped 
C      cross-reference tables that contain the source characteristics.
C      
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/99 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER*16 :: PROGNAME = 'ALOCCHRT' ! program name

C***********************************************************************
C   begin body of subroutine ALOCCHRT

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( CHRT02 ) ) THEN

            DEALLOCATE( CHRT02, CHRT03, CHRT04, CHRT05, CHRT06 )
            DEALLOCATE( CHRT07, CHRT08, CHRT09, CHRT10, CHRT11 )
            DEALLOCATE( CHRT12, CHRT13, CHRT14, CHRT15, CHRT16 )

        END IF

        J = ICSIZE( 2 )                               ! SCC=left, FIP=0
        ALLOCATE( CHRT02( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT02', PROGNAME )

        J = ICSIZE( 3 )                               ! SCC=all, FIP=0
        ALLOCATE( CHRT03( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT03', PROGNAME )
                
        J = ICSIZE( 4 )                               ! SCC=0, FIP=state
        ALLOCATE( CHRT04( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT04', PROGNAME )
            
        J = ICSIZE( 5 )                               ! SCC=left, FIP=state
        ALLOCATE( CHRT05( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT05', PROGNAME )
            
        J = ICSIZE( 6 )                               ! SCC=all, FIP=state
        ALLOCATE( CHRT06( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT06', PROGNAME )
                        
        J = ICSIZE( 7 )                               ! SCC=0, FIP=all
        ALLOCATE( CHRT07( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT07', PROGNAME )
            
        J = ICSIZE( 8 )                               ! SCC=left, FIP=all
        ALLOCATE( CHRT08( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT08', PROGNAME )
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( CHRT09( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT09', PROGNAME )
            
        J = ICSIZE( 10 )                              ! PLANT=non-blank, SCC=0
        ALLOCATE( CHRT10( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT10', PROGNAME )
            
        J = ICSIZE( 11 )                              ! PLANT=non-blank, SCC=all     
        ALLOCATE( CHRT11( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT11', PROGNAME )
            
        J = ICSIZE( 12 )                              ! CHAR1=non-blank, SCC=all     
        ALLOCATE( CHRT12( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT12', PROGNAME )
            
        J = ICSIZE( 13 )                              ! CHAR2=non-blank, SCC=all
        ALLOCATE( CHRT13( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT13', PROGNAME )
            
        J = ICSIZE( 14 )                              ! CHAR3=non-blank, SCC=all
        ALLOCATE( CHRT14( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT14', PROGNAME )
          
        J = ICSIZE( 15 )                              ! CHAR4=non-blank, SCC=all
        ALLOCATE( CHRT15( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT15', PROGNAME )
            
        J = ICSIZE( 16 )                              ! CHAR5=non-blank, SCC=all
        ALLOCATE( CHRT16( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT16', PROGNAME )
            
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCCHRT
