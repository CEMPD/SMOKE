
        SUBROUTINE ALOCCTBL( NIPPA, ICSIZE )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the control 
C      cross-reference tables that contain the index to the control data, and 
C      it initializes these to missing.  The subroutine arguments are the number
C      of inventory pollutants + activities and an array that contains the 
C      dimensions for each of the different groups of the cross-reference.
C      Note that these tables are used multiple times in the same program for 
C      different control packets, which are processed one at a time.
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
        INTEGER, INTENT(IN) :: NIPPA        ! number of pollutants + activities
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER*16 :: PROGNAME = 'ALOCCTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCCTBL

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( ICTL01 ) ) THEN

            DEALLOCATE( ICTL01, ICTL02, ICTL03, ICTL04, ICTL05, ICTL06 )
            DEALLOCATE( ICTL07, ICTL08, ICTL09, ICTL10, ICTL11 )
            DEALLOCATE( ICTL12, ICTL13, ICTL14, ICTL15, ICTL16 )
            DEALLOCATE( ICTL02A, ICTL02B, ICTL02C )
            DEALLOCATE( ICTL05A, ICTL05B, ICTL05C )
            DEALLOCATE( ICTL08A, ICTL08B, ICTL08C )

        END IF

        ALLOCATE( ICTL01( NIPPA ), STAT=IOS )         ! SCC=0, FIP=0
        CALL CHECKMEM( IOS, 'ICTL01', PROGNAME )
        ICTL01 = IMISS3

        J = ICSIZE( 2 )                               ! SCC=left, FIP=0
        ALLOCATE( ICTL02( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL02', PROGNAME )
        ICTL02 = IMISS3

        J = ICSIZE( 3 )                               ! SCC=all, FIP=0
        ALLOCATE( ICTL03( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL03', PROGNAME )
        ICTL03 = IMISS3
  
        J = ICSIZE( 4 )                               ! SCC=0, FIP=state
        ALLOCATE( ICTL04( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL04', PROGNAME )
        ICTL04 = IMISS3

        J = ICSIZE( 5 )                               ! SCC=left, FIP=state
        ALLOCATE( ICTL05( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL05', PROGNAME )
        ICTL05 = IMISS3
            
        J = ICSIZE( 6 )                               ! SCC=all, FIP=state
        ALLOCATE( ICTL06( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL06', PROGNAME )
        ICTL06 = IMISS3
                        
        J = ICSIZE( 7 )                               ! SCC=0, FIP=all
        ALLOCATE( ICTL07( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL07', PROGNAME )
        ICTL07 = IMISS3
            
        J = ICSIZE( 8 )                               ! SCC=left, FIP=all
        ALLOCATE( ICTL08( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL08', PROGNAME )
        ICTL08 = IMISS3
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( ICTL09( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL09', PROGNAME )
        ICTL09 = IMISS3
            
        J = ICSIZE( 10 )                              ! PLANT=non-blank, SCC=0
        ALLOCATE( ICTL10( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL10', PROGNAME )
        ICTL10 = IMISS3
            
        J = ICSIZE( 11 )                              ! PLANT=non-blank, SCC=all     
        ALLOCATE( ICTL11( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL11', PROGNAME )
        ICTL11 = IMISS3
            
        J = ICSIZE( 12 )                              ! CHAR1=non-blank, SCC=all     
        ALLOCATE( ICTL12( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL12', PROGNAME )
        ICTL12 = IMISS3
            
        J = ICSIZE( 13 )                              ! CHAR2=non-blank, SCC=all
        ALLOCATE( ICTL13( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL13', PROGNAME )
        ICTL13 = IMISS3
            
        J = ICSIZE( 14 )                              ! CHAR3=non-blank, SCC=all
        ALLOCATE( ICTL14( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL14', PROGNAME )
        ICTL14 = IMISS3
          
        J = ICSIZE( 15 )                              ! CHAR4=non-blank, SCC=all
        ALLOCATE( ICTL15( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL15', PROGNAME )
        ICTL15 = IMISS3
            
        J = ICSIZE( 16 )                              ! CHAR5=non-blank, SCC=all
        ALLOCATE( ICTL16( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL16', PROGNAME )
        ICTL16 = IMISS3
            
C.........  NOTE- Added later            
        J = MAX( 1, ICSIZE( 17 ) )                     ! SCC=level 1, FIP=0
        ALLOCATE( ICTL02A( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL02A', PROGNAME )

        J = MAX( 1, ICSIZE( 18 ) )                     ! SCC=level 2, FIP=0
        ALLOCATE( ICTL02B( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL02B', PROGNAME )

        J = MAX( 1, ICSIZE( 19 ) )                     ! SCC=level 3, FIP=0
        ALLOCATE( ICTL02C( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL02C', PROGNAME )

        J = MAX( 1, ICSIZE( 20 ) )                     ! SCC=level 1, FIP=state
        ALLOCATE( ICTL05A( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL05A', PROGNAME )

        J = MAX( 1, ICSIZE( 21 ) )                     ! SCC=level 2, FIP=state
        ALLOCATE( ICTL05B( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL05B', PROGNAME )

        J = MAX( 1, ICSIZE( 22 ) )                     ! SCC=level 3, FIP=state
        ALLOCATE( ICTL05C( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL05C', PROGNAME )

        J = MAX( 1, ICSIZE( 23 ) )                     ! SCC=level 1, FIP=all
        ALLOCATE( ICTL08A( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL08A', PROGNAME )

        J = MAX( 1, ICSIZE( 24 ) )                     ! SCC=level 2, FIP=all
        ALLOCATE( ICTL08B( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL08B', PROGNAME )

        J = MAX( 1, ICSIZE( 25 ) )                     ! SCC=level 3, FIP=all
        ALLOCATE( ICTL08C( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL08C', PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCCTBL
