
        SUBROUTINE ALOCSTBL( NIPOL, ICSIZE )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the speciation
C      cross-reference tables that contain the speciation profile numbers, and
C      it initializes these to missing.  The subroutine arguments are the number
C      of inventory pollutants and an array that contains the dimensions for
C      each of the different groups of the cross-reference.
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
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   This module is for cross reference tables
        USE MODXREF

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: NIPOL        ! number of pollutants and emis type
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER*16 :: PROGNAME = 'ALOCSTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCSTBL

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( CSPT02 ) ) THEN

            DEALLOCATE( CSPT02, CSPT03, CSPT04, CSPT05, CSPT06 )
            DEALLOCATE( CSPT07, CSPT08, CSPT09, CSPT10, CSPT11 )
            DEALLOCATE( CSPT12, CSPT13, CSPT14, CSPT15, CSPT16 )

        END IF

        ALLOCATE( CSPT01( NIPOL ), STAT=IOS )         ! SCC=0, FIP=0
        CALL CHECKMEM( IOS, 'CSPT01', PROGNAME )
        CSPT01 = EMCMISS3

        J = ICSIZE( 2 )                               ! SCC=left, FIP=0
        ALLOCATE( CSPT02( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT02', PROGNAME )
        CSPT02 = EMCMISS3

        J = ICSIZE( 3 )                               ! SCC=all, FIP=0
        ALLOCATE( CSPT03( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT03', PROGNAME )
        CSPT03 = EMCMISS3
                
        J = ICSIZE( 4 )                               ! SCC=0, FIP=state
        ALLOCATE( CSPT04( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT04', PROGNAME )
        CSPT04 = EMCMISS3
            
        J = ICSIZE( 5 )                               ! SCC=left, FIP=state
        ALLOCATE( CSPT05( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT05', PROGNAME )
        CSPT05 = EMCMISS3
            
        J = ICSIZE( 6 )                               ! SCC=all, FIP=state
        ALLOCATE( CSPT06( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT06', PROGNAME )
        CSPT06 = EMCMISS3
                        
        J = ICSIZE( 7 )                               ! SCC=0, FIP=all
        ALLOCATE( CSPT07( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT07', PROGNAME )
        CSPT07 = EMCMISS3
            
        J = ICSIZE( 8 )                               ! SCC=left, FIP=all
        ALLOCATE( CSPT08( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT08', PROGNAME )
        CSPT08 = EMCMISS3
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( CSPT09( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT09', PROGNAME )
        CSPT09 = EMCMISS3
            
        J = ICSIZE( 10 )                              ! PLANT=non-blank, SCC=0
        ALLOCATE( CSPT10( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT10', PROGNAME )
        CSPT10 = EMCMISS3
            
        J = ICSIZE( 11 )                              ! PLANT=non-blank, SCC=all     
        ALLOCATE( CSPT11( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT11', PROGNAME )
        CSPT11 = EMCMISS3
            
        J = ICSIZE( 12 )                              ! CHAR1=non-blank, SCC=all     
        ALLOCATE( CSPT12( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT12', PROGNAME )
        CSPT12 = EMCMISS3
         
        J = ICSIZE( 13 )                              ! CHAR2=non-blank, SCC=all
        ALLOCATE( CSPT13( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT13', PROGNAME )
        CSPT13 = EMCMISS3
            
        J = ICSIZE( 14 )                              ! CHAR3=non-blank, SCC=all
        ALLOCATE( CSPT14( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT14', PROGNAME )
        CSPT14 = EMCMISS3
          
        J = ICSIZE( 15 )                              ! CHAR4=non-blank, SCC=all
        ALLOCATE( CSPT15( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT15', PROGNAME )
        CSPT15 = EMCMISS3
            
        J = ICSIZE( 16 )                              ! CHAR5=non-blank, SCC=all
        ALLOCATE( CSPT16( J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT16', PROGNAME )
        CSPT16 = EMCMISS3
            
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCSTBL
