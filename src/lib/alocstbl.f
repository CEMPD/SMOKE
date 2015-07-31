
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

C...........   This module is for cross reference tables
        USE MODXREF, ONLY: CSPT01, CSPT02, CSPT03, CSPT04, CSPT05,
     &          CSPT06, CSPT07, CSPT08, CSPT09, CSPT10, CSPT11,
     &          CSPT12, CSPT13, CSPT14, CSPT15, CSPT16,
     &          CSPT26, CSPT27, CSPT28, CSPT29, CSPT30, CSPT31,
     &          CSPT32, CSPT33, CSPT34, CSPT35, CSPT36, CSPT37

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: NIPOL        ! number of pollutants and emis type
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER(16) :: PROGNAME = 'ALOCSTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCSTBL

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( CSPT02 ) ) THEN

            DEALLOCATE( CSPT02, CSPT03, CSPT04, CSPT05, CSPT06 )
            DEALLOCATE( CSPT07, CSPT08, CSPT09, CSPT10, CSPT11 )
            DEALLOCATE( CSPT12, CSPT13, CSPT14, CSPT15, CSPT16 )
            DEALLOCATE( CSPT26, CSPT27, CSPT28, CSPT29, CSPT30, CSPT31 )
            DEALLOCATE( CSPT32, CSPT33, CSPT34, CSPT35, CSPT36, CSPT37 )

        END IF

        ALLOCATE( CSPT01( NIPOL ), STAT=IOS )         ! SCC=0, FIP=0
        CALL CHECKMEM( IOS, 'CSPT01', PROGNAME )
        CSPT01 = EMCMISS3

        J = ICSIZE( 2 )                               ! SCC=left, FIP=0
        ALLOCATE( CSPT02( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT02', PROGNAME )
        CSPT02 = EMCMISS3

        J = ICSIZE( 3 )                               ! SCC=all, FIP=0
        ALLOCATE( CSPT03( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT03', PROGNAME )
        CSPT03 = EMCMISS3
                
        J = ICSIZE( 4 )                               ! SCC=0, FIP=state
        ALLOCATE( CSPT04( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT04', PROGNAME )
        CSPT04 = EMCMISS3
            
        J = ICSIZE( 5 )                               ! SCC=left, FIP=state
        ALLOCATE( CSPT05( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT05', PROGNAME )
        CSPT05 = EMCMISS3
            
        J = ICSIZE( 6 )                               ! SCC=all, FIP=state
        ALLOCATE( CSPT06( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT06', PROGNAME )
        CSPT06 = EMCMISS3
                        
        J = ICSIZE( 7 )                               ! SCC=0, FIP=all
        ALLOCATE( CSPT07( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT07', PROGNAME )
        CSPT07 = EMCMISS3
            
        J = ICSIZE( 8 )                               ! SCC=left, FIP=all
        ALLOCATE( CSPT08( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT08', PROGNAME )
        CSPT08 = EMCMISS3
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( CSPT09( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT09', PROGNAME )
        CSPT09 = EMCMISS3
            
        J = ICSIZE( 10 )                              ! PLANT=non-blank, SCC=0
        ALLOCATE( CSPT10( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT10', PROGNAME )
        CSPT10 = EMCMISS3
            
        J = ICSIZE( 11 )                              ! PLANT=non-blank, SCC=all     
        ALLOCATE( CSPT11( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT11', PROGNAME )
        CSPT11 = EMCMISS3
            
        J = ICSIZE( 12 )                              ! CHAR1=non-blank, SCC=all     
        ALLOCATE( CSPT12( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT12', PROGNAME )
        CSPT12 = EMCMISS3
         
        J = ICSIZE( 13 )                              ! CHAR2=non-blank, SCC=all
        ALLOCATE( CSPT13( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT13', PROGNAME )
        CSPT13 = EMCMISS3
            
        J = ICSIZE( 14 )                              ! CHAR3=non-blank, SCC=all
        ALLOCATE( CSPT14( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT14', PROGNAME )
        CSPT14 = EMCMISS3
          
        J = ICSIZE( 15 )                              ! CHAR4=non-blank, SCC=all
        ALLOCATE( CSPT15( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT15', PROGNAME )
        CSPT15 = EMCMISS3
            
        J = ICSIZE( 16 )                              ! CHAR5=non-blank, SCC=all
        ALLOCATE( CSPT16( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT16', PROGNAME )
        CSPT16 = EMCMISS3

C.........  SIC code matches

        J = ICSIZE( 26 )                              ! FIP=0, SIC=left
        ALLOCATE( CSPT26( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT26', PROGNAME )
        CSPT26 = EMCMISS3
        
        J = ICSIZE( 27 )                              ! FIP=0, SIC=all
        ALLOCATE( CSPT27( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT27', PROGNAME )
        CSPT27 = EMCMISS3
        
        J = ICSIZE( 28 )                              ! FIP=state, SIC=left
        ALLOCATE( CSPT28( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT28', PROGNAME )
        CSPT28 = EMCMISS3
        
        J = ICSIZE( 29 )                              ! FIP=state, SIC=all
        ALLOCATE( CSPT29( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT29', PROGNAME )
        CSPT29 = EMCMISS3
        
        J = ICSIZE( 30 )                              ! FIP=all, SIC=left
        ALLOCATE( CSPT30( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT30', PROGNAME )
        CSPT30 = EMCMISS3
        
        J = ICSIZE( 31 )                              ! FIP=all, SIC=all
        ALLOCATE( CSPT31( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT31', PROGNAME )
        CSPT31 = EMCMISS3

C.........  MACT code matches
        
        J = ICSIZE( 32 )                              ! FIP=0, SCC=0, MACT=all
        ALLOCATE( CSPT32( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT32', PROGNAME )
        CSPT32 = EMCMISS3
        
        J = ICSIZE( 33 )                              ! FIP=0, SCC=all, MACT=all
        ALLOCATE( CSPT33( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT33', PROGNAME )
        CSPT33 = EMCMISS3
        
        J = ICSIZE( 34 )                              ! FIP=state, SCC=0, MACT=all
        ALLOCATE( CSPT34( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT34', PROGNAME )
        CSPT34 = EMCMISS3
        
        J = ICSIZE( 35 )                              ! FIP=state, SCC=all, MACT=all
        ALLOCATE( CSPT35( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT35', PROGNAME )
        CSPT35 = EMCMISS3
        
        J = ICSIZE( 36 )                              ! FIP=all, SCC=0, MACT=all
        ALLOCATE( CSPT36( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT36', PROGNAME )
        CSPT36 = EMCMISS3
        
        J = ICSIZE( 37 )                              ! FIP=all, SCC=all, MACT=all
        ALLOCATE( CSPT37( -1:J,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPT37', PROGNAME )
        CSPT37 = EMCMISS3
            
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCSTBL
