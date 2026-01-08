
        SUBROUTINE ALOCTTBL( NIPPA, ICSIZE )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the temporal 
C      cross-reference tables that contain the temporal profile numbers, and it
C      initializes these to missing.  The subroutine arguments are the number
C      of inventory pollutants + activities and an array that contains the 
C      dimensions for each of the different groups of the cross-reference.
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

C...........   This module is for cross reference tables
        USE M3UTILIO

        USE MODXREF, ONLY: MPRT01, WPRT01, DPRT01,
     &                     MPRT02, WPRT02, DPRT02,
     &                     MPRT03, WPRT03, DPRT03,
     &                     MPRT04, WPRT04, DPRT04,
     &                     MPRT05, WPRT05, DPRT05,
     &                     MPRT06, WPRT06, DPRT06,
     &                     MPRT07, WPRT07, DPRT07,
     &                     MPRT08, WPRT08, DPRT08,
     &                     MPRT09, WPRT09, DPRT09,
     &                     MPRT10, WPRT10, DPRT10,
     &                     MPRT11, WPRT11, DPRT11,
     &                     MPRT12, WPRT12, DPRT12,
     &                     MPRT13, WPRT13, DPRT13,
     &                     MPRT14, WPRT14, DPRT14,
     &                     MPRT15, WPRT15, DPRT15,
     &                     MPRT16, WPRT16, DPRT16

        IMPLICIT NONE

C...........   INCLUDES

C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: NIPPA        ! number of pollutants + activities
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER(16) :: PROGNAME = 'ALOCTTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCTTBL

        ALLOCATE( MPRT01( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT01', PROGNAME )
        ALLOCATE( WPRT01( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT01', PROGNAME )
        ALLOCATE( DPRT01( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT01', PROGNAME )
        MPRT01 = IMISS3 ! arrays
        WPRT01 = IMISS3
        DPRT01 = IMISS3

        J = ICSIZE( 2 )                                       ! SCC=left, FIP=0
        ALLOCATE( MPRT02( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT02', PROGNAME )
        ALLOCATE( WPRT02( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT02', PROGNAME )
        ALLOCATE( DPRT02( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT02', PROGNAME )
        IF( J .GT. 0 ) THEN
            MPRT02 = IMISS3 ! arrays
            WPRT02 = IMISS3
            DPRT02 = IMISS3
        ENDIF

        J = ICSIZE( 3 )                                   ! SCC=all, FIP=0
        ALLOCATE( MPRT03( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT03', PROGNAME )
        ALLOCATE( WPRT03( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT03', PROGNAME )
        ALLOCATE( DPRT03( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT03', PROGNAME )
        IF( J .GT. 0 ) THEN
            MPRT03 = IMISS3 ! arrays
            WPRT03 = IMISS3
            DPRT03 = IMISS3
        ENDIF
                
        J = ICSIZE( 4 )                                 ! SCC=0, FIP=state
        ALLOCATE( MPRT04( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT04', PROGNAME )
        ALLOCATE( WPRT04( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT04', PROGNAME )
        ALLOCATE( DPRT04( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT04', PROGNAME )
        IF( J .GT. 0 ) THEN
            MPRT04 = IMISS3 ! arrays
            WPRT04 = IMISS3
            DPRT04 = IMISS3
        ENDIF
            
        J = ICSIZE( 5 )                                 ! SCC=left, FIP=state
        ALLOCATE( MPRT05( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT05', PROGNAME )
        ALLOCATE( WPRT05( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT05', PROGNAME )
        ALLOCATE( DPRT05( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT05', PROGNAME )
        IF( J .GT. 0 ) THEN
            MPRT05 = IMISS3 ! arrays
            WPRT05 = IMISS3
            DPRT05 = IMISS3
        ENDIF
            
        J = ICSIZE( 6 )                          ! SCC=all, FIP=state
        ALLOCATE( MPRT06( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT06', PROGNAME )
        ALLOCATE( WPRT06( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT06', PROGNAME )
        ALLOCATE( DPRT06( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT06', PROGNAME )
        IF( J .GT. 0 ) THEN
            MPRT06 = IMISS3 ! arrays
            WPRT06 = IMISS3
            DPRT06 = IMISS3
        ENDIF
                        
        J = ICSIZE( 7 )                             ! SCC=0, FIP=all
        ALLOCATE( MPRT07( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT07', PROGNAME )
        ALLOCATE( WPRT07( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT07', PROGNAME )
        ALLOCATE( DPRT07( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT07', PROGNAME )
        IF( J .GT. 0 ) THEN
            MPRT07 = IMISS3 ! arrays
            WPRT07 = IMISS3
            DPRT07 = IMISS3
        ENDIF
            
        J = ICSIZE( 8 )                         ! SCC=left, FIP=all
        ALLOCATE( MPRT08( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT08', PROGNAME )
        ALLOCATE( WPRT08( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT08', PROGNAME )
        ALLOCATE( DPRT08( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT08', PROGNAME )
        IF( J .GT. 0 ) THEN
            MPRT08 = IMISS3 ! arrays
            WPRT08 = IMISS3
            DPRT08 = IMISS3
        ENDIF
                        
        J = ICSIZE( 9 )                         ! SCC=all, FIP=all
        ALLOCATE( MPRT09( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT09', PROGNAME )
        ALLOCATE( WPRT09( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT09', PROGNAME )
        ALLOCATE( DPRT09( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT09', PROGNAME )
        IF( J .GT. 0 ) THEN 
            MPRT09 = IMISS3 ! arrays
            WPRT09 = IMISS3
            DPRT09 = IMISS3
        ENDIF
            
        J = ICSIZE( 10 )                       ! PLANT=non-blank, SCC=0
        ALLOCATE( MPRT10( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT10', PROGNAME )
        ALLOCATE( WPRT10( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT10', PROGNAME )
        ALLOCATE( DPRT10( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT10', PROGNAME )

        IF( J .GT. 0 ) THEN
            MPRT10 = IMISS3 ! arrays
            WPRT10 = IMISS3
            DPRT10 = IMISS3
        ENDIF
            
        J = ICSIZE( 11 )                          ! PLANT=non-blank, SCC=all    
        ALLOCATE( MPRT11( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT11', PROGNAME )
        ALLOCATE( WPRT11( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT11', PROGNAME )
        ALLOCATE( DPRT11( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT11', PROGNAME )
        IF( J .GT. 0 ) THEN 
            MPRT11 = IMISS3 ! arrays
            WPRT11 = IMISS3
            DPRT11 = IMISS3
        ENDIF
            
        J = ICSIZE( 12 )                         ! CHAR1=non-blank, SCC=all     
        ALLOCATE( MPRT12( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT12', PROGNAME )
        ALLOCATE( WPRT12( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT12', PROGNAME )
        ALLOCATE( DPRT12( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT12', PROGNAME )
        IF( J .GT. 0 ) THEN
            MPRT12 = IMISS3 ! arrays
            WPRT12 = IMISS3
            DPRT12 = IMISS3
        ENDIF
            
        J = ICSIZE( 13 )                        ! CHAR2=non-blank, SCC=all
        ALLOCATE( MPRT13( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT13', PROGNAME )
        ALLOCATE( WPRT13( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT13', PROGNAME )
        ALLOCATE( DPRT13( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT13', PROGNAME )
        IF( J .GT. 0 ) THEN
            MPRT13 = IMISS3 ! arrays
            WPRT13 = IMISS3
            DPRT13 = IMISS3
        ENDIF
            
        J = ICSIZE( 14 )                     ! CHAR3=non-blank, SCC=all
        ALLOCATE( MPRT14( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT14', PROGNAME )
        ALLOCATE( WPRT14( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT14', PROGNAME )
        ALLOCATE( DPRT14( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT14', PROGNAME )
        IF( J .GT. 0 ) THEN 
            MPRT14 = IMISS3 ! arrays
            WPRT14 = IMISS3
            DPRT14 = IMISS3
        ENDIF
            
        J = ICSIZE( 15 )                     ! CHAR4=non-blank, SCC=all
        ALLOCATE( MPRT15( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT15', PROGNAME )
        ALLOCATE( WPRT15( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT15', PROGNAME )
        ALLOCATE( DPRT15( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT15', PROGNAME )
        IF( J .GT. 0 ) THEN 
            MPRT15 = IMISS3 ! arrays
            WPRT15 = IMISS3
            DPRT15 = IMISS3
        ENDIF
            
        J = ICSIZE( 16 )                     ! CHAR5=non-blank, SCC=all
        ALLOCATE( MPRT16( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRT16', PROGNAME )
        ALLOCATE( WPRT16( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRT16', PROGNAME )
        ALLOCATE( DPRT16( J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRT16', PROGNAME )
        IF( J .GT. 0 ) THEN 
            MPRT16 = IMISS3 ! arrays
            WPRT16 = IMISS3
            DPRT16 = IMISS3
        ENDIF
            
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCTTBL
