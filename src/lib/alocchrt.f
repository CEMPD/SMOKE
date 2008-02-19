
        SUBROUTINE ALOCCHRT( ICSIZE )

C***********************************************************************
C  subroutine body starts at line 60
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
C***************************************************************************
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
        USE MODXREF, ONLY: CHRT02, CHRT03, CHRT04, CHRT05, CHRT06,
     &      CHRT07, CHRT08, CHRT09, CHRT10, CHRT11,
     &      CHRT12, CHRT13, CHRT14, CHRT15, CHRT16,
     &      CHRT02A, CHRT02B, CHRT02C,
     &      CHRT05A, CHRT05B, CHRT05C,
     &      CHRT08A, CHRT08B, CHRT08C,
     &      CHRT26, CHRT27, CHRT28, CHRT29, CHRT30, CHRT31,
     &      CHRT32, CHRT33, CHRT34, CHRT35, CHRT36, CHRT37, CHRT38

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER(16) :: PROGNAME = 'ALOCCHRT' ! program name

C***********************************************************************
C   begin body of subroutine ALOCCHRT

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( CHRT02 ) ) THEN

            DEALLOCATE( CHRT02, CHRT03, CHRT04, CHRT05, CHRT06 )
            DEALLOCATE( CHRT07, CHRT08, CHRT09, CHRT10, CHRT11 )
            DEALLOCATE( CHRT12, CHRT13, CHRT14, CHRT15, CHRT16 )
            DEALLOCATE( CHRT02A, CHRT02B, CHRT02C )
            DEALLOCATE( CHRT05A, CHRT05B, CHRT05C )
            DEALLOCATE( CHRT08A, CHRT08B, CHRT08C )
            DEALLOCATE( CHRT26, CHRT27, CHRT28, CHRT29, CHRT30, CHRT31 )
            DEALLOCATE( CHRT32, CHRT33, CHRT34, CHRT35, CHRT36, CHRT37 )
            DEALLOCATE( CHRT38 )

        END IF

        J = MAX( 1, ICSIZE( 2 ) )                     ! SCC=left, FIP=0
        ALLOCATE( CHRT02( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT02', PROGNAME )

        J = MAX( 1, ICSIZE( 3 ) )                     ! SCC=all, FIP=0
        ALLOCATE( CHRT03( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT03', PROGNAME )
                
        J = MAX( 1, ICSIZE( 4 ) )                     ! SCC=0, FIP=state
        ALLOCATE( CHRT04( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT04', PROGNAME )
            
        J = MAX( 1, ICSIZE( 5 ) )                     ! SCC=left, FIP=state
        ALLOCATE( CHRT05( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT05', PROGNAME )
            
        J = MAX( 1, ICSIZE( 6 ) )                     ! SCC=all, FIP=state
        ALLOCATE( CHRT06( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT06', PROGNAME )
                        
        J = MAX( 1, ICSIZE( 7 ) )                     ! SCC=0, FIP=all
        ALLOCATE( CHRT07( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT07', PROGNAME )
            
        J = MAX( 1, ICSIZE( 8 ) )                     ! SCC=left, FIP=all
        ALLOCATE( CHRT08( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT08', PROGNAME )
                        
        J = MAX( 1, ICSIZE( 9 ) )                     ! SCC=all, FIP=all
        ALLOCATE( CHRT09( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT09', PROGNAME )
            
        J = MAX( 1, ICSIZE( 10 ) )                    ! PLANT=non-blank, SCC=0
        ALLOCATE( CHRT10( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT10', PROGNAME )
            
        J = MAX( 1, ICSIZE( 11 ) )                    ! PLANT=non-blank, SCC=all
        ALLOCATE( CHRT11( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT11', PROGNAME )
            
        J = MAX( 1, ICSIZE( 12 ) )                    ! CHAR1=non-blank, SCC=all
        ALLOCATE( CHRT12( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT12', PROGNAME )
            
        J = MAX( 1, ICSIZE( 13 ) )                    ! CHAR2=non-blank, SCC=all
        ALLOCATE( CHRT13( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT13', PROGNAME )
            
        J = MAX( 1, ICSIZE( 14 ) )                    ! CHAR3=non-blank, SCC=all
        ALLOCATE( CHRT14( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT14', PROGNAME )
          
        J = MAX( 1, ICSIZE( 15 ) )                    ! CHAR4=non-blank, SCC=all
        ALLOCATE( CHRT15( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT15', PROGNAME )
            
        J = MAX( 1, ICSIZE( 16 ) )                    ! CHAR5=non-blank, SCC=all
        ALLOCATE( CHRT16( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT16', PROGNAME )

C.........  NOTE- Added later            
        J = MAX( 1, ICSIZE( 17 ) )                     ! SCC=level 1, FIP=0
        ALLOCATE( CHRT02A( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT02A', PROGNAME )

        J = MAX( 1, ICSIZE( 18 ) )                     ! SCC=level 2, FIP=0
        ALLOCATE( CHRT02B( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT02B', PROGNAME )

        J = MAX( 1, ICSIZE( 19 ) )                     ! SCC=level 3, FIP=0
        ALLOCATE( CHRT02C( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT02C', PROGNAME )

        J = MAX( 1, ICSIZE( 20 ) )                     ! SCC=level 1, FIP=state
        ALLOCATE( CHRT05A( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT05A', PROGNAME )

        J = MAX( 1, ICSIZE( 21 ) )                     ! SCC=level 2, FIP=state
        ALLOCATE( CHRT05B( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT05B', PROGNAME )

        J = MAX( 1, ICSIZE( 22 ) )                     ! SCC=level 3, FIP=state
        ALLOCATE( CHRT05C( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT05C', PROGNAME )

        J = MAX( 1, ICSIZE( 23 ) )                     ! SCC=level 1, FIP=all
        ALLOCATE( CHRT08A( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT08A', PROGNAME )

        J = MAX( 1, ICSIZE( 24 ) )                     ! SCC=level 2, FIP=all
        ALLOCATE( CHRT08B( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT08B', PROGNAME )

        J = MAX( 1, ICSIZE( 25 ) )                     ! SCC=level 3, FIP=all
        ALLOCATE( CHRT08C( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT08C', PROGNAME )

        J = MAX( 1, ICSIZE( 26 ) )                     ! SIC=2-digit, FIP=0
        ALLOCATE( CHRT26( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT26', PROGNAME )

        J = MAX( 1, ICSIZE( 27 ) )                     ! SIC=all, FIP=0
        ALLOCATE( CHRT27( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT27', PROGNAME )

        J = MAX( 1, ICSIZE( 28 ) )                     ! SIC=2-digit, FIP=state
        ALLOCATE( CHRT28( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT28', PROGNAME )

        J = MAX( 1, ICSIZE( 29 ) )                     ! SIC=all, FIP=state
        ALLOCATE( CHRT29( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT29', PROGNAME )

        J = MAX( 1, ICSIZE( 30 ) )                     ! SIC=2-digit, FIP=all
        ALLOCATE( CHRT30( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT30', PROGNAME )

        J = MAX( 1, ICSIZE( 31 ) )                     ! SIC=all, FIP=all
        ALLOCATE( CHRT31( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT31', PROGNAME )

C.........  Added most recently
        J = MAX( 1, ICSIZE( 32 ) )                     ! MACT=all, FIP=0, SCC=0
        ALLOCATE( CHRT32( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT32', PROGNAME )
        
        J = MAX( 1, ICSIZE( 33 ) )                     ! MACT=all, FIP=0, SCC=all
        ALLOCATE( CHRT33( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT33', PROGNAME )
        
        J = MAX( 1, ICSIZE( 34 ) )                     ! MACT=all, FIP=state, SCC=0
        ALLOCATE( CHRT34( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT34', PROGNAME )
        
        J = MAX( 1, ICSIZE( 35 ) )                     ! MACT=all, FIP=state, SCC=all
        ALLOCATE( CHRT35( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT35', PROGNAME )
        
        J = MAX( 1, ICSIZE( 36 ) )                     ! MACT=all, FIP=all, SCC=0
        ALLOCATE( CHRT36( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT36', PROGNAME )
        
        J = MAX( 1, ICSIZE( 37 ) )                     ! MACT=all, FIP=all, SCC=all
        ALLOCATE( CHRT37( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT37', PROGNAME )

        J = MAX( 1, ICSIZE( 38 ) )                     ! PLANT=non-blank, MACT=all
        ALLOCATE( CHRT38( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRT38', PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCCHRT
