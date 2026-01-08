
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

        USE MODXREF, ONLY: ICTL01, ICTL02, ICTL03, ICTL04, ICTL05,
     &          ICTL06, ICTL07, ICTL08, ICTL09, ICTL10,
     &          ICTL11, ICTL12, ICTL13, ICTL14, ICTL15, ICTL16,
     &          ICTL02A, ICTL02B, ICTL02C,
     &          ICTL05A, ICTL05B, ICTL05C,
     &          ICTL08A, ICTL08B, ICTL08C,
     &          ICTL26, ICTL27, ICTL28, ICTL29, ICTL30, ICTL31,
     &          ICTL32, ICTL33, ICTL34, ICTL35, ICTL36, ICTL37, ICTL38

        IMPLICIT NONE

C...........   INCLUDES
C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: NIPPA        ! number of pollutants + activities
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER(16) :: PROGNAME = 'ALOCCTBL' ! program name

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
            DEALLOCATE( ICTL26, ICTL27, ICTL28, ICTL29, ICTL30, ICTL31 )
            DEALLOCATE( ICTL32, ICTL33, ICTL34, ICTL35, ICTL36, ICTL37 )
            DEALLOCATE( ICTL38 )

        END IF

        ALLOCATE( ICTL01( NIPPA ), STAT=IOS )         ! SCC=0, FIP=0
        CALL CHECKMEM( IOS, 'ICTL01', PROGNAME )
        ICTL01 = IMISS3

        J = ICSIZE( 2 )                               ! SCC=left, FIP=0
        ALLOCATE( ICTL02( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL02', PROGNAME )
        ICTL02 = IMISS3

        J = ICSIZE( 3 )                               ! SCC=all, FIP=0
        ALLOCATE( ICTL03( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL03', PROGNAME )
        ICTL03 = IMISS3
  
        J = ICSIZE( 4 )                               ! SCC=0, FIP=state
        ALLOCATE( ICTL04( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL04', PROGNAME )
        ICTL04 = IMISS3

        J = ICSIZE( 5 )                               ! SCC=left, FIP=state
        ALLOCATE( ICTL05( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL05', PROGNAME )
        ICTL05 = IMISS3
            
        J = ICSIZE( 6 )                               ! SCC=all, FIP=state
        ALLOCATE( ICTL06( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL06', PROGNAME )
        ICTL06 = IMISS3
                        
        J = ICSIZE( 7 )                               ! SCC=0, FIP=all
        ALLOCATE( ICTL07( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL07', PROGNAME )
        ICTL07 = IMISS3
            
        J = ICSIZE( 8 )                               ! SCC=left, FIP=all
        ALLOCATE( ICTL08( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL08', PROGNAME )
        ICTL08 = IMISS3
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( ICTL09( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL09', PROGNAME )
        ICTL09 = IMISS3
            
        J = ICSIZE( 10 )                              ! PLANT=non-blank, SCC=0
        ALLOCATE( ICTL10( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL10', PROGNAME )
        ICTL10 = IMISS3
            
        J = ICSIZE( 38 )                              ! PLANT=non-blank, MACT=all
        ALLOCATE( ICTL38( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL38', PROGNAME )
        ICTL38 = IMISS3
            
        J = ICSIZE( 11 )                              ! PLANT=non-blank, SCC=all     
        ALLOCATE( ICTL11( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL11', PROGNAME )
        ICTL11 = IMISS3
            
        J = ICSIZE( 12 )                              ! CHAR1=non-blank, SCC=all     
        ALLOCATE( ICTL12( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL12', PROGNAME )
        ICTL12 = IMISS3
            
        J = ICSIZE( 13 )                              ! CHAR2=non-blank, SCC=all
        ALLOCATE( ICTL13( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL13', PROGNAME )
        ICTL13 = IMISS3
            
        J = ICSIZE( 14 )                              ! CHAR3=non-blank, SCC=all
        ALLOCATE( ICTL14( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL14', PROGNAME )
        ICTL14 = IMISS3
          
        J = ICSIZE( 15 )                              ! CHAR4=non-blank, SCC=all
        ALLOCATE( ICTL15( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL15', PROGNAME )
        ICTL15 = IMISS3
            
        J = ICSIZE( 16 )                              ! CHAR5=non-blank, SCC=all
        ALLOCATE( ICTL16( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL16', PROGNAME )
        ICTL16 = IMISS3
            
C.........  NOTE- Added later            
        J = MAX( 1, ICSIZE( 17 ) )                     ! SCC=level 1, FIP=0
        ALLOCATE( ICTL02A( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL02A', PROGNAME )
        ICTL02A = IMISS3

        J = MAX( 1, ICSIZE( 18 ) )                     ! SCC=level 2, FIP=0
        ALLOCATE( ICTL02B( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL02B', PROGNAME )
        ICTL02B = IMISS3

        J = MAX( 1, ICSIZE( 19 ) )                     ! SCC=level 3, FIP=0
        ALLOCATE( ICTL02C( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL02C', PROGNAME )
        ICTL02C = IMISS3

        J = MAX( 1, ICSIZE( 20 ) )                     ! SCC=level 1, FIP=state
        ALLOCATE( ICTL05A( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL05A', PROGNAME )
        ICTL05A = IMISS3

        J = MAX( 1, ICSIZE( 21 ) )                     ! SCC=level 2, FIP=state
        ALLOCATE( ICTL05B( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL05B', PROGNAME )
        ICTL05B = IMISS3

        J = MAX( 1, ICSIZE( 22 ) )                     ! SCC=level 3, FIP=state
        ALLOCATE( ICTL05C( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL05C', PROGNAME )
        ICTL05C = IMISS3

        J = MAX( 1, ICSIZE( 23 ) )                     ! SCC=level 1, FIP=all
        ALLOCATE( ICTL08A( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL08A', PROGNAME )
        ICTL08A = IMISS3

        J = MAX( 1, ICSIZE( 24 ) )                     ! SCC=level 2, FIP=all
        ALLOCATE( ICTL08B( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL08B', PROGNAME )
        ICTL08B = IMISS3

        J = MAX( 1, ICSIZE( 25 ) )                     ! SCC=level 3, FIP=all
        ALLOCATE( ICTL08C( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL08C', PROGNAME )
        ICTL08C = IMISS3

        J = MAX( 1, ICSIZE( 26 ) )                     ! SIC=2-digit, FIP=0
        ALLOCATE( ICTL26( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL26', PROGNAME )
        ICTL26 = IMISS3

        J = MAX( 1, ICSIZE( 27 ) )                     ! SIC=all, FIP=0
        ALLOCATE( ICTL27( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL27', PROGNAME )
        ICTL27 = IMISS3

        J = MAX( 1, ICSIZE( 28 ) )                     ! SIC=2-digit, FIP=state
        ALLOCATE( ICTL28( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL28', PROGNAME )
        ICTL28 = IMISS3

        J = MAX( 1, ICSIZE( 29 ) )                     ! SIC=all, FIP=state
        ALLOCATE( ICTL29( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL29', PROGNAME )
        ICTL29 = IMISS3

        J = MAX( 1, ICSIZE( 30 ) )                     ! SIC=2-digit, FIP=county
        ALLOCATE( ICTL30( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL30', PROGNAME )
        ICTL30 = IMISS3

        J = MAX( 1, ICSIZE( 31 ) )                     ! SIC=all, FIP=county
        ALLOCATE( ICTL31( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL31', PROGNAME )
        ICTL31 = IMISS3

!.........  MACT code matches
        J = MAX( 1, ICSIZE( 32 ) )                     ! FIP=0, SCC=0, MACT=all
        ALLOCATE( ICTL32( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL32', PROGNAME )
        ICTL32 = IMISS3

        J = MAX( 1, ICSIZE( 33 ) )                     ! FIP=0, SCC=all, MACT=all
        ALLOCATE( ICTL33( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL33', PROGNAME )
        ICTL33 = IMISS3

        J = MAX( 1, ICSIZE( 34 ) )                     ! FIP=state, SCC=0, MACT=all
        ALLOCATE( ICTL34( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL34', PROGNAME )
        ICTL34 = IMISS3

        J = MAX( 1, ICSIZE( 35 ) )                     ! FIP=state, SCC=all, MACT=all
        ALLOCATE( ICTL35( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL35', PROGNAME )
        ICTL35 = IMISS3

        J = MAX( 1, ICSIZE( 36 ) )                     ! FIP=all, SCC=0, MACT=all
        ALLOCATE( ICTL36( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL36', PROGNAME )
        ICTL36 = IMISS3

        J = MAX( 1, ICSIZE( 37 ) )                     ! FIP=all, SCC=all, MACT=all
        ALLOCATE( ICTL37( -1:J,NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTL37', PROGNAME )
        ICTL37 = IMISS3

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCCTBL
