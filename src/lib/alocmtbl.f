
        SUBROUTINE ALOCMTBL( ICSIZE )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the VMT mix
C      table that contain the various mobile-source characteristics. The  
C      subroutine argument is an array that contains the dimensions for each 
C      of the different groups of the cross-reference.
C      
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 2/2000 by M. Houyoux
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

        USE MODXREF, ONLY: IMVS01, IMVS02, IMVS03, IMVS04, IMVS05,
     &                     IMVS06, IMVS07, IMVS08, IMVS09, IMVS10,
     &                     IMVS11, IMVS12

        IMPLICIT NONE

C...........   INCLUDES:
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER(16) :: PROGNAME = 'ALOCMTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCMTBL

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( IMVS02 ) ) THEN

            DEALLOCATE( IMVS02, IMVS03, IMVS04, IMVS05 )
            DEALLOCATE( IMVS06, IMVS07, IMVS08, IMVS09 )

        END IF

        IMVS01 = IMISS3

        J = ICSIZE( 2 )                               ! SCC=roadtype, FIP=0
        ALLOCATE( IMVS02( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS02', PROGNAME )
        IMVS02 = IMISS3    ! array

        J = ICSIZE( 3 )                               ! SCC=all, FIP=0
        ALLOCATE( IMVS03( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS03', PROGNAME )
        IMVS03 = IMISS3    ! array
                
        J = ICSIZE( 4 )                               ! SCC=0, FIP=state
        ALLOCATE( IMVS04( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS04', PROGNAME )
        IMVS04 = IMISS3    ! array
            
        J = ICSIZE( 5 )                               ! SCC=roadtype, FIP=state
        ALLOCATE( IMVS05( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS05', PROGNAME )
        IMVS05 = IMISS3    ! array
            
        J = ICSIZE( 6 )                               ! SCC=all, FIP=state
        ALLOCATE( IMVS06( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS06', PROGNAME )
        IMVS06 = IMISS3    ! array
                        
        J = ICSIZE( 7 )                               ! SCC=0, FIP=all
        ALLOCATE( IMVS07( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS07', PROGNAME )
        IMVS07 = IMISS3    ! array
            
        J = ICSIZE( 8 )                               ! SCC=roadtype, FIP=all
        ALLOCATE( IMVS08( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS08', PROGNAME )
        IMVS08 = IMISS3    ! array
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( IMVS09( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS09', PROGNAME )
        IMVS09 = IMISS3    ! array
                        
        J = ICSIZE( 10 )                              ! LINK=non-blank, SCC=0
        ALLOCATE( IMVS10( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS10', PROGNAME )
        IMVS09 = IMISS3    ! array
                        
        J = ICSIZE( 11 )                              ! LINK=non-blank, SCC=rtyp
        ALLOCATE( IMVS11( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS11', PROGNAME )
        IMVS09 = IMISS3    ! array
                        
        J = ICSIZE( 12 )                              ! LINK=non-blank, SCC=all
        ALLOCATE( IMVS12( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS12', PROGNAME )
        IMVS09 = IMISS3    ! array
                        
       RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCMTBL
