
        SUBROUTINE ALOCGTBL( ICSIZE )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the gridding
C      cross-reference tables that contain the gridding surrogate code numbers,
C      and it initializes these to missing.  The subroutine argument is 
C      an array that contains the dimensions for each of the different groups
C      of the cross-reference.
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

        USE MODXREF, ONLY: ISRG01, ISRG02, ISRG03, ISRG04, ISRG05,
     &                     ISRG06, ISRG07, ISRG08, ISRG09

        IMPLICIT NONE

C...........   INCLUDES:
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER(16) :: PROGNAME = 'ALOCGTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCGTBL

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( ISRG02 ) ) THEN

            DEALLOCATE( ISRG02, ISRG03, ISRG04, ISRG05 )
            DEALLOCATE( ISRG06, ISRG07, ISRG08, ISRG09 )

        END IF

        ISRG01 = IMISS3

        J = ICSIZE( 2 )                               ! SCC=left, FIP=0
        ALLOCATE( ISRG02( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISRG02', PROGNAME )
        ISRG02 = IMISS3    ! array

        J = ICSIZE( 3 )                               ! SCC=all, FIP=0
        ALLOCATE( ISRG03( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISRG03', PROGNAME )
        ISRG03 = IMISS3    ! array
                
        J = ICSIZE( 4 )                               ! SCC=0, FIP=state
        ALLOCATE( ISRG04( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISRG04', PROGNAME )
        ISRG04 = IMISS3    ! array
            
        J = ICSIZE( 5 )                               ! SCC=left, FIP=state
        ALLOCATE( ISRG05( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISRG05', PROGNAME )
        ISRG05 = IMISS3    ! array
            
        J = ICSIZE( 6 )                               ! SCC=all, FIP=state
        ALLOCATE( ISRG06( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISRG06', PROGNAME )
        ISRG06 = IMISS3    ! array
                        
        J = ICSIZE( 7 )                               ! SCC=0, FIP=all
        ALLOCATE( ISRG07( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISRG07', PROGNAME )
        ISRG07 = IMISS3    ! array
            
        J = ICSIZE( 8 )                               ! SCC=left, FIP=all
        ALLOCATE( ISRG08( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISRG08', PROGNAME )
        ISRG08 = IMISS3    ! array
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( ISRG09( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISRG09', PROGNAME )
        ISRG09 = IMISS3    ! array
                        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCGTBL
