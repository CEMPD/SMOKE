
       SUBROUTINE ALOCPTBL( ICSIZE )
       
C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the speed
C      cross-reference tables that contain the speed profile numbers, and
C      it initializes these to missing.  The subroutine argument is
C      an array that contains the dimensions for each of the different 
C      groups of the cross-reference.
C      
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 1/03 by C. Seppanen
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
        USE MODXREF, ONLY: ISPD01, ISPD02, ISPD03, ISPD04, ISPD05,
     &                     ISPD06, ISPD07, ISPD08, ISPD09

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER*16 :: PROGNAME = 'ALOCPTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCPTBL

C.........  First deallocate if these have previously been allocated
        IF( ALLOCATED( ISPD02 ) ) THEN

            DEALLOCATE( ISPD02, ISPD03, ISPD04, ISPD05 )
            DEALLOCATE( ISPD06, ISPD07, ISPD08, ISPD09 )

        END IF
        
        ISPD01 = IMISS3

        J = ICSIZE( 2 )                               ! SCC=left, FIP=0
        ALLOCATE( ISPD02( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPD02', PROGNAME )
        ISPD02 = IMISS3    ! array

        J = ICSIZE( 3 )                               ! SCC=all, FIP=0
        ALLOCATE( ISPD03( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPD03', PROGNAME )
        ISPD03 = IMISS3    ! array
                
        J = ICSIZE( 4 )                               ! SCC=0, FIP=state
        ALLOCATE( ISPD04( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPD04', PROGNAME )
        ISPD04 = IMISS3    ! array
            
        J = ICSIZE( 5 )                               ! SCC=left, FIP=state
        ALLOCATE( ISPD05( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPD05', PROGNAME )
        ISPD05 = IMISS3    ! array
            
        J = ICSIZE( 6 )                               ! SCC=all, FIP=state
        ALLOCATE( ISPD06( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPD06', PROGNAME )
        ISPD06 = IMISS3    ! array
                        
        J = ICSIZE( 7 )                               ! SCC=0, FIP=all
        ALLOCATE( ISPD07( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPD07', PROGNAME )
        ISPD07 = IMISS3    ! array
            
        J = ICSIZE( 8 )                               ! SCC=left, FIP=all
        ALLOCATE( ISPD08( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPD08', PROGNAME )
        ISPD08 = IMISS3    ! array
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( ISPD09( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPD09', PROGNAME )
        ISPD09 = IMISS3    ! array
                        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCPTBL
