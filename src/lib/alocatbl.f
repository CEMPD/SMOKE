
        SUBROUTINE ALOCATBL( ICSIZE )

C***********************************************************************
C  subroutine body starts at line 46
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the area-
C      to-point file that contains the table number, row number, and count,
C      and it initializes these to missing.  The subroutine argument is 
C      an array that contains the dimensions for each of the different groups
C      of the cross-reference. Only group 9 (full FIPS, full SCC) is 
C      currently valid for the area-to-point assignments.
C      
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 11/02 by M. Houyoux
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

        USE MODXREF, ONLY: ARPT08, ARPT09

        IMPLICIT NONE

C...........   INCLUDES:
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER(16) :: PROGNAME = 'ALOCATBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCATBL

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( ARPT08 ) ) DEALLOCATE( ARPT08, ARPT09 )
                        
        J = ICSIZE( 8 )                               ! SCC=left, FIP=all
        ALLOCATE( ARPT08( J,3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ARPT08', PROGNAME )
        ARPT08 = IMISS3    ! array
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( ARPT09( J,3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ARPT09', PROGNAME )
        ARPT09 = IMISS3    ! array
                        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCATBL
