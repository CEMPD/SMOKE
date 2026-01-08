
        SUBROUTINE ALOCETBL( NACTV, ICSIZE )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the emission factor 
C      cross-reference tables that contain the indices to the unsorted PSIs, and
C      it initializes these to missing.  The subroutine arguments are the number
C      of inventory pollutants and an array that contains the dimensions for 
C      each of the different groups of the cross-reference.
C      
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 6/99 by M. Houyoux
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

        USE MODXREF, ONLY: IEFS01, IEFS02, IEFS03, IEFS04, IEFS05, 
     &                     IEFS06, IEFS07, IEFS08, IEFS09, IEFS10,
     &                     IEFS11, IEFS12, IEFS13, IEFS14, IEFS15, 
     &                     IEFS16

        IMPLICIT NONE

C...........   INCLUDES
C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: NACTV        ! number of pollutants
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       J     ! counter and indices

        INTEGER       IOS              ! i/o status

        CHARACTER(16) :: PROGNAME = 'ALOCETBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCETBL

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( IEFS02 ) ) THEN

            DEALLOCATE( IEFS02, IEFS03, IEFS04, IEFS05, IEFS06 )
            DEALLOCATE( IEFS07, IEFS08, IEFS09, IEFS10, IEFS11 )
            DEALLOCATE( IEFS12, IEFS13, IEFS14, IEFS15, IEFS16 )

        END IF

        ALLOCATE( IEFS01( NACTV ), STAT=IOS )         ! SCC=0, FIP=0
        CALL CHECKMEM( IOS, 'IEFS01', PROGNAME )
        IEFS01 = IMISS3

        J = ICSIZE( 2 )                               ! SCC=left, FIP=0
        ALLOCATE( IEFS02( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS02', PROGNAME )
        IEFS02 = IMISS3

        J = ICSIZE( 3 )                               ! SCC=all, FIP=0
        ALLOCATE( IEFS03( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS03', PROGNAME )
        IEFS03 = IMISS3
  
        J = ICSIZE( 4 )                               ! SCC=0, FIP=state
        ALLOCATE( IEFS04( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS04', PROGNAME )
        IEFS04 = IMISS3

        J = ICSIZE( 5 )                               ! SCC=left, FIP=state
        ALLOCATE( IEFS05( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS05', PROGNAME )
        IEFS05 = IMISS3
            
        J = ICSIZE( 6 )                               ! SCC=all, FIP=state
        ALLOCATE( IEFS06( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS06', PROGNAME )
        IEFS06 = IMISS3
                        
        J = ICSIZE( 7 )                               ! SCC=0, FIP=all
        ALLOCATE( IEFS07( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS07', PROGNAME )
        IEFS07 = IMISS3
            
        J = ICSIZE( 8 )                               ! SCC=left, FIP=all
        ALLOCATE( IEFS08( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS08', PROGNAME )
        IEFS08 = IMISS3
                        
        J = ICSIZE( 9 )                               ! SCC=all, FIP=all
        ALLOCATE( IEFS09( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS09', PROGNAME )
        IEFS09 = IMISS3
            
        J = ICSIZE( 10 )                              ! PLANT=non-blank, SCC=0
        ALLOCATE( IEFS10( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS10', PROGNAME )
        IEFS10 = IMISS3
            
        J = ICSIZE( 11 )                              ! PLANT=non-blank, SCC=all     
        ALLOCATE( IEFS11( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS11', PROGNAME )
        IEFS11 = IMISS3
            
        J = ICSIZE( 12 )                              ! CHAR1=non-blank, SCC=all     
        ALLOCATE( IEFS12( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS12', PROGNAME )
        IEFS12 = IMISS3
            
        J = ICSIZE( 13 )                              ! CHAR2=non-blank, SCC=all
        ALLOCATE( IEFS13( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS13', PROGNAME )
        IEFS13 = IMISS3
            
        J = ICSIZE( 14 )                              ! CHAR3=non-blank, SCC=all
        ALLOCATE( IEFS14( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS14', PROGNAME )
        IEFS14 = IMISS3
          
        J = ICSIZE( 15 )                              ! CHAR4=non-blank, SCC=all
        ALLOCATE( IEFS15( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS15', PROGNAME )
        IEFS15 = IMISS3
            
        J = ICSIZE( 16 )                              ! CHAR5=non-blank, SCC=all
        ALLOCATE( IEFS16( J,NACTV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IEFS16', PROGNAME )
        IEFS16 = IMISS3
            
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALOCETBL
