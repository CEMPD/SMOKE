
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE EFCOMBO( FNAME, NTEMP, MPTYPE, 
     &                      PERCENT, LIST, JDATE )

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C      This subroutine creates emission factors using multiple previously 
C      calculated emission factors and percentages for each
C
C  PRECONDITIONS REQUIRED:
C      Must already have pure emission factors from which can calculate
C      combination factors.  Must have arrays that specify how many and
C      which pure EFs and percentages for each.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 2/96 by MRH
C
C***********************************************************************

        IMPLICIT NONE

        INCLUDE 'CHDIMS3.EXT'
        INCLUDE 'MBDIMS3.EXT'
        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures.

        INTEGER JDATE, MPTYPE, NTEMP
        REAL    PERCENT( MPTYPE )
        INTEGER LIST   ( MPTYPE )
        CHARACTER*(*) FNAME

C.........  Read-in emission factors storage
   
        REAL    EF( MAXTMP, NVTYPE, 11 )  ! Emis factor storage

C.........  Output emission factors
        REAL    EXH_CO ( MAXTMP, NVTYPE )  ! Emis factor: CO exhaust
        REAL    EXH_NOX( MAXTMP, NVTYPE )  ! Emis factor: NOX exhaust
        REAL    EXH_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC exhaust
        REAL    EVP_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC evaporative
        REAL    RNL_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC running
        REAL    RST_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC resting
        REAL    RF1_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC refuel /mile
        REAL    RF2_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC refuel /gallon
        REAL    HOT_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC hot soak
        REAL    CRC_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC crank case
        REAL    RHR_VOC( MAXTMP, NVTYPE )  ! Emis factor: VOC resting /hour

        COMMON / EFS / EXH_CO , EXH_NOX, EXH_VOC, EVP_VOC, RNL_VOC,
     &                 RST_VOC, RF1_VOC, RF2_VOC, HOT_VOC, CRC_VOC,
     &                 RHR_VOC

        INTEGER I, V, T
        CHARACTER*256 MESG

C.........  Initialize EF sum variables
        DO 22 V = 1, NVTYPE
            DO 11 T= 1, NTEMP

                    EXH_CO ( T,V ) = 0.
                    EXH_NOX( T,V ) = 0.
                    EXH_VOC( T,V ) = 0.
                    EVP_VOC( T,V ) = 0.
                    RNL_VOC( T,V ) = 0.
                    RST_VOC( T,V ) = 0.
                    RF1_VOC( T,V ) = 0.
                    RF2_VOC( T,V ) = 0.
                    HOT_VOC( T,V ) = 0.
                    CRC_VOC( T,V ) = 0.
                    RHR_VOC( T,V ) = 0.

   11       CONTINUE 
   22   CONTINUE 

C.........  Loop over components of combo EFs
        DO 55 I = 1, MPTYPE

C.............  Read in each required existing EF set
            IF( .NOT. READ3
     &          ( FNAME, 'ALL', 1, JDATE, LIST( I ), EF ) ) THEN

                MESG = 'Could not read EFs from ' // FNAME
                CALL M3EXIT( 'EFCOMBO', 0, 0, MESG, 2 )

            ENDIF 

C.............  Loop over vehicle classes
            DO 44 V = 1, NVTYPE
 
C.................  Loop over temperature range
                DO 33 T = 1, NTEMP
                    
                    EXH_CO ( T,V ) = EXH_CO ( T,V ) + 
     &                               EF( T,V,1 ) * PERCENT( I )
                   
                    EXH_NOX( T,V ) = EXH_NOX( T,V ) + 
     &                               EF( T,V,2 ) * PERCENT( I )
                  
                    EXH_VOC( T,V ) = EXH_VOC( T,V ) + 
     &                               EF( T,V,3 ) * PERCENT( I )
                 
                    EVP_VOC( T,V ) = EVP_VOC( T,V ) + 
     &                               EF( T,V,4 ) * PERCENT( I )

                    RNL_VOC( T,V ) = RNL_VOC( T,V ) +
     &                               EF( T,V,5 ) * PERCENT( I )

                    RST_VOC( T,V ) = RST_VOC( T,V ) + 
     &                               EF( T,V,6 ) * PERCENT( I )
              
                    RF1_VOC( T,V ) = RF1_VOC( T,V ) + 
     &                               EF( T,V,7 ) * PERCENT( I )
             
                    RF2_VOC( T,V ) = RF2_VOC( T,V ) + 
     &                               EF( T,V,8 ) * PERCENT( I )
            
                    HOT_VOC( T,V ) = HOT_VOC( T,V ) + 
     &                               EF( T,V,9 ) * PERCENT( I )
           
                    CRC_VOC( T,V ) = CRC_VOC( T,V ) + 
     &                               EF( T,V,10 ) * PERCENT( I )
          
                    RHR_VOC( T,V ) = RHR_VOC( T,V ) + 
     &                               EF( T,V,11 ) * PERCENT( I )
 
   33           CONTINUE  ! Temperature
 
   44       CONTINUE      ! Vehicle type
 
   55   CONTINUE

        RETURN 
        END
