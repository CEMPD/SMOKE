
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE EFCOMBO( FNAME1, FNAME2, NTEMP, MPTYPE, 
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

C.........  Subroutine call variables

        REAL          PERCENT( MPTYPE )  ! percentage of each PSI in combo

        INTEGER       JDATE              ! date ID for reading EF files
        INTEGER       LIST   ( MPTYPE )  ! list of PSIs in combo
        INTEGER       MPTYPE             ! number of factors in combination
        INTEGER       NTEMP              ! non-diurnal number of temperatures

        CHARACTER*16  FNAME1             ! File name for non-diurnal EF file
        CHARACTER*16  FNAME2             ! File name for diurnal EF file

C.........  Read-in emission factors storage
   
        REAL    NDEF( MAXTMP , NVTYPE, MXNDET )  ! Non-diurnal Emis fac storage
        REAL    DEF ( NMINMAX, NVTYPE, MXDET  )  ! Diurnal     Emis fac storage

C...........   Master Non-Diurnal EFs list variables

        REAL    EXH_CO ( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: CO exhaust
        REAL    EXH_NOX( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: NOX exhaust
        REAL    EXH_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC exhaust
        REAL    EVP_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC evaporative
        REAL    RNL_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC running
        REAL    RST_VOC( MAXTMP, NVTYPE )  ! Non-Diurn Emis fac: VOC resting
 
        COMMON / EFS1 / EXH_CO , EXH_NOX, EXH_VOC, EVP_VOC, RNL_VOC,
     &                  RST_VOC
 
C...........   Master Diurnal EFs list variables
 
        REAL    WDL_VOC( MXTMMI, NVTYPE )  ! Diurnal Emis fac: VOC weighted
        REAL    DNL_VOC( MXTMMI, NVTYPE )  ! Diurnal Emis fac: VOC diurnal
        REAL    HOT_VOC( MXTMMI, NVTYPE )  ! Diurnal Emis fac: VOC hot soak
        REAL    CRC_VOC( MXTMMI, NVTYPE )  ! Diurnal Emis fac: VOC crank case

        COMMON / EFS2 / WDL_VOC, DNL_VOC, HOT_VOC, CRC_VOC

        INTEGER I, V, T      ! Counters and pointers
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

   11       CONTINUE 

            DO 13 T= 1, NMINMAX

                WDL_VOC( T,V ) = 0.
                DNL_VOC( T,V ) = 0.
                HOT_VOC( T,V ) = 0.
                CRC_VOC( T,V ) = 0.

   13       CONTINUE 

   22   CONTINUE 

C.........  Loop over components of combo EFs
        DO 55 I = 1, MPTYPE

C.............  Read in each required existing Non-diurnal EF set
            IF( .NOT. READ3
     &          ( FNAME1, 'ALL', 1, JDATE, LIST( I ), NDEF ) ) THEN

                MESG = 'Could not read non-diurnal EFs from ' // FNAME1
                CALL M3EXIT( 'EFCOMBO', JDATE, LIST( I ), MESG, 2 )

            ENDIF 

C.............  Read in each required existing Diurnal EF set
            IF( .NOT. READ3
     &          ( FNAME2, 'ALL', 1, JDATE, LIST( I ), DEF ) ) THEN

                MESG = 'Could not read diurnal EFs from ' // FNAME2
                CALL M3EXIT( 'EFCOMBO', JDATE, LIST( I ), MESG, 2 )

            ENDIF 

C.............  Loop over vehicle classes
            DO 44 V = 1, NVTYPE
 
C.................  Loop over non-diurnal temperature range
                DO 33 T = 1, NTEMP
                    
                    EXH_CO ( T,V ) = EXH_CO ( T,V ) + 
     &                               NDEF( T,V,1 ) * PERCENT( I )
                   
                    EXH_NOX( T,V ) = EXH_NOX( T,V ) + 
     &                               NDEF( T,V,2 ) * PERCENT( I )
                  
                    EXH_VOC( T,V ) = EXH_VOC( T,V ) + 
     &                               NDEF( T,V,3 ) * PERCENT( I )
                 
                    EVP_VOC( T,V ) = EVP_VOC( T,V ) + 
     &                               NDEF( T,V,4 ) * PERCENT( I )

                    RNL_VOC( T,V ) = RNL_VOC( T,V ) +
     &                               NDEF( T,V,5 ) * PERCENT( I )

                    RST_VOC( T,V ) = RST_VOC( T,V ) + 
     &                               NDEF( T,V,6 ) * PERCENT( I )
              
   33           CONTINUE  ! Non-diurnal Temperatures

C.................  Loop over diurnal temperature range

                DO 35 T = 1, NMINMAX

                    WDL_VOC( T,V ) = WDL_VOC( T,V ) + 
     &                               DEF( T,V,1 ) * PERCENT( I )
           
                    DNL_VOC( T,V ) = DNL_VOC( T,V ) + 
     &                               DEF( T,V,2 ) * PERCENT( I )
           
                    HOT_VOC( T,V ) = HOT_VOC( T,V ) + 
     &                               DEF( T,V,3 ) * PERCENT( I )
           
                    CRC_VOC( T,V ) = CRC_VOC( T,V ) + 
     &                               DEF( T,V,4 ) * PERCENT( I )
          
   35           CONTINUE  ! Diurnal Temperatures
 
   44       CONTINUE      ! Vehicle type
 
   55   CONTINUE

        RETURN 
        END
