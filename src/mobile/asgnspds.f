
       SUBROUTINE ASGNSPDS
       
C***********************************************************************
C  subroutine body starts at line 107
C
C  DESCRIPTION:
C      For each source, find the most specific match to determine
C      if the source uses a speed profile or uses the inventory speed. 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 12/02 by C. Seppanen (based on asgnnhapx.f)
C
C************************************************************************
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

C...........   MODULES for public variables   
C...........   This module contains the source ararys
        USE MODSOURC, ONLY: CSOURC, CSCC
        
C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: CHRT02, CHRT03, CHRT04, CHRT05, CHRT06,
     &                     CHRT07, CHRT08, CHRT09, ISPD01,
     &                     ISPD02, ISPD03, ISPD04, ISPD05, ISPD06,
     &                     ISPD07, ISPD08, ISPD09, SPDPROFID, TXCNT
        
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, LSCCEND

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         FINDC
        
        EXTERNAL        FINDC
        
C.........  Other local variables
        INTEGER         S                     ! counter
        
        INTEGER          F0, F1, F2, F3, F4, F5  ! tmp find indices
        INTEGER         IOS                   ! I/O status
        
        LOGICAL       :: EFLAG    = .FALSE.
        
        CHARACTER(10)         RWTFMT   ! formt to write rdway type to string
        CHARACTER(10)         VIDFMT   ! format to write veh ID to string
        CHARACTER(STALEN3)    CSTA     ! tmp Country/state code
        CHARACTER(STSLEN3)    CSTASCC  ! tmp Country/state code // SCC
        CHARACTER(STSLEN3)    CSTASL   ! tmp Country/state code // left SCC
        CHARACTER(SCCLEN3)    TSCCL    ! tmp left digits of TSCC
        CHARACTER(SRCLEN3)    CSRC     ! tmp source chars string
        CHARACTER(RWTLEN3)    CRWT     !  buffer for roadway type
        CHARACTER(FIPLEN3)    CFIP     ! tmp (character) FIPS code
        CHARACTER(FPSLEN3)    CFIPSCC  ! tmp FIPS code // SCC
        CHARACTER(FPSLEN3)    CFIPSL   ! tmp FIPS code // left SCC
        CHARACTER(SCCLEN3)    TSCC     ! tmp 10-digit SCC
        CHARACTER(VIDLEN3)    CVID     ! buffer for vehicle type ID
        CHARACTER(300)        MESG     ! message buffer

        CHARACTER(16) :: PROGNAME = 'ASGNSPDS' ! program name

C***********************************************************************
C   begin body of subroutine ASGNSPDS

C.........  Allocate memory for source-based array from MODXREF
        ALLOCATE( SPDPROFID( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPDPROFID', PROGNAME )
        SPDPROFID = -1   ! array

C.........  Set up formats
        WRITE( RWTFMT, '("(I",I2.2,".",I2.2,")")' ) RWTLEN3, RWTLEN3
        WRITE( VIDFMT, '("(I",I2.2,".",I2.2,")")' ) VIDLEN3, VIDLEN3
                
C.........  Loop through the sources
        DO S = 1, NSRC
        
C.............  Create selection
            CSRC    = CSOURC( S )
            TSCC    = CSCC  ( S )
            CFIP    = CSRC  ( 1:FIPLEN3 )
            CSTA    = CFIP  ( 1:STALEN3 )
            
C            WRITE( CRWT, RWTFMT ) IRCLAS( S )
C            WRITE( CVID, VIDFMT ) IVTYPE( S )
C            TSCC = CRWT // CVID
C            CALL PADZERO( TSCC )
            
            TSCCL   = TSCC( 1:LSCCEND )
            CFIPSCC = CFIP // TSCC
            CFIPSL  = CFIP // TSCCL
            CSTASCC = CSTA // TSCC
            CSTASL  = CSTA // TSCCL
            

C.................  Try for FIPS code & SCC match; then
C                           FIPS code & left SCC match; then
C                           Cy/st code & SCC match; then
C                           Cy/st code & left SCC match; then
C                           SCC match; then
C                           left SCC match

            F5 = FINDC( CFIPSCC, TXCNT( 9 ), CHRT09 ) 
            F4 = FINDC( CFIPSL , TXCNT( 8 ), CHRT08 ) 
            F3 = FINDC( CSTASCC, TXCNT( 6 ), CHRT06 ) 
            F2 = FINDC( CSTASL , TXCNT( 5 ), CHRT05 ) 
            F1 = FINDC( TSCC   , TXCNT( 3 ), CHRT03 ) 
            F0 = FINDC( TSCCL  , TXCNT( 2 ), CHRT02 ) 

            IF( F5 .GT. 0 ) THEN
                SPDPROFID( S ) = ISPD09( F5 ) 
                CYCLE                       !  to end of sources-loop

            ELSEIF( F4 .GT. 0 ) THEN
                SPDPROFID( S ) = ISPD08( F4 ) 
                CYCLE                       !  to end of sources-loop

            ELSEIF( F3 .GT. 0 ) THEN
                SPDPROFID( S ) = ISPD06( F3 ) 
                CYCLE                       !  to end of sources-loop

            ELSEIF( F2 .GT. 0 ) THEN
                SPDPROFID( S ) = ISPD05( F2 ) 
                CYCLE                       !  to end of sources-loop

            ELSEIF( F1 .GT. 0 ) THEN
                SPDPROFID( S ) = ISPD03( F1 ) 
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 ) THEN
                SPDPROFID( S ) = ISPD02( F0 ) 
                CYCLE                       !  to end of sources-loop

            END IF

C.............  Try for any FIPS code match
            F0 = FINDC( CFIP, TXCNT( 7 ), CHRT07 ) 

            IF( F0 .GT. 0 ) THEN
                SPDPROFID( S ) = ISPD07( F0 ) 
                CYCLE                       !  to end of sources-loop
            END IF

C.............  Try for any country/state code match (not, pol-specific)
            F0 = FINDC( CSTA, TXCNT( 4 ), CHRT04 ) 

            IF( F0 .GT. 0 ) THEN
                SPDPROFID( S ) = ISPD04( F0 ) 
                CYCLE                       !  to end of sources-loop
            END IF

C.............  Use default if available
            IF( ISPD01 /= IMISS3 ) THEN
                SPDPROFID( S ) = ISPD01
                CYCLE
            END IF
            
        END DO  ! end loop on sources
        
        IF( EFLAG ) THEN
            MESG = 'Problem assigning speed profiles to sources'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        RETURN
        
        END SUBROUTINE ASGNSPDS
        
