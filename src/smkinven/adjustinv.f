
        SUBROUTINE ADJUSTINV( NRAWBP, UDEV, YDEV, CDEV, LDEV )

C**************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine reads the area-to-point file and processes the
C      sources as needed. It also reads and assigns the nonHAP exclusions.

C  PRECONDITIONS REQUIRED:
C      CSOURC, POLVAL, and NPCNT arrays allocated and populated in sorted order
C      NSRC is set.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 11/02 by C. Seppanen
C
C**************************************************************************
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

C...........   MODULES for public variables
        
        IMPLICIT NONE

C...........   INCLUDES

C...........   EXTERNAL FUNCTIONS and their descriptions

C...........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (IN) :: NRAWBP  ! no. raw records by pollutant
        INTEGER , INTENT (IN) :: UDEV    ! unit no. for non-HAP exclusions
        INTEGER , INTENT (IN) :: YDEV    ! unit no. for ar-to-point
        INTEGER , INTENT (IN) :: CDEV    ! SCC descriptions unit no.
        INTEGER , INTENT (IN) :: LDEV    ! log file unit no.

C...........   Local pointers

C...........   Other local variables

        CHARACTER*16 :: PROGNAME = 'ADJUSTINV' ! program name

C***********************************************************************
C   begin body of subroutine ADJUSTINV

C..........  If area-to-point factors file is present...
        IF( YDEV .GT. 0 ) THEN

C.............  Read and preprocess area-to-point factors file
C.............  Result of this call is that the NAR2PT and AR2PTABL 
C               arrays from MODAR2PT and the CHRT09 and ARPT09 arrays
C               from MODLISTS will be populated.
            CALL RDAR2PT( YDEV, CDEV, LDEV )

C.............  Assign area-to-point cross-reference entries to sources
C.............  Result of this call is that the AR2PTTBL, AR2PTIDX, and
C               AR2PTCNT arrays from MODLISTS will be populated
            CALL ASGNAR2PT

C.............  Process area-to-point sources
            CALL PROCAR2PT( NRAWBP )
            
        END IF

C..........  If non-HAP exclusions file is present...
        IF( UDEV .GT. 0 ) THEN

C.............  Read and preprocess NONHAPVOC exclusions x-ref
C.............  Only the CHRT* arrays of the MODXREF will be populated,
C               because we only need to identify the sources, not assign
C               anything to them.
            CALL RDXCLUDE( UDEV )

C.............   Assign array for non-HAP exclusions
            CALL ASGNNHAPX

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE ADJUSTINV
