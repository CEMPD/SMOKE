
        SUBROUTINE INITSTCY

C***********************************************************************
C  subroutine INITSTCY body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to initialize the necessary fields
C      for performing state and county totals.  The first call sets up the
C      indices from each source to each county.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 8/99 by M. Houyoux
C
C***********************************************************************
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: NMSRC, MCFIP,     ! no. of sources by category
     &                      LREPSTA,          ! output state total emissions flag
     &                      LREPCNY,          ! output county total emissions flag
     &                      LREPSCC,          ! output SCC total emissions flag
     &                      MEBSTA, ! state total speciated emissions
     &                      MEBCNY, ! county total speciated emissions
     &                      MEBSUM, ! source total speciated emissions total
     &                      MEBSRC, ! source total speciated emissions by hour
     &                      MEBSCC, ! SCC total speciated emissions
     &                      MEBSTC  ! state-SCC total speciated emissions
     
C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: MICNY, NCOUNTY, CNTYCOD

C...........   This module is the source inventory arrays
        USE MODSOURC, ONLY: CIFIP, CSCC

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: MISCC, DISFLTYP, DISFL, GASFLTYP, GASFL, ETHFLTYP, ETHFL,
     &                       CNGFLTYP, CNGFL, LPGFLTYP, LPGFL
C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVSCC, INVSCC, NINVIFIP, INVCFIP

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1  
        INTEGER         FINDC, ENVINT

        EXTERNAL   FIND1, FINDC, ENVINT

C...........   Other local variables

        INTEGER          IOS      ! i/o status
        INTEGER          J, K, S  ! counter
        INTEGER          IDISFL    ! diesel fuel type
        INTEGER          IGASFL    ! gasoline fuel type
        INTEGER          IETHFL    ! ethanol fuel type
        INTEGER          ICNGFL    ! CNG fuel type
        INTEGER          ILPGFL    ! LPG fuel type
 
        LOGICAL, SAVE :: FIRSTIME = .TRUE. ! true: first time routine called

        CHARACTER(FIPLEN3) CFIP      ! tmp cy/st/co code
        CHARACTER(FIPLEN3) PFIP      ! previous cy/st/co code

        CHARACTER(300)   MESG     ! message buffer

        CHARACTER(16) :: PROGNAME = 'INITSTCY' ! program name

C***********************************************************************
C   begin body of subroutine INITSTCY
        
        IF( FIRSTIME ) THEN

C.............  Allocate memory for indices from Co/st/cy codes to counties
            ALLOCATE( MICNY( NMSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MICNY', PROGNAME )

C.............  Allocate memory for index from master list of SCCs to source SCC
            ALLOCATE( MISCC( NMSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MISCC', PROGNAME )
            ALLOCATE( MCFIP( NMSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MCFIP', PROGNAME )
            ALLOCATE( GASFL( NMSRC ), STAT=IOS )       ! arry for gasoline fuel type or not
            CALL CHECKMEM( IOS, 'GASFL', PROGNAME )
            ALLOCATE( DISFL( NMSRC ), STAT=IOS )       ! arry for diesel fuel type or not
            CALL CHECKMEM( IOS, 'DISFL', PROGNAME )
            ALLOCATE( CNGFL( NMSRC ), STAT=IOS )       ! arry for CNG fuel type or not
            CALL CHECKMEM( IOS, 'CNGFL', PROGNAME )
            ALLOCATE( LPGFL( NMSRC ), STAT=IOS )       ! arry for LPG fuel type or not
            CALL CHECKMEM( IOS, 'LPGFL', PROGNAME )
            ALLOCATE( ETHFL( NMSRC ), STAT=IOS )       ! arry for ethnol fuel type or not
            CALL CHECKMEM( IOS, 'ETHFL', PROGNAME )
            DISFL = .FALSE.
            GASFL = .FALSE.
            CNGFL = .FALSE.
            LPGFL = .FALSE.
            ETHFL = .FALSE.

            IGASFL = ENVINT( 'GASOLINE_FUEL_CODE', 'Gasoline fuel ' //
     &                      'type code [ex: 1]', 1, IOS )
            WRITE( GASFLTYP,'( I2.2)' ) IGASFL
 
            IDISFL = ENVINT( 'DIESEL_FUEL_CODE', 'Diesel fuel ' //
     &                      'type code [ex: 2]', 2, IOS )
            WRITE( DISFLTYP,'( I2.2)' ) IDISFL

            ICNGFL = ENVINT( 'CNG_FUEL_CODE', 'CNG fuel ' //
     &                      'type code [ex: 3]', 3, IOS )
            WRITE( CNGFLTYP,'( I2.2)' ) ICNGFL

            ILPGFL = ENVINT( 'LPG_FUEL_CODE', 'LPG fuel ' //
     &                      'type code [ex: 4]', 4, IOS )
            WRITE( LPGFLTYP,'( I2.2)' ) ILPGFL

            IETHFL = ENVINT( 'ETHANOL_FUEL_CODE', 'Ethanol fuel ' //
     &                      'type code [ex: 5]', 5, IOS )
            WRITE( ETHFLTYP,'( I2.2)' ) IETHFL

C.............  Create indices to counties from Co/st/cy codes and for SCCs
            PFIP = ' ' 
            DO S = 1, NMSRC
            
                CFIP = CIFIP( S )
                
                IF( CFIP .NE. PFIP ) THEN
                
                    J = MAX( FINDC( CFIP, NCOUNTY, CNTYCOD ), 0 )
                    K = MAX( FINDC( CFIP, NINVIFIP,INVCFIP ), 0 )
                    PFIP = CFIP
                    
                END IF
                
                MICNY( S ) = J
                MCFIP( S ) = K 
                MISCC( S ) = MAX( FINDC( CSCC( S ), NINVSCC, INVSCC ), 0 )

C..............   Check whether diseel fuel type or not
                IF( CSCC(S)( 11:12 ) == '22' ) THEN
                    IF( CSCC(S)( 13:14 ) == GASFLTYP ) GASFL( S ) = .TRUE.
                    IF( CSCC(S)( 13:14 ) == DISFLTYP ) DISFL( S ) = .TRUE.
                    IF( CSCC(S)( 13:14 ) == CNGFLTYP ) CNGFL( S ) = .TRUE.
                    IF( CSCC(S)( 13:14 ) == LPGFLTYP ) LPGFL( S ) = .TRUE.
                    IF( CSCC(S)( 13:14 ) == ETHFLTYP ) ETHFL( S ) = .TRUE.
                END IF
 
            END DO
            
            FIRSTIME = .FALSE.
        
        END IF

C.........  Initialize totals to zero...
C.........  SCC totals...
        IF( LREPSCC ) THEN
            MEBSCC = 0.
        END IF

C.........  State totals...
        IF( LREPSTA ) THEN
            MEBSTA = 0.
            IF( LREPSCC ) THEN
                MEBSTC = 0.
            END IF
        END IF

C.........  County totals...
        IF( LREPCNY ) THEN
            MEBCNY = 0.
        END IF

C.........  Source totals...
        MEBSUM = 0.

        RETURN

        END SUBROUTINE INITSTCY
