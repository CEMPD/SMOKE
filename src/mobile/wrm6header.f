
        SUBROUTINE WRM6HEADER( MDEV )

C***********************************************************************
C  subroutine body starts at line 71
C
C  DESCRIPTION:
C       Writes MOBILE6 header commands to the input file
C
C  PRECONDITIONS REQUIRED:
C       MDEV has been opened
C       MEPROC pollutant names must match those in this file
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
C
C***********************************************************************
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
C***********************************************************************

C...........   This module contains emission factor tables and related
        USE MODEMFAC, ONLY: EMTPOL, NEPOL, OUTPUTHC
        
        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        
C...........   EXTERNAL FUNCTIONS
        CHARACTER*2, EXTERNAL :: CRLF
                
C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: MDEV    ! M6 input file unit no.        

C...........   Local variables
        INTEGER           I             ! counter

        LOGICAL           WROTE_OCARB   ! true: wrote OCARBON to pol array
        LOGICAL           WROTE_ECARB   ! true: wrote ECARBON to pol array
        LOGICAL           WROTE_GASPM   ! true: wrote GASPM to pol array
        LOGICAL           WROTE_BRAKE   ! true: wrote BRAKE to pol array
        LOGICAL           WROTE_TIRE    ! true: wrote TIRE to pol array
        LOGICAL           ISUSERHAP     ! true: at least one user-defined HAP 

        CHARACTER(LEN=50)  BASICPOL     ! basic pollutants to calculate
        CHARACTER(LEN=50)  PMPOL        ! particulates to calculate
        CHARACTER(LEN=50)  AIRPOL       ! air toxics to calculate
        CHARACTER(LEN=300) MESG         ! message buffer 
        
        CHARACTER*16 :: PROGNAME = 'WRM6HEADER'   ! program name

C***********************************************************************
C   begin body of subroutine WRM6HEADER

C.........  Initialize arrays and logical variables
        BASICPOL = ' '
        PMPOL    = ' '
        AIRPOL   = ' '
        
        WROTE_OCARB = .FALSE.
        WROTE_ECARB = .FALSE.
        WROTE_GASPM = .FALSE.
        WROTE_BRAKE = .FALSE.
        WROTE_TIRE  = .FALSE.
        ISUSERHAP   = .FALSE.

C.........  Write basic MOBILE6 header info
        WRITE( MDEV,93000 ) 'MOBILE6 INPUT FILE :'
        WRITE( MDEV,93000 ) 'NO DESC OUTPUT     :'
        WRITE( MDEV,93000 ) 'DATABASE OUTPUT    :'

C.........  Loop through pollutants and generate MOBILE6 inputs        
        DO I = 1, NEPOL
            SELECT CASE ( EMTPOL( I ) )
            CASE( 'CO' )
                BASICPOL = TRIM( BASICPOL ) // ' CO'
            CASE( 'NOX' )
                BASICPOL = TRIM( BASICPOL ) // ' NOX'
            CASE( 'VOC', 'THC', 'NMHC', 'TOG', 'NMOG' )
                BASICPOL = TRIM( BASICPOL ) // ' HC'
            CASE( 'SO4' )
                PMPOL = TRIM( PMPOL ) // ' SO4'
            CASE( 'OCARB25', 'OCARBPMC' )
                IF( .NOT. WROTE_OCARB ) THEN
                    PMPOL = TRIM( PMPOL ) // ' OCARBON'
                    WROTE_OCARB = .TRUE.
                END IF
            CASE( 'ECARB25', 'ECARBPMC' )
                IF( .NOT. WROTE_ECARB ) THEN
                    PMPOL = TRIM( PMPOL ) // ' ECARBON'
                    WROTE_ECARB = .TRUE.
                END IF
            CASE( 'GASPM25', 'GASPMC' )
                IF( .NOT. WROTE_GASPM ) THEN
                    PMPOL = TRIM( PMPOL ) // ' GASPM'
                    WROTE_GASPM = .TRUE.
                END IF
            CASE( 'SO2' )
                PMPOL = TRIM( PMPOL ) // ' SO2'
            CASE( 'NH3' )
                PMPOL = TRIM( PMPOL ) // ' NH3'
            CASE( 'BRAKE25', 'BRAKEPMC' )
                IF( .NOT. WROTE_BRAKE ) THEN
                    PMPOL = TRIM( PMPOL ) // ' BRAKE'
                    WROTE_BRAKE = .TRUE.
                END IF
            CASE( 'TIRE25', 'TIREPMC' )
                IF( .NOT. WROTE_TIRE ) THEN
                    PMPOL = TRIM( PMPOL ) // ' TIRE'
                    WROTE_TIRE = .TRUE.
                END IF
            CASE( 'BENZENE' )
                AIRPOL = TRIM( AIRPOL ) // ' BENZ'
            CASE( 'MTBE' )
                AIRPOL = TRIM( AIRPOL ) // ' MTBE'
            CASE( 'BUTADIENE' )
                AIRPOL = TRIM( AIRPOL ) // ' BUTA'
            CASE( 'FORMALDEHYD' )
                AIRPOL = TRIM( AIRPOL ) // ' FORM'
            CASE( 'ACETALD' )
                AIRPOL = TRIM( AIRPOL ) // ' ACET'
            CASE( 'ACROLEIN' )
                AIRPOL = TRIM( AIRPOL ) // ' ACRO'
            CASE DEFAULT
C.................  Check if current pollutant is output hydrocarbon
                IF( EMTPOL( I ) == OUTPUTHC ) THEN
                    BASICPOL = TRIM( BASICPOL ) // ' HC'
                ELSE
                    MESG = 'WARNING: Pollutant ' // 
     &                     TRIM( EMTPOL( I ) ) // ' in MEPROC ' //
     &                     'file is not a MOBILE6 intrisic ' //
     &                     'pollutant.' // CRLF() // BLANK10 // 
     &                     'Assuming that ' // TRIM( EMTPOL( I ) ) // 
     &                     ' is a user-defined HAP.'
                    CALL M3MSG2( MESG )
                    ISUSERHAP = .TRUE.
                END IF
            END SELECT
        END DO

C.........  Exit if trying to generate user-defined HAPs without intrinsic toxics
        IF( ISUSERHAP .AND. AIRPOL == ' ' ) THEN
            MESG = 'ERROR: Cannot run MOBILE6 with user-defined ' //
     &             'HAPs and no intrinsic toxic pollutants.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( BASICPOL /= ' ' ) THEN
            WRITE( MDEV,93000 ) 'POLLUTANTS         : ' // 
     &                          TRIM( BASICPOL )
        END IF
        
        IF( PMPOL /= ' ' ) THEN
            WRITE( MDEV,93000 ) 'PARTICULATES       : ' // 
     &                          TRIM( PMPOL )
        END IF
        
        IF( AIRPOL /= ' ' ) THEN
            WRITE( MDEV,93000 ) 'AIR TOXICS         : ' // 
     &                          TRIM( AIRPOL )
        END IF

C.........  Finish MOBILE6 header
        WRITE( MDEV,93000 ) 'RUN DATA           :'
        WRITE( MDEV,93000 ) ' '

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
        
        END SUBROUTINE WRM6HEADER
