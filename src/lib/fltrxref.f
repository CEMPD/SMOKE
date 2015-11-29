
        SUBROUTINE FLTRXREF( CFIP, CSIC, TSCC, CPOA, CMCT, IXSIC, IXSCC,
     &                       IXPOA, SKIPPOA, SKIPREC )

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      The purpose of subroutine FLTRXREF is be to post-process x-ref 
C      information to scan for '-9', pad with zeros, compare SCC to master list,
C      compare SIC to master list, and compare pollutant/activity name with
C      master list.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     
C
C  REVISION  HISTORY:
C      Started 3/99 by M. Houyoux
C
C**************************************************************************
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

C.........  MODULES for public variables
C...........   This module contains the source arrays
        USE MODSOURC, ONLY: CISIC, CMACT

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: INVSCC, INVSCL, INVSIC, INVSIC2, INVMACT,
     &                      NINVSCC, NINVSCL, NINVSIC, NINVSIC2,
     &                      NINVMACT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, LSCCEND, RSCCBEG, NIPPA, EANAM

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        INTEGER       FIND1
        INTEGER       FINDC
        INTEGER       INDEX1
        INTEGER       STR2INT
        LOGICAL       SETSCCTYPE, CHKEXPSCC, CHKEXPSIC

        EXTERNAL      FIND1, FINDC, INDEX1, STR2INT, SETSCCTYPE

C...........   SUBROUTINE ARGUMENTS:
        CHARACTER(*), INTENT(IN OUT) :: CFIP   ! cntry/state/county code
        CHARACTER(*), INTENT(IN OUT) :: CSIC   ! standard indust. code
        CHARACTER(*), INTENT(IN OUT) :: TSCC   ! source category code
        CHARACTER(*), INTENT(IN OUT) :: CPOA   ! pollutant/activity name
        CHARACTER(*), INTENT(IN OUT) :: CMCT   ! MACT code
        INTEGER     , INTENT   (OUT) :: IXSIC  ! index of SIC in SIC list
        INTEGER     , INTENT   (OUT) :: IXSCC  ! index of SCC in SCC list
        INTEGER     , INTENT   (OUT) :: IXPOA  ! index of pol/act in master list
        LOGICAL     , INTENT   (OUT) :: SKIPPOA! true: skipped rec is pollutant
                                               ! or activity specific
        LOGICAL     , INTENT   (OUT) :: SKIPREC! true: skip record in caller

C...........   Other local variables

        INTEGER          IXACT   ! index to master activity names array
        INTEGER          L       ! tmp string length

        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: 1st time subroutine called 
        LOGICAL, SAVE :: PFLAG    = .FALSE.  ! true: point or area sources
        LOGICAL, SAVE :: MFLAG    = .FALSE.  ! true: have MACT codes
        LOGICAL          SCCFLAG             ! true: SCC type is different from previous

        CHARACTER(300)     MESG         ! message buffer
        CHARACTER(SCCLEN3) SCCL         ! left part of SCC
        CHARACTER(SCCLEN3) SCCR         ! right part of SCC
        CHARACTER(SICLEN3) CSIC2        ! 2-digit SIC
        CHARACTER(SCCLEN3), SAVE :: SCRZERO  ! zero right digits of TSCC
        CHARACTER(SCCLEN3), SAVE :: SCCZERO  ! zero SCC
        CHARACTER(SICLEN3), SAVE :: SICZERO  ! zero SIC

        CHARACTER(16) :: PROGNAME = 'FLTRXREF' ! program name

C***********************************************************************
C   Begin body of subroutine FLTRXREF

C.........  For the first time the subroutine is called...
        IF( FIRSTIME ) THEN

C.............  Make sure unique list of SCC and left-SCCs for this inventory
C               was already generated.
            IF( .NOT. ALLOCATED( INVSCC ) ) THEN
                MESG = 'INTERNAL ERROR: Routine GENUSLST needs to be '//
     &                 'called before the first call to ' // PROGNAME
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF

C.............  Set flags indicating which source category is being processed
            PFLAG = ASSOCIATED( CISIC )
            MFLAG = ASSOCIATED( CMACT )

C.............  Check length of SCC string
            L = LEN( TSCC )
            IF( L .NE. SCCLEN3 ) THEN
                MESG = 'INTERNAL ERROR: Length of SCC string not '//
     &                 'correct in call to ' // PROGNAME
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END IF

C.............  Set up zero strings for FIPS code of zero and SCC code of zero
            SCCZERO = REPEAT( '0', SCCLEN3 )
            SCRZERO = REPEAT( '0', SCCLEN3 - LSCCEND )
            SICZERO = REPEAT( '0', SICLEN3 )

            FIRSTIME = .FALSE.

        END IF

C.........  Initialize call subroutine outputs
        IXSIC = 0
        IXSCC = 0
        IXPOA = 0
        SKIPREC = .FALSE.
        SKIPPOA = .FALSE.

C.........  Smart interpretation of country/state/county code
        CALL FLTRNEG( CFIP )     ! Filter 0 and -9 to blank
        CALL PADZERO( CFIP )     ! Pad with zeros

C.........  Smart interpretation of pollutant/activity name
        CALL FLTRNEG( CPOA )     ! Filter 0 and -9 to blank

C.........  Smart interpretation of SIC
        CALL FLTRNEG( CSIC )     ! Filter 0 and -9 to blank
        CALL PADZERO( CSIC )     ! Pad with zeros

C.........  Smart interpretation of SCC
        CALL FLTRNEG( TSCC )     ! Filter 0 and -9 to blank
        CALL PADZERO( TSCC )     ! Pad with zeros

C.........  Set type of SCC                
        SCCFLAG = SETSCCTYPE( TSCC )
        
C.........  If SCC type has changed, reset zero string for right SCC
        IF( SCCFLAG ) SCRZERO = REPEAT( '0', SCCLEN3 - LSCCEND )

C.........  Smart interpretation of MACT code
        CALL FLTRNEG( CMCT )     ! Filter 0 and -9 to blank
        CALL PADZERO( CMCT )     ! Pad with zeros

C.........  Set left and right portions of SCC
        SCCL = TSCC(       1:LSCCEND )
        SCCR = TSCC( RSCCBEG:SCCLEN3 )

C......... If SCC is non-blank, compare with master SCC list.
        IF( TSCC .NE. SCCZERO ) THEN

            IF( .NOT. CHKEXPSCC( TSCC ) .AND. SCCR .EQ. SCRZERO ) THEN
                IXSCC = FINDC( SCCL, NINVSCL, INVSCL )

            ELSE
                IXSCC = FINDC( TSCC, NINVSCC, INVSCC )

            END IF

C.................  Mobile sources can have zeros for vehicle types and not 
C                   road classes, so check to make sure that this isn't the
C                   case.
            IF( CATEGORY    .EQ. 'MOBILE' .AND.
     &          TSCC( SCCEXPLEN3+1:SCCEXPLEN3+2 ) .EQ. '22'     .AND.
     &          TSCC( SCCEXPLEN3+3:SCCEXPLEN3+6 ) .EQ. '0000'         ) THEN

                SKIPREC = .FALSE.

            ELSE

                SKIPREC = ( IXSCC .LE. 0 )  

            END IF

        END IF

C.........  Check SIC with inventory SIC list.  The record might not match
C           based on SCC, but maybe by SIC.
        IF( ( TSCC == SCCZERO .OR. SKIPREC ) .AND. 
     &        PFLAG .AND. CSIC .NE. SICZERO ) THEN

            IF( .NOT. CHKEXPSIC( CSIC ) .AND. 
     &          CSIC( SICLEN3-1:SICLEN3 ) == '00' ) THEN
                CSIC2 = CSIC( 1:SICLEN3-2 )
                IXSIC = FINDC( CSIC2, NINVSIC2, INVSIC2 )

            ELSE
                IXSIC = FINDC( CSIC, NINVSIC, INVSIC )

            END IF

            SKIPREC = ( IXSIC .LE. 0 )

        END IF

C.........  Check MACT with inventory MACT list
        IF( MFLAG .AND. CMCT /= '000000' ) THEN
            SKIPREC = ( FINDC( CMCT, NINVMACT, INVMACT ) <= 0 )
        END IF 

C.........  Filter the case where the pollutant/activity code is not present
        IF( CPOA .EQ. ' ' ) THEN

            IXPOA = 0

C.........  Ensure that pol/act is in master list of pol/act or
C           skip the pollutant/activity-specific entry
        ELSE

            IXPOA = INDEX1( CPOA, NIPPA, EANAM )

            IF( IXPOA .LE. 0 ) THEN
            
                    SKIPPOA = .TRUE.   ! indicates skipped pol/act-spec entries
                    SKIPREC = .TRUE.   ! indicates skip this rec in calling prgm
   
            END IF

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )

       END SUBROUTINE FLTRXREF
