
        SUBROUTINE FLTRXREF( CFIP, CSIC, TSCC, CPOA, IXSIC, IXSCC,
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
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C...........   This module contains the source arrays
        USE MODSOURC

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
        INTEGER       FIND1
        INTEGER       FINDC
        INTEGER       INDEX1
        INTEGER       STR2INT

        EXTERNAL      FIND1, FINDC, INDEX1, STR2INT

C...........   SUBROUTINE ARGUMENTS:
        CHARACTER(*), INTENT(IN OUT) :: CFIP   ! cntry/state/county code
        CHARACTER(*), INTENT(IN OUT) :: CSIC   ! standard indust. code
        CHARACTER(*), INTENT(IN OUT) :: TSCC   ! source category code
        CHARACTER(*), INTENT(IN OUT) :: CPOA   ! pollutant/activity name
        INTEGER     , INTENT   (OUT) :: IXSIC  ! index of SIC in SIC list
        INTEGER     , INTENT   (OUT) :: IXSCC  ! index of SCC in SCC list
        INTEGER     , INTENT   (OUT) :: IXPOA  ! index of pol/act in master list
        LOGICAL     , INTENT   (OUT) :: SKIPPOA! true: skipped rec is pollutant
                                               ! or activity specific
        LOGICAL     , INTENT   (OUT) :: SKIPREC! true: skip record in caller

C...........   Other local variables

        INTEGER          IXACT   ! index to master activity names array
        INTEGER          L       ! tmp string length
        INTEGER          SIC     ! tmp standard industrial code

        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: 1st time subroutine called 
        LOGICAL, SAVE :: MFLAG    = .FALSE.  ! true: mobile sources
        LOGICAL, SAVE :: PFLAG    = .FALSE.  ! true: point sources

        CHARACTER*300          MESG         ! message buffer
        CHARACTER(LEN=SCCLEN3) SCCL         ! left part of SCC
        CHARACTER(LEN=SCCLEN3) SCCR         ! right part of SCC
        CHARACTER(LEN=SCCLEN3), SAVE :: SCRZERO  ! zero right digits of TSCC
        CHARACTER(LEN=SCCLEN3), SAVE :: SCCZERO  ! zero SCC
 
        CHARACTER*16 :: PROGNAME = 'FLTRXREF' ! program name

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
            MFLAG = ALLOCATED( IRCLAS )
            PFLAG = ALLOCATED( ISIC )

C.............  Check length of SCC string
            L = LEN( TSCC )
            IF( L .NE. SCCLEN3 ) THEN
                MESG = 'INTERNAL ERROR: Length of SCC string not'//
     &                 'correct in call to ' // PROGNAME
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END IF

C.............  Set up zero strings for FIPS code of zero and SCC code of zero
            SCCZERO = REPEAT( '0', SCCLEN3 )
            SCRZERO = REPEAT( '0', SCCLEN3 - LSCCEND )

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
        SIC = STR2INT( CSIC )    ! Convert to integer

C.........  Smart interpretation of SCC
        CALL FLTRNEG( TSCC )     ! Filter 0 and -9 to blank
        CALL PADZERO( TSCC )     ! Pad with zeros

C......... Set left and right portions of SCC
        SCCL = TSCC(       1:LSCCEND )
        SCCR = TSCC( RSCCBEG:SCCLEN3 )

C......... Convert character SCC field to integer SCC number while
C          allowing for case that SCC is blank.  If non-blank, compare
C          with master SCC list for area and point sources.
        IF( TSCC .NE. SCCZERO .AND. .NOT. MFLAG ) THEN

            IF( SCCR .EQ. SCRZERO ) THEN
                IXSCC = FINDC( SCCL, NINVSCL, INVSCL )

            ELSE
                IXSCC = FINDC( TSCC, NINVSCC, INVSCC )

            END IF

            SKIPREC = ( IXSCC .LE. 0 )  

        END IF

C.........  Check SIC with inventory SIC list.  The record might not match
C           based on SCC, but maybe by SIC.
        IF( SKIPREC .AND. PFLAG ) THEN

            IXSIC = FIND1( SIC, NINVSIC, INVSIC )

            SKIPREC = ( IXSIC .LE. 0 )

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
