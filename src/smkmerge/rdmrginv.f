
        SUBROUTINE RDMRGINV

C***********************************************************************
C  subroutine RDMRGINV body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to read in inventory variables needed
C      for the merge program
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE M3UTILIO

        USE MODMERGE, ONLY: AFLAG,  MFLAG,  PFLAG,      ! source flags by category
     &                      AENAME, MENAME, PENAME,     ! inventory file names
     &                      ASDEV,  MSDEV,  PSDEV,      ! inventory file names
     &                      AIFIP,  MIFIP,  PIFIP,      ! country/state/county codes
     &                      NASRC,  NMSRC,  NPSRC       ! no. of sources

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP

C...........   This module contains the source arrays
        USE MODSOURC, ONLY: CIFIP

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.........  Local allocatable parameters
        CHARACTER(FIPLEN3), ALLOCATABLE :: AFIPS ( : )  ! area source co/st/cy codes
        CHARACTER(FIPLEN3), ALLOCATABLE :: MFIPS ( : )  ! mobile source co/st/cy codes
        CHARACTER(FIPLEN3), ALLOCATABLE :: PFIPS ( : )  ! point source co/st/cy codes
        CHARACTER(FIPLEN3), ALLOCATABLE :: TMPFIP( : )  ! all unsort source co/st/cy codes
        INTEGER, ALLOCATABLE :: SIDX  ( : )  ! sorting index

C.........  Local parameters
        INTEGER       J, K, L, N, S      ! pointers and counters

        INTEGER       IOS          ! i/o status
        INTEGER       MXNFIP       ! max no. FIPS codes

        INTEGER    :: NAFIP = 0
        INTEGER    :: NMFIP = 0
        INTEGER    :: NPFIP = 0
        
        CHARACTER(FIPLEN3)    PFIP         ! tmp fip from previous iteration

        CHARACTER(33), PARAMETER :: PART1 = 
     &                             'Error reading variable CIFIP from '
        CHARACTER(15), PARAMETER :: PART3 = ' INVENTORY file'

C...........   Other local variables

        CHARACTER(300)   MESG     ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDMRGINV' ! program name

C***********************************************************************
C   begin body of subroutine RDMRGINV

C.........  Read in area source country, state, county code
        IF( AFLAG ) THEN

            ALLOCATE( AIFIP( NASRC ), STAT=IOS )   ! country/state/county codes
            CALL CHECKMEM( IOS, 'AIFIP', PROGNAME )

            CALL RDINVCHR( 'AREA', AENAME, ASDEV, NASRC, 1, 'CIFIP' )
            AIFIP = CIFIP

C.............  Count area source country, state, and county codes
            CALL COUNT_FIPS( NASRC, AIFIP, NAFIP )

C.............  Allocate memory for codes
            ALLOCATE( AFIPS( NAFIP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AFIPS', PROGNAME )

C.............  Create temporary area-source FIPs code list
            CALL CREATE_FIPS( NASRC, NAFIP, AIFIP, AFIPS )

        END IF


C.........  Read in mobile source country, state, county code
        IF( MFLAG ) THEN
            ALLOCATE( MIFIP( NMSRC ), STAT=IOS )   ! country/state/county codes
            CALL CHECKMEM( IOS, 'MIFIP', PROGNAME )

            CALL RDINVCHR( 'MOBILE', MENAME, MSDEV, NMSRC, 1, 'CIFIP' )
            MIFIP = CIFIP

C.............  Count mobile source country, state, and county codes
            CALL COUNT_FIPS( NMSRC, MIFIP, NMFIP )

C.............  Allocate memory for codes
            ALLOCATE( MFIPS( NMFIP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MFIPS', PROGNAME )

C.............  Create temporary mobile-source FIPs code list
            CALL CREATE_FIPS( NMSRC, NMFIP, MIFIP, MFIPS )
        END IF

C.........  Create temporary mobile-source FIPs code list

C.........  Read in point source country, state, county code
        IF( PFLAG ) THEN
            ALLOCATE( PIFIP( NPSRC ), STAT=IOS )   ! country/state/county codes
            CALL CHECKMEM( IOS, 'PIFIP', PROGNAME )

            CALL RDINVCHR( 'POINT', PENAME, PSDEV, NPSRC, 1, 'CIFIP' )
            PIFIP = CIFIP

C.............  Count point source country, state, and county codes
            CALL COUNT_FIPS( NPSRC, PIFIP, NPFIP )

C.............  Allocate memory for codes
            ALLOCATE( PFIPS( NPFIP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PFIPS', PROGNAME )

C.............  Create temporary point-source FIPs code list
            CALL CREATE_FIPS( NPSRC, NPFIP, PIFIP, PFIPS )

        END IF

C.........  Combine temporary lists of FIPs codes into single sorted list for 
C           reading state and counties file
        MXNFIP = NAFIP + NMFIP + NPFIP
        ALLOCATE( TMPFIP( MXNFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPFIP', PROGNAME )
        ALLOCATE( INVCFIP( MXNFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVCFIP', PROGNAME )
        ALLOCATE( SIDX( MXNFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SIDX', PROGNAME )

C.........  Add area sources to complete FIPS list
        N = 0
        DO K = 1, NAFIP
            N = N + 1
            SIDX  ( N ) = N
            TMPFIP( N ) = AFIPS( K )
        END DO

C.........  Add mobile sources to complete FIPS list
        DO K = 1, NMFIP
            N = N + 1
            SIDX  ( N ) = N
            TMPFIP( N ) = MFIPS( K )
        END DO

C.........  Add point sources to complete FIPS list
        DO K = 1, NPFIP
            N = N + 1
            SIDX  ( N ) = N
            TMPFIP( N ) = PFIPS( K )
        END DO

        MXNFIP = N

C.........  Sort all FIPS codes
        CALL SORTIC( MXNFIP, SIDX, TMPFIP )

C.........  Store unique list across all source categories
        N = 0
        PFIP = ' '
        DO K = 1, MXNFIP

            J = SIDX( K )
            IF( TMPFIP( J ) .NE. PFIP ) THEN
                N = N + 1
                INVCFIP( N ) = TMPFIP( J )
            END IF

            PFIP = TMPFIP( J )

        END DO

        NINVIFIP = N

C.........  Deallocate local memory
    	IF( ALLOCATED( AFIPS ) ) DEALLOCATE( AFIPS )
    	IF( ALLOCATED( MFIPS ) ) DEALLOCATE( MFIPS )
    	IF( ALLOCATED( PFIPS ) ) DEALLOCATE( PFIPS )
    	IF( ALLOCATED( TMPFIP ) ) DEALLOCATE( TMPFIP )
    	IF( ALLOCATED( SIDX ) ) DEALLOCATE( SIDX )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C*****************  INTERNAL SUBPROGRAMS  ******************************

        CONTAINS

C.............  This subprogram counts the FIPS codes in the sorted 
C               input array of codes
            SUBROUTINE COUNT_FIPS( NSRC, CIFIP, NFIP )

C.............  Subprogram arguments
            INTEGER,            INTENT (IN) :: NSRC
            CHARACTER(FIPLEN3), INTENT (IN) :: CIFIP( NSRC )
            INTEGER,            INTENT(OUT) :: NFIP

C.............  Local variables
            INTEGER  K

            CHARACTER(FIPLEN3) PFIP   ! from previous iteration

C----------------------------------------------------------------------

            K = 0
            PFIP = ' '
            DO S = 1, NSRC

               IF( CIFIP( S ) .NE. PFIP ) K = K + 1
               PFIP = CIFIP( S )

            END DO
            NFIP = K

            END SUBROUTINE COUNT_FIPS

C....................................................................
C....................................................................

C.............  This subprogram creates a list of unique FIPS codes
C               from the sorted input list
            SUBROUTINE CREATE_FIPS( NSRC, NFIP, CIFIP, CFIPS )

C.............  Subprogram arguments
            INTEGER,            INTENT (IN) :: NSRC
            INTEGER,            INTENT (IN) :: NFIP
            CHARACTER(FIPLEN3), INTENT (IN) :: CIFIP( NSRC )
            CHARACTER(FIPLEN3), INTENT(OUT) :: CFIPS( NFIP )

C.............  Local variables
            INTEGER  K

            CHARACTER(FIPLEN3)  PFIP   ! from previous iteration

C----------------------------------------------------------------------

            K = 0
            PFIP = ' '
            DO S = 1, NSRC

               IF( CIFIP( S ) .NE. PFIP ) THEN
                   K = K + 1
                   CFIPS( K ) = CIFIP( S )
               END IF

               PFIP = CIFIP( S )

            END DO

            END SUBROUTINE CREATE_FIPS

        END SUBROUTINE RDMRGINV
