
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  Local allocatable parameters
        INTEGER, ALLOCATABLE :: AFIPS ( : )  ! area source co/st/cy codes
        INTEGER, ALLOCATABLE :: MFIPS ( : )  ! mobile source co/st/cy codes
        INTEGER, ALLOCATABLE :: PFIPS ( : )  ! point source co/st/cy codes
        INTEGER, ALLOCATABLE :: TMPFIP( : )  ! all unsort source co/st/cy codes
        INTEGER, ALLOCATABLE :: SIDX  ( : )  ! sorting index

C.........  Local parameters
        INTEGER       J, K, L, N, S      ! pointers and counters

        INTEGER       IOS          ! i/o status
        INTEGER       MXNFIP       ! max no. FIPS codes
        INTEGER    :: NAFIP = 0
        INTEGER    :: NMFIP = 0
        INTEGER    :: NPFIP = 0
        INTEGER       PFIP         ! tmp fip from previous iteration

        CHARACTER*33 ,PARAMETER :: PART1 = 
     &                             'Error reading variable IFIP from '
        CHARACTER*15 ,PARAMETER :: PART3 = ' INVENTORY file'

C...........   Other local variables

        CHARACTER*300    MESG     ! message buffer

        CHARACTER*16  :: PROGNAME = 'RDMRGINV' ! program name

C***********************************************************************
C   begin body of subroutine RDMRGINV

        IF( ALLOCATED( INVIFIP ) ) DEALLOCATE( INVIFIP )
        IF( ALLOCATED( AIFIP   ) ) DEALLOCATE( AIFIP )
        IF( ALLOCATED( AFIPS   ) ) DEALLOCATE( AFIPS  )
        IF( ALLOCATED( MFIPS   ) ) DEALLOCATE( MFIPS  )
        IF( ALLOCATED( PFIPS   ) ) DEALLOCATE( PFIPS  )
        IF( ALLOCATED( TMPFIP  ) ) DEALLOCATE( TMPFIP )
        IF( ALLOCATED( SIDX    ) ) DEALLOCATE( SIDX   )

C.........  Read in area source country, state, county code
        IF( AFLAG ) THEN

            ALLOCATE( AIFIP( NASRC ), STAT=IOS )   ! country/state/county codes
            CALL CHECKMEM( IOS, 'AIFIP', PROGNAME )

            IF( .NOT. READ3( AENAME,'IFIP',ALLAYS3,0,0,AIFIP ) ) THEN
                MESG = PART1 // 'AREA' // PART3
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

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

            IF( .NOT. READ3( MENAME,'IFIP',ALLAYS3,0,0,MIFIP ) ) THEN
                MESG = PART1 // 'MOBILE' // PART3
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

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

            IF( .NOT. READ3( PENAME,'IFIP',ALLAYS3,0,0,PIFIP ) ) THEN
                MESG = PART1 // 'POINT' // PART3
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

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
        ALLOCATE( INVIFIP( MXNFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVIFIP', PROGNAME )
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
        CALL SORTI1( MXNFIP, SIDX, TMPFIP )

C.........  Store unique list across all source categories
        N = 0
        PFIP =  -9
        DO K = 1, MXNFIP

            J = SIDX( K )
            IF( TMPFIP( J ) .NE. PFIP ) THEN
                N = N + 1
                INVIFIP( N ) = TMPFIP( J )
            END IF

            PFIP = TMPFIP( J )

        END DO

        NINVIFIP = N

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C*****************  INTERNAL SUBPROGRAMS  ******************************

        CONTAINS

C.............  This subprogram counts the FIPS codes in the sorted 
C               input array of codes
            SUBROUTINE COUNT_FIPS( NSRC, IFIP, NFIP )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: NSRC
            INTEGER, INTENT (IN) :: IFIP( NSRC )
            INTEGER, INTENT(OUT) :: NFIP

C.............  Local variables
            INTEGER  K
            INTEGER  PFIP   ! from previous iteration

C----------------------------------------------------------------------

            K = 0
            PFIP = -9
            DO S = 1, NSRC

               IF( IFIP( S ) .NE. PFIP ) K = K + 1
               PFIP = IFIP( S )

            END DO
            NFIP = K

            END SUBROUTINE COUNT_FIPS

C....................................................................
C....................................................................

C.............  This subprogram creates a list of unique FIPS codes
C               from the sorted input list
            SUBROUTINE CREATE_FIPS( NSRC, NFIP, IFIP, FIPS )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: NSRC
            INTEGER, INTENT (IN) :: NFIP
            INTEGER, INTENT (IN) :: IFIP( NSRC )
            INTEGER, INTENT(OUT) :: FIPS( NFIP )

C.............  Local variables
            INTEGER  K
            INTEGER  PFIP   ! from previous iteration

C----------------------------------------------------------------------

            K = 0
            PFIP = -9
            DO S = 1, NSRC

               IF( IFIP( S ) .NE. PFIP ) THEN
                   K = K + 1
                   FIPS( K ) = IFIP( S )
               END IF

               PFIP = IFIP( S )

            END DO

            END SUBROUTINE CREATE_FIPS

        END SUBROUTINE RDMRGINV
