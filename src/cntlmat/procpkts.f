
        SUBROUTINE PROCPKTS( ADEV, CDEV, GDEV, LDEV, CPYEAR,
     &                       PKTTYP, ENAME, SFLAG )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine is responsible for processing control packet data
C      that has already been grouped into the control cross-reference tables
C      and control data tables.  For control packets that do not depend
C      on pollutants, it calls routines to assign the control to the
C      appropriate sources and generate the appropriate matrix.  For control
C      packets that do depend on pollutants, it creates the pollutant 
C      group structure and processes the packet for all pollutant groups.  The
C      index to the control data tables is stored for all pollutants in
C      temporary files for later use, after the output files have been
C      opened.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Started 3/99 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C.........  This module is for cross reference tables
        USE MODXREF

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2   CRLF
        INTEGER       INDEX1
        INTEGER       STR2INT
        REAL          STR2REAL

        EXTERNAL      CRLF, INDEX1, STR2INT, STR2REAL

C...........   SUBROUTINE ARGUMENTS:

        INTEGER     , INTENT (IN) :: ADEV      ! file unit no. for tmp ADD file
        INTEGER     , INTENT (IN) :: CDEV      ! file unit no. for tmp CTL file 
        INTEGER     , INTENT (IN) :: GDEV      ! file unit no. for tmp CTG file
        INTEGER     , INTENT (IN) :: LDEV      ! file unit no. for tmp ALW file
        INTEGER     , INTENT (IN) :: CPYEAR    ! year to project to 
        CHARACTER(*), INTENT (IN) :: PKTTYP    ! packet type
        CHARACTER(*), INTENT (IN) :: ENAME     ! inventory file name
        LOGICAL     , INTENT(OUT) :: SFLAG     ! true: at least one packet done

C.........  Reshaped inventory pollutants and associated variables
        INTEGER         NGRP                ! number of pollutant groups 
        INTEGER         NGSZ                ! number of pollutants per group 
        INTEGER               , ALLOCATABLE:: IPSTAT ( : )   ! pol status (0|1)
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: EINAM2D( :,: )

C...........   Other local variables

        INTEGER         I, J, K, L, L1, L2, N      ! counters and indices

        INTEGER         IOS                   ! i/o error status
        INTEGER         VIDXMULT( NIPPA )     ! pollutant flags

        LOGICAL       :: EFLAG    = .FALSE.   ! error flag
        LOGICAL, SAVE :: FIRSTIME = .TRUE.    ! true: first time routine called
        LOGICAL       :: OFLAG(NPACKET) = .FALSE.   ! true: tmp file has not
                                                    ! been opened

        CHARACTER*5     CPOS        ! tmp sorted position of pol
        CHARACTER*300   LINE        ! read buffer for a line
        CHARACTER*300   MESG        ! message buffer

        CHARACTER(LEN=IOVLEN3), SAVE :: RPOL ! pol name for reactivity controls

        CHARACTER*16 :: PROGNAME = 'PROCPKTS' ! program name

C***********************************************************************
C   Begin body of subroutine PROCPKTS

        IF( FIRSTIME ) THEN

            MESG = 'Pollutant for creating reactivity matrix'
            CALL ENVSTR( 'REACTIVITY_POL', MESG, 'VOC', RPOL, IOS )

            SFLAG = .FALSE.     ! Initialize status as no packets applied

            FIRSTIME = .FALSE.
   
        END IF

        SELECT CASE( PKTTYP )

C.........  Reactivity packet...
        CASE( 'REACTIVITY' )

C.............  Make sure that the pollutant for the reactivity packet is 
C               in the inventory
            J = INDEX1( RPOL, NIPOL, EINAM )
            IF( J .LE. 0 ) THEN

                MESG = 'Environment variable "REACTIVITY_POL" is set '//
     &                 'to pollutant "' // RPOL( 1:LEN_TRIM( RPOL ) ) //
     &                 '",' // CRLF() // BLANK10 // 'but this ' //
     &                 'pollutant is not in the inventory!'
                CALL M3MSG2( MESG )

                MESG = 'Reactivity matrix creation skipped.'
                CALL M3MSG2( MESG )
                RETURN

            END IF

C.............  Generate reactivity matrices 
            CALL GENREACT( NSRC, NIPOL, BYEAR, CPYEAR, ENAME, 
     &                     RPOL, EINAM )

            SFLAG = .TRUE.

        CASE( 'PROJECTION' )

C.............  Generate projection matrix

            CALL GENPROJ( NSRC, NIPPA, BYEAR, CPYEAR, ENAME, EANAM )

            SFLAG = .TRUE.

        CASE DEFAULT

C.............  It is important that all major arrays must be allocated by this
C               point because the next memory allocation step is going to pick
C               a data structure that will fit within the limits of the host.
C.............  Note that this routine only determines the allocation the
C               first time it is called.
            CALL ALOCCMAT( NGRP, NGSZ )

            IF( .NOT. ALLOCATED( EINAM2D ) ) THEN

C.................  Create 2-dimensional arrays of pollutant names
                ALLOCATE( EINAM2D( NGSZ, NGRP ), STAT=IOS )
                CALL CHECKMEM( IOS, 'EINAM2D', PROGNAME )
                EINAM2D = ' '
                EINAM2D = RESHAPE( EINAM, (/ NGSZ, NGRP /) )

C.................  Create array for indicating the status of pollutants at each
C                   iteration
                ALLOCATE( IPSTAT( NGSZ ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IPSTAT', PROGNAME )

            END IF

C.................  Create array for names of pollutants that receive controls
                ALLOCATE( PNAMMULT( NIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PNAMMULT', PROGNAME )

C.................  Create array of flags indicating which controls are
C                   applied to each pollutant receiving at least one type
C                   of control
                ALLOCATE( PCTLFLAG( NIPPA, 4 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PCTLFALG', PROGNAME )

C.............  Initialize pollutant indictor to zero
            IPSTAT = 0          ! Array

C.............  Initialize control flags to false
            PCTLFLAG = .FALSE.  ! Array
 
C.............  Loop through the pollutant groups...

C.............  Apply the current packets to appropriate sources and pollutants
C               in the inventory.
C.............  Since the table indices are stored separately for each type of
C               control, the assignment subroutine ASGNCNTL must be called using
C               different source-based arrays.
C.............  Write temporary ASCII files containing the indices to the
C               control data tables.  This is because we only want to write out
C               the I/O API control matrices for the pollutants that are 
C               actually affected by controls, but don't know which pollutants
C               to open the output file(s) with until all of the pollutants have 
C               been processed. 

            DO N = 1, NGRP

C.................  Write message stating the pollutants are being processed
                CALL POLMESG( NGSZ, EINAM2D( 1,N ) )

                SELECT CASE( PKTTYP )

                CASE( 'CTG' )

                    CTGIDX = 0   ! array
                    VIDXMULT = 0 ! array
                    CALL ASGNCNTL( NSRC, NGSZ, PKTTYP, EINAM2D( 1,N ), 
     &                             IPSTAT, CTGIDX )

                    CALL UPDATE_POLLIST( N, NGSZ, CTGIDX,
     &                                   VIDXMULT, 2 )
                    IF ( .NOT. OFLAG(1) ) THEN
                       CALL OPENCTMP( PKTTYP, GDEV )
                       OFLAG(1) = .TRUE.
                    END IF
                    CALL WRCTMP( GDEV, N, NGSZ, CTGIDX, VIDXMULT )

                    SFLAG = .TRUE.

                CASE( 'CONTROL' )

                    CTLIDX = 0   ! array
                    VIDXMULT = 0 ! array
                    CALL ASGNCNTL( NSRC, NGSZ, PKTTYP, EINAM2D( 1,N ), 
     &                             IPSTAT, CTLIDX )

                    CALL UPDATE_POLLIST( N, NGSZ, CTLIDX,
     &                                   VIDXMULT, 1 )
                    IF ( .NOT. OFLAG(2) ) THEN
                       CALL OPENCTMP( PKTTYP, CDEV )
                       OFLAG(2) = .TRUE.
                    END IF
                    CALL WRCTMP( CDEV, N, NGSZ, CTLIDX, VIDXMULT )

                    SFLAG = .TRUE.

                CASE( 'ALLOWABLE' )

                    ALWIDX = 0   ! array
                    VIDXMULT = 0 ! array
                    CALL ASGNCNTL( NSRC, NGSZ, PKTTYP, EINAM2D( 1,N ), 
     &                             IPSTAT, ALWIDX )

                    CALL UPDATE_POLLIST( N, NGSZ, ALWIDX,
     &                                   VIDXMULT, 3 )
                    IF ( .NOT. OFLAG(3) ) THEN
                       CALL OPENCTMP( PKTTYP, LDEV )
                       OFLAG(3) = .TRUE.
                    END IF
                    CALL WRCTMP( LDEV, N, NGSZ, ALWIDX, VIDXMULT )

                    SFLAG = .TRUE.

                CASE( 'ADD' )

                    ADDIDX = 0   ! array
                    VIDXMULT = 0 ! array
                    CALL ASGNCNTL( NSRC, NGSZ, PKTTYP, EINAM2D( 1,N ), 
     &                             IPSTAT, ADDIDX )

c                    CALL UPDATE_POLLIST( N, NGSZ, ADDIDX,
c     &                                   VIDXMULT, 4 )
c                    IF ( .NOT. OFLAG(4) ) THEN
c                       CALL OPENCTMP( PKTTYPE, ADEV )
c                       OFLAG(4) = .TRUE.
c                    END IF
c                    CALL WRCTMP( ADEV, N, NGSZ, ADDIDX, VIDXMULT )

                    SFLAG = .TRUE.

                END SELECT

            END DO   ! End loop on pollutant groups

        END SELECT   ! End select on pol-specific packet or not

C...........   Rewind tmp files

        REWIND( ADEV )
        REWIND( CDEV )
        REWIND( GDEV )
        REWIND( LDEV )

        RETURN
       
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000  FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )

C STOPPED HERE:  Need to add UPDATE_POLLIST to maintain the lists of 
C       variable names for the mult control matrix and addative control matrix

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram flags pollutants that have any
C               controls applied.
            SUBROUTINE UPDATE_POLLIST( IGRP, NGSZ, IDX,
     &                                 VIDX, CFLG )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: IGRP            ! pollutant group
            INTEGER     , INTENT (IN) :: NGSZ            ! # pollutants/group
            INTEGER     , INTENT (IN) :: IDX(NSRC,NGSZ)  ! index to data tables
            INTEGER , INTENT (IN OUT) :: VIDX(NIPOL)     ! pollutant flags
            INTEGER     , INTENT (IN) :: CFLG            ! control flag

C.............  Local variables
            INTEGER   I,K,S   ! counters and indices
            LOGICAL   ::      SRCLOOP= .TRUE.   ! true: pollutant 'I' does not
                                                ! have controls applied


C----------------------------------------------------------------------

C.............  Initialize counters and indices
        K = 0

C.............  Loop through all pollutants in the current group
        DO I = 1,NGSZ

           SRCLOOP = .TRUE.
           S = 0
           K = NGSZ*IGRP - NGSZ + I  ! compute index to master pollutant
                                     ! list for current pollutant

C.............  For current pollutant, loop through sources until a source with
C               controls is encountered. Terminate loop when all sources have
C               been examined.            
           DO WHILE( SRCLOOP .AND. S .LT. NSRC )

              S = S + 1
              IF ( IDX(S,I) .GT. 0 ) SRCLOOP = .FALSE.  ! controls encountered,
                                                        ! exit source loop
           END DO  ! end source loop

C.............  Check to see if current pollutant has controls applied
           IF ( SRCLOOP ) THEN       ! no controls
              VIDX(K) = 0
           ELSE                      ! controls
              NVCMULT  = NVCMULT + 1
              VIDX(K)  = 1
              PCTLFLAG( NVCMULT, CFLG ) = .TRUE.
              PNAMMULT( NVCMULT ) = EINAM( K )
           END IF

        END DO  ! end pollutant group loop

        RETURN

        END SUBROUTINE UPDATE_POLLIST

C----------------------------------------------------------------------

       END SUBROUTINE PROCPKTS
