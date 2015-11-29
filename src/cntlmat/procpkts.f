
        SUBROUTINE PROCPKTS( PDEV, CDEV, GDEV, LDEV, MDEV, WDEV, CPYEAR,
     &                       PKTTYP, ENAME, LPSASGN, USEPOL, SFLAG,
     &                       LPTMP, LCTMP )

C***********************************************************************
C  subroutine body starts at line 116
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

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: CSOURC

C.........  This module is for cross reference tables
        USE MODXREF, ONLY: ASGNINDX

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: PNAMMULT, PNAMPROJ, FACTOR, 
     &                      BACKOUT, DATVAL, GRPINDX, GRPFLAG, GRPSTIDX,
     &                      GRPCHAR, GRPINEM, GRPOUTEM, POLSFLAG, 
     &                      NVPROJ, NVCMULT, PCTLFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NSRC, NIPPA, NIPOL, NIACT, NPPOL,
     &                     EINAM, EANAM, ACTVTY

        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        CHARACTER(2)  CRLF
        INTEGER       INDEX1
        INTEGER       STR2INT
        REAL          STR2REAL

        EXTERNAL      CRLF, INDEX1, STR2INT, STR2REAL

C...........   SUBROUTINE ARGUMENTS:

        INTEGER     , INTENT (IN) :: PDEV      ! file unit no. for tmp PROJ file
        INTEGER     , INTENT (IN) :: CDEV      ! file unit no. for tmp CTL file 
        INTEGER     , INTENT (IN) :: GDEV      ! file unit no. for tmp CTG file
        INTEGER     , INTENT (IN) :: LDEV      ! file unit no. for tmp ALW file
        INTEGER     , INTENT (IN) :: MDEV      ! file unit no. for tmp MACT file
        INTEGER     , INTENT (IN) :: WDEV      ! file unit no. for warnings/error file
        INTEGER     , INTENT (IN) :: CPYEAR    ! year to project to 
        CHARACTER(*), INTENT (IN) :: PKTTYP    ! packet type
        CHARACTER(*), INTENT (IN) :: ENAME     ! inventory file name
        LOGICAL     , INTENT (IN) :: LPSASGN   ! true: matrix needs to be by pollutant
        LOGICAL , INTENT (IN OUT) :: USEPOL( NIPPA ) ! true: pol in current pkt
        LOGICAL     , INTENT(OUT) :: SFLAG     ! true: at least one packet done
        LOGICAL     , INTENT(OUT) :: LPTMP     ! true: projection tmp file written
        LOGICAL     , INTENT(OUT) :: LCTMP     ! true: control tmp file written

C.........  Reshaped inventory pollutants and associated variables
c        INTEGER         NGRP                ! number of pollutant groups 
c        INTEGER         NGSZ                ! number of pollutants per group 
c        INTEGER               , ALLOCATABLE:: IPSTAT ( : )   ! pol status (0|1)

C...........   Other local variables

        INTEGER         I, J, K, L, L1, L2, S, V      ! counters and indices

        INTEGER         IOS                   ! i/o error status
        INTEGER         NGRP                  ! number of reporting groups
        INTEGER         SCCBEG                ! begining of SCC in CSOURC string
        INTEGER         SCCEND                ! end of SCC in CSOURC string
        INTEGER         VIDXMULT( NIPPA )     ! pollutant flags for control
        INTEGER         VIDXPROJ( NIPPA )     ! pollutant flags for projection

        LOGICAL       :: DSFLAG   = .FALSE.   ! true: this call to ASGNCNTL has data-specific match
        LOGICAL       :: EFLAG    = .FALSE.   ! error flag
        LOGICAL, SAVE :: FIRSTIME = .TRUE.    ! true: first time routine called
        LOGICAL, SAVE :: OFLAG(NPACKET) = .FALSE.   ! true: tmp file has not been opened
        LOGICAL       :: CCHEK                ! temporary control tmp file status per packet

        CHARACTER(5)    CPOS        ! tmp sorted position of pol
        CHARACTER(256)  LINE        ! read buffer for a line
        CHARACTER(256)  MESG        ! message buffer

        CHARACTER(IOVLEN3), SAVE :: RPOL ! pol name for reactivity controls
        CHARACTER(FPLLEN3) CPLT       ! tmp point src info through plant
        CHARACTER(FPLLEN3) PPLT       ! previous CPLT
        CHARACTER(STALEN3) CSTA       ! tmp char state
        CHARACTER(STALEN3) PSTA       ! previous char state
        CHARACTER(SCCLEN3) TSCC       ! tmp SCC
        CHARACTER(SCCLEN3) PSCC       ! previous SCC

        CHARACTER(16) :: PROGNAME = 'PROCPKTS' ! program name

C***********************************************************************
C   Begin body of subroutine PROCPKTS

        IF( FIRSTIME ) THEN

            SFLAG = .FALSE.     ! Initialize status as no packets applied
            FIRSTIME = .FALSE.
   
        END IF

C.........  Reactivity packet...
        IF( PKTTYP .EQ. 'REACTIVITY' ) THEN

C.............  Get environment variable setting for reactivity pollutant
            MESG = 'Pollutant for creating reactivity matrix'
            CALL ENVSTR( 'REACTIVITY_POL', MESG, 'VOC', RPOL, IOS )

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
            USEPOL = .TRUE.  ! array
            CALL GENREACT( CPYEAR, ENAME, RPOL, USEPOL )

            SFLAG = .TRUE.

C.........  If projection or control packet, then allocate memory for 
C           needed arrays, pollutant-specific arrays, and output matrices.
        ELSE

C............  Allocate arrays for managing output pollutants/activities
            IF( .NOT. ALLOCATED( PNAMMULT ) ) THEN

C................  Create array for names of pollutants that receive controls
                ALLOCATE( PNAMMULT( NIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PNAMMULT', PROGNAME )
                PNAMMULT = ' '    ! array

C................  Create array for names of pollutants that receive projections
                ALLOCATE( PNAMPROJ( NIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PNAMPROJ', PROGNAME )
                PNAMPROJ = ' '    ! array

C.................  Create array of flags indicating which controls are
C                   applied to each pollutant receiving at least one type
C                   of control or projection
                ALLOCATE( PCTLFLAG( NIPPA, 4 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PCTLFALG', PROGNAME )
                PCTLFLAG = .FALSE.    ! array

            END IF

C...........  Allocate for work output factor for both projection and control
            IF( .NOT. ALLOCATED( FACTOR ) ) THEN
                ALLOCATE( ASGNINDX( NSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ASGNINDX', PROGNAME )
                ALLOCATE( FACTOR( NSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'FACTOR', PROGNAME )
            END IF
            ASGNINDX = 0
            FACTOR = 1.  !  array

C...........  Allocate for work arrays for multiplicative controls
            IF( PKTTYP .NE. 'PROJECTION' ) THEN

C...............  Allocate main work arrays
                IF( .NOT. ALLOCATED( BACKOUT ) ) THEN
                    ALLOCATE( BACKOUT( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'BACKOUT', PROGNAME )
                    ALLOCATE( DATVAL( NSRC,NPPOL ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'DATVAL', PROGNAME )
                END IF
                BACKOUT = 0.   ! array
                DATVAL  = 0.   ! array

C...............  Allocate first set of reporting arrays
                IF( .NOT. ALLOCATED( GRPINDX ) ) THEN
                    ALLOCATE( GRPINDX( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'GRPINDX', PROGNAME )
                END IF
                IF( .NOT. ALLOCATED( GRPSTIDX ) ) THEN
                    ALLOCATE( GRPSTIDX( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'GRPSTIDX', PROGNAME )
                    ALLOCATE( GRPCHAR( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'GRPCHAR', PROGNAME )
                    GRPINDX  = 0  ! array
                END IF

C...............  If haven't already, get set up for group reporting...
                IF ( .NOT. ALLOCATED( GRPFLAG ) ) THEN

C..................  Count the number of groups in the inventory
                    IF( CATEGORY .EQ. 'POINT' ) THEN

                        PPLT = ' '
                        NGRP = 0
                        DO S = 1, NSRC
                            CPLT = CSOURC( S )( 1:FPLLEN3 )
                            IF( CPLT .NE. PPLT ) THEN
                                NGRP = NGRP + 1
                                PPLT = CPLT
                            END IF
                            GRPINDX ( S ) = NGRP
                            GRPSTIDX( S ) = S     ! Needed for loops, but not used to sort
                        END DO

                    ELSE 

                        IF( CATEGORY .EQ. 'AREA' ) THEN
                            SCCBEG = ARBEGL3( 2 )
                            SCCEND = ARENDL3( 2 )
                        ELSE            ! MOBILE
                            SCCBEG = MBBEGL3( 5 )
                            SCCEND = MBENDL3( 5 )
                        END IF

C......................  Build and sort source array for SCC-state grouping
                        DO S = 1, NSRC
                            CSTA = CSOURC( S )( 1     :STALEN3 )
                            TSCC = CSOURC( S )( SCCBEG:SCCEND  )

                            GRPSTIDX( S ) = S  
                            GRPCHAR ( S ) = CSTA // TSCC
                        END DO

                        CALL SORTIC( NSRC, GRPSTIDX, GRPCHAR )

C......................  Count the number of state/SCCs in the domain
                        PSTA = ' '
                        PSCC = ' '
                        SCCBEG = STALEN3 + 1
                        SCCEND = STALEN3 + SCCLEN3
                        DO S = 1, NSRC
                            J = GRPSTIDX( S )
                            CSTA = GRPCHAR( J )( 1     :STALEN3 )
                            TSCC = GRPCHAR( J )( SCCBEG:SCCEND  )
                            IF( CSTA .NE. PSTA .OR. 
     &                          TSCC .NE. PSCC       ) THEN
                                NGRP = NGRP + 1
                                PSTA = CSTA
                                PSCC = TSCC
                            END IF
                            GRPINDX( J ) = NGRP
                        END DO

                    END IF    ! If point or non-point sources 

                    IF( .NOT. ALLOCATED( GRPFLAG ) ) THEN
                        ALLOCATE( GRPFLAG( NGRP ), STAT=IOS )
                        CALL CHECKMEM( IOS, 'GRPFLAG', PROGNAME )
                        ALLOCATE( GRPINEM( NGRP, NIPPA ), STAT=IOS )
                        CALL CHECKMEM( IOS, 'GRPINEM', PROGNAME )
                        ALLOCATE( GRPOUTEM( NGRP, NIPPA ), STAT=IOS )
                        CALL CHECKMEM( IOS, 'GRPOUTEM', PROGNAME )
                    END IF
                END IF    ! If group information not previously allocated

                GRPINEM  = 0. ! array
                GRPOUTEM = 0. ! array
                GRPFLAG  = .FALSE.  ! array

            END IF   ! For multiplicative control packets

C.................  Create array for indicating the status of pollutants at each
C                   iteration
c                ALLOCATE( IPSTAT( NGSZ ), STAT=IOS )
c                CALL CHECKMEM( IOS, 'IPSTAT', PROGNAME )

c            IPSTAT = 0          ! Array 

        END IF       ! For reactivity or not

C.........  If the packet does not have pol/act-specific entries, then
C              create using the "all" keyword.
C.........  For projection matrix only for now...
        IF ( PKTTYP .EQ. 'PROJECTION' .AND. .NOT. LPSASGN ) THEN

            CALL ASGNCNTL( NSRC, WDEV, PKTTYP, 'all', DSFLAG, 
     &                     ASGNINDX )
            
            IF ( .NOT. OFLAG(1) ) THEN
               CALL OPENCTMP( PKTTYP, PDEV )
               OFLAG(1) = .TRUE.
            END IF
            CALL WRCTMP( PDEV, 1, ASGNINDX, 1, LPTMP )

            SFLAG = .TRUE. 

C.........  For projection for specific pollutants and multiplicative controls...
        ELSE IF ( PKTTYP .NE. 'REACTIVITY' ) THEN

C.............  Loop through the pollutant groups...

C.............  Apply the current packets to appropriate sources and pollutants
C               in the inventory.
C.............  Write temporary ASCII files containing the indices to the
C               control data tables.  This is because we only want to write out
C               the I/O API control matrices for the pollutants that are 
C               actually affected by controls, but don't know which pollutants
C               to open the output file(s) with until all of the pollutants have 
C               been processed. 

            DO V = 1, NIPPA

C...............  Skip pollutants that do not apply to this packet.
                IF( .NOT. USEPOL( V ) ) CYCLE

                SELECT CASE( PKTTYP )

                CASE( 'PROJECTION' )

C......................  Generate projection matrix
                    CALL ASGNCNTL( NSRC, WDEV, PKTTYP, EANAM( V ),
     &                             DSFLAG, ASGNINDX )
                    IF( DSFLAG ) POLSFLAG = .TRUE.

                    VIDXPROJ = 0 ! array
                    CALL UPDATE_POLLIST( V, ASGNINDX, 0, VIDXPROJ,
     &                                   NVPROJ, PNAMPROJ )

                    IF ( .NOT. OFLAG(1) ) THEN
                       CALL OPENCTMP( PKTTYP, PDEV )
                       OFLAG(1) = .TRUE.
                    END IF
                    CALL WRCTMP( PDEV, V, ASGNINDX, VIDXPROJ, LPTMP )

                    SFLAG = .TRUE. 

                CASE( 'CTG' )

                    VIDXMULT = 0 ! array
                    CALL ASGNCNTL( NSRC, WDEV, PKTTYP, EANAM( V ), 
     &                             DSFLAG, ASGNINDX )

                    CALL UPDATE_POLLIST( V, ASGNINDX, 2, VIDXMULT,
     &                                   NVCMULT, PNAMMULT )
                    IF ( .NOT. OFLAG(2) ) THEN
                       CALL OPENCTMP( PKTTYP, GDEV )
                       OFLAG(2) = .TRUE.
                    END IF
                    
                    CALL WRCTMP( GDEV, V, ASGNINDX, VIDXMULT, CCHEK )
                    LCTMP = ( LCTMP .OR. CCHEK )

                    SFLAG = .TRUE.

                CASE( 'CONTROL' )

C...................  Skip activities because they
C                     do not have the base-year control effectiveness
                    J = INDEX1( EANAM( V ), NIACT, ACTVTY )
                    IF ( J .GT. 0 ) THEN
                        MESG = 'Skipping activity "' //
     &                             TRIM( EANAM( V ) )// '" since '//
     &                             'CONTROL packet cannot apply.'
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF

                    VIDXMULT = 0 ! array
                    CALL ASGNCNTL( NSRC, WDEV, PKTTYP, EANAM( V ),
     &                             DSFLAG, ASGNINDX )

                    CALL UPDATE_POLLIST( V, ASGNINDX, 1, VIDXMULT,
     &                                   NVCMULT, PNAMMULT )
                    IF ( .NOT. OFLAG(3) ) THEN
                       CALL OPENCTMP( PKTTYP, CDEV )
                       OFLAG(3) = .TRUE.
                    END IF
                    CALL WRCTMP( CDEV, V, ASGNINDX, VIDXMULT, CCHEK )
                    LCTMP = ( LCTMP .OR. CCHEK )

                    SFLAG = .TRUE.

                CASE( 'ALLOWABLE' )

                    VIDXMULT = 0 ! array
                    CALL ASGNCNTL( NSRC, WDEV, PKTTYP, EANAM( V ), 
     &                             DSFLAG, ASGNINDX )

                    CALL UPDATE_POLLIST( V, ASGNINDX, 3, VIDXMULT,
     &                                   NVCMULT, PNAMMULT )
                    IF ( .NOT. OFLAG(4) ) THEN
                       CALL OPENCTMP( PKTTYP, LDEV )
                       OFLAG(4) = .TRUE.
                    END IF
                    CALL WRCTMP( LDEV, V, ASGNINDX, VIDXMULT, CCHEK )
                    LCTMP = ( LCTMP .OR. CCHEK )

                    SFLAG = .TRUE.

                CASE( 'MACT' )
                
                    VIDXMULT = 0 ! array
                    CALL ASGNCNTL( NSRC, WDEV, PKTTYP, EANAM( V ),
     &                             DSFLAG, ASGNINDX )
                                   
                    CALL UPDATE_POLLIST( V, ASGNINDX, 4, VIDXMULT,
     &                                   NVCMULT, PNAMMULT )
     
                    IF ( .NOT. OFLAG(5) ) THEN
                        CALL OPENCTMP( PKTTYP, MDEV )
                        OFLAG(5) = .TRUE.
                    END IF
                    CALL WRCTMP( MDEV, V, ASGNINDX, VIDXMULT, CCHEK )
                    LCTMP = ( LCTMP .OR. CCHEK )
                    
                    SFLAG = .TRUE.

                END SELECT

            END DO   ! End loop on pollutants and activities

        END IF       ! End select on pol-specific packet or not

C...........   Rewind tmp files
        IF( PDEV .GT. 0 ) REWIND( PDEV )
        IF( CDEV .GT. 0 ) REWIND( CDEV )
        IF( GDEV .GT. 0 ) REWIND( GDEV )
        IF( LDEV .GT. 0 ) REWIND( LDEV )
        IF( MDEV .GT. 0 ) REWIND( MDEV )

        RETURN
       
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000  FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram flags pollutants that have any
C               controls applied.
            SUBROUTINE UPDATE_POLLIST( POLID, IDX, CFLG, VIDX,
     &                                 NPCNT, PNAMES )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: POLID            ! pollutant index
            INTEGER     , INTENT (IN) :: IDX(NSRC)  ! index to data tables
            INTEGER     , INTENT (IN) :: CFLG            ! control flag
            INTEGER , INTENT (IN OUT) :: VIDX(NIPPA)     ! pollutant flags
            INTEGER , INTENT (IN OUT) :: NPCNT           ! pollutant/act count
            CHARACTER(*), INTENT (IN OUT) :: PNAMES(NIPPA)   ! pollutant/act names

C.............  Local variables
            INTEGER   I, S   ! counters and indices
            LOGICAL   ::      SRCLOOP= .TRUE.   ! true: pollutant 'I' does not
                                                ! have controls applied

C----------------------------------------------------------------------

C.............  For current pollutant, loop through sources until a source with
C               controls is encountered. Terminate loop when all sources have
C               been examined.
        SRCLOOP = .TRUE.
        S = 0
        DO WHILE( SRCLOOP .AND. S .LT. NSRC )

           S = S + 1
           IF ( IDX(S) .GT. 0 ) SRCLOOP = .FALSE.  ! controls encountered,
                                                     ! exit source loop
        END DO  ! end source loop

C.............  Check to see if current pollutant has controls applied
        IF ( SRCLOOP ) THEN       ! no controls
           VIDX( POLID ) = 0

        ELSE                      ! controls

C............  See if pollutant is already in list
           I = INDEX1( EANAM( POLID ), NPCNT, PNAMES )

C............ If not, add it
           IF( I .LE. 0 ) THEN
              NPCNT  = NPCNT + 1
              PNAMES( NPCNT ) = EANAM( POLID )
              I = NPCNT
           END IF

C..............  Update pol/act flag for this packet
           VIDX( POLID )  = 1

C..............  Update control type flag if control type is set
           IF( CFLG .GT. 0 ) PCTLFLAG( I, CFLG ) = .TRUE.

        END IF

        RETURN

        END SUBROUTINE UPDATE_POLLIST

C----------------------------------------------------------------------

       END SUBROUTINE PROCPKTS
