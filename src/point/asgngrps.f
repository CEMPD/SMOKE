
        SUBROUTINE ASGNGRPS( NSP, NCRIT, MXCHK, CRITVALS, 
     &                       CRITYPES, NINVGRP )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION: 
C     This routine assigns groups to the inventory sources based on comparison
C     criteria provided as subroutine arguments.
C
C  PRECONDITIONS REQUIRED:
C     Comparison criteria must be in a specific format in which the "ORs" are
C     on separate rows and the "ANDs" are on the same row.  Modules MODINFO
C     and MODSOURC must be populated, and MODELEV unallocated.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     EVALCRIT is responsible for testing the formulas
C
C  REVISION  HISTORY:
C     Written 7/2001 by M. Houyoux
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

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   ! emissions constant parameters
        INCLUDE 'PARMS3.EXT'    ! I/O API constants
        INCLUDE 'CONST3.EXT'    ! physical and mathematical constants
        INCLUDE 'FLTERR.EXT'    ! error filter statement function

C...........   ARGUMENTS and their descriptions:
        INTEGER     , INTENT (IN) :: NSP     ! no. of variables in formulas
        INTEGER     , INTENT (IN) :: NCRIT   ! no. of comparison criteria
        INTEGER     , INTENT (IN) :: MXCHK   ! max. no. checks per NSP and NCRIT
        REAL        , INTENT (IN) :: CRITVALS( NCRIT, MXCHK, NSP ) ! formula values
        CHARACTER(*), INTENT (IN) :: CRITYPES( NCRIT, MXCHK, NSP ) ! formula types
        INTEGER     , INTENT(OUT) :: NINVGRP ! no. of groups

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2 CRLF
        LOGICAL     EVALCRIT

        EXTERNAL    CRLF, EVALCRIT

C...........   LOCAL PARAMETERS and their descriptions:
        INTEGER, PARAMETER :: MXLOCGRP = 2000  ! Max number groups per facility
        INTEGER, PARAMETER :: MXSPGRP  = 500   ! Max number sources per group

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: TGRPCNT( : ) ! tmp no. sources in group
        REAL   , ALLOCATABLE :: TGRPHT ( : ) ! tmp stack heights [m]
        REAL   , ALLOCATABLE :: TGRPDM ( : ) ! tmp stack diameters [m]
        REAL   , ALLOCATABLE :: TGRPTK ( : ) ! tmp stack exit temperatures [K]
        REAL   , ALLOCATABLE :: TGRPVE ( : ) ! tmp stack exit velocities [m/s]
        REAL   , ALLOCATABLE :: TGRPFL ( : ) ! tmp stack flow rates [m^3/s]

        REAL   , ALLOCATABLE :: VALS   ( : ) ! tmp source stack parameter array
        REAL   , ALLOCATABLE :: REFS   ( : ) ! tmp group stack parameter array

        LOGICAL, ALLOCATABLE :: GPSTAT( :,:,: ) ! not used except EVALCRIT call
        CHARACTER(LEN=PLTLEN3), ALLOCATABLE :: CHRS( : )        ! dummy test character strings
        CHARACTER(LEN=PLTLEN3), ALLOCATABLE :: COMCHRS( :,:,: ) ! dummy formula strings

C...........   Local fixed arrays
        INTEGER     GSRC  ( MXLOCGRP, MXSPGRP ) ! source IDs for local groups
        INTEGER     GCNT  ( MXLOCGRP )  ! local group count
        INTEGER     LOCGID( MXLOCGRP )  ! local group ID
        REAL        G_DM  ( MXLOCGRP )  ! local group diameter
        REAL        G_FL  ( MXLOCGRP )  ! local group exit flow rate
        REAL        G_HT  ( MXLOCGRP )  ! local group stack height
        REAL        G_TK  ( MXLOCGRP )  ! local group exit temperature
        REAL        G_VE  ( MXLOCGRP )  ! local group exit velocity

C...........   OTHER LOCAL VARIABLES and their descriptions:

        INTEGER     G, I, J, L2, S, S2      ! counters and indices

        INTEGER     FIP              ! tmp FIPS code
        INTEGER     GLM              ! local G max value
        INTEGER     IOS              ! i/o status
        INTEGER     LOCG             ! local G
        INTEGER     MXGCNT           ! max of all GCNT entries
        INTEGER  :: NLOCGRP = 1      ! number of local (facility) groups
        INTEGER     PFIP             ! previous FIPS code
        INTEGER     PRVG             ! G for previous interation

        REAL        DM               ! tmp stack diameter [m]
        REAL        FL               ! tmp stack exit flow [m]
        REAL        HT               ! tmp stack height [m]
        REAL        TK               ! tmp stack exit temperature [K]
        REAL        VE               ! tmp stack exit velocity [m/s]
        REAL        W1, W2           ! tmp multipler weights

        LOGICAL :: LANYGRP = .FALSE.  ! true: one or more groups in inventory
        LOGICAL :: LGROUP  = .FALSE.  ! true: group found for this source
        LOGICAL :: LGRPALL = .FALSE.  ! true: group found for previous facility
        LOGICAL :: STATUS  = .FALSE.  ! true: tmp evaluation status

        CHARACTER*300   BUFFER       ! src description buffers
        CHARACTER*300   MESG         ! msg buffer

        CHARACTER(LEN=CHRLEN3) PLT   ! tmp facility ID
        CHARACTER(LEN=CHRLEN3) PPLT  ! tmp facility ID

        CHARACTER*16 :: PROGNAME = 'ASGNGRPS'   !  program name

C***********************************************************************
C   begin body of subroutine ASGNGRPS

C.........  Write status message
        MESG = 'Assigning stack groups...'
        CALL M3MSG2( MESG )

C.........  Allocate temporary arrays to store group stack parameters
        ALLOCATE( TGRPCNT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TGRPCNT', PROGNAME )
        ALLOCATE( TGRPHT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TGRPHT', PROGNAME )
        ALLOCATE( TGRPDM( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TGRPDM', PROGNAME )
        ALLOCATE( TGRPTK( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TGRPTK', PROGNAME ) 
        ALLOCATE( TGRPVE( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TGRPVE', PROGNAME )
        ALLOCATE( TGRPFL( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TGRPFL', PROGNAME )
        TGRPCNT = 0
        TGRPHT = 0.
        TGRPDM = 0.
        TGRPTK = 0.
        TGRPVE = 0.
        TGRPFL = 0.

C.........  Allocate temporary arrays for interpreting formulas
        ALLOCATE( VALS( NSP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VALS', PROGNAME )
        ALLOCATE( CHRS( NSP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHRS', PROGNAME )
        ALLOCATE( REFS( NSP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'REFS', PROGNAME )
        ALLOCATE( GPSTAT( NSP, NCRIT, MXCHK ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GPSTAT', PROGNAME )
        ALLOCATE( COMCHRS( NSP, NCRIT, MXCHK ), STAT=IOS )
        CALL CHECKMEM( IOS, 'COMCHRS', PROGNAME )  
        CHRS    = ' '   ! array
        COMCHRS = ' '   ! array

C.........  Loop through sources to establish groups
        G    = 0
        PFIP = 0
        PPLT = ' '
        DO S = 1, NSRC

            FIP = IFIP  ( S )
            PLT = CSOURC( S )( SC_BEGP( 2 ):SC_ENDP( 2 ) )
            HT  = STKHT ( S )
            DM  = STKDM ( S )
            TK  = STKTK ( S )
            VE  = STKVE ( S )
            FL  = 0.25 * PI * DM * DM * VE

            GINDEX( S ) = S   ! Set sorting index

C.............  For a new facility, assume that a group will form, but store
C               source ID to reset group if none materialized.
            IF( FIP .NE. PFIP .OR. PLT .NE. PPLT ) THEN

C.................  If there was one or more groups for the previous facility,
C                   compute the group numbers and store info by source
                IF( LGRPALL ) THEN

                    GLM = 0
                    DO I = 1, NLOCGRP

C.........................  If there are multiple stacks in the local group,
C                           then it's a global group.  Update global arrays.
                        IF ( GCNT( I ) .GT. 1 ) THEN

                            DO J = 1, GCNT( I )

                                S2 = GSRC( I,J )
                                GROUPID( S2 ) = G + LOCGID( I )
                                TGRPCNT( S2 ) = GCNT( I )
                                TGRPHT ( S2 ) = G_HT( I )
                                TGRPDM ( S2 ) = G_DM( I )
                                TGRPTK ( S2 ) = G_TK( I )
                                TGRPVE ( S2 ) = G_VE( I )
                                TGRPFL ( S2 ) = G_FL( I )

                                GLM = MAX( GLM, GROUPID( S2 ) )  ! Max for facility

                            END DO

                        END IF

                    END DO

C.....................  Update global counter. When get to the end, we'll have
C                       to reset the groups so that there aren't any holes.
                    G = GLM

                END IF

C.................  Initialize info for possible current group
                PFIP    = FIP
                PPLT    = PLT
                LGRPALL = .FALSE.
                NLOCGRP = 1
                GSRC  ( 1,1 ) = S
                GCNT  ( 1 ) = 1
                LOCGID( 1 ) = 1
                G_HT  ( 1 ) = HT
                G_DM  ( 1 ) = DM
                G_TK  ( 1 ) = TK
                G_VE  ( 1 ) = VE
                G_FL  ( 1 ) = FL

C.............  For the same facility, compare stack parameters to the
C               group stack parameters using the tolerances.  
C.............  If there is a match in stack parameters, set flag current source
C               with the group number, 
            ELSE

                VALS( HT_IDX ) = HT
                VALS( DM_IDX ) = DM
                VALS( TK_IDX ) = TK
                VALS( VE_IDX ) = VE
                VALS( FL_IDX ) = FL

C.................  Compare this source to all other groups/sources for this
C                   plant
                LGROUP = .FALSE.
                DO I = 1, NLOCGRP

                    IF( I .LE. MXLOCGRP ) THEN
                        REFS( HT_IDX ) = G_HT( I )  
                        REFS( DM_IDX ) = G_DM( I )
                        REFS( TK_IDX ) = G_TK( I )
                        REFS( VE_IDX ) = G_VE( I )
                        REFS( FL_IDX ) = G_FL( I )
                    END IF

C.....................  Check tolerances. Use REFS for RANK field since RANK
C                       can't be used to group stacks (it makes no sense)
                    STATUS = EVALCRIT( NSP, NCRIT, MXCHK, VALS, REFS,
     &                                 REFS, CHRS, CRITVALS, COMCHRS,
     &                                 CRITYPES, GPSTAT)

C.....................  If stack parameters meet the criteria, then recompute
C                       the group stack parameters as the weighted average.
C.....................  Exit from loop if a match has been found.
                    IF( STATUS ) THEN

                        GCNT( I ) = GCNT( I ) + 1                
                        W1 = REAL( GCNT( I ) - 1 ) / REAL( GCNT( I ) )
                        W2 = 1. / REAL( GCNT( I ) )
                        G_HT( I ) = G_HT( I ) * W1 + HT * W2
                        G_DM( I ) = G_DM( I ) * W1 + DM * W2
                        G_TK( I ) = G_TK( I ) * W1 + TK * W2
                        G_VE( I ) = G_VE( I ) * W1 + VE * W2
                        G_FL( I ) = G_FL( I ) * W1 + FL * W2

                        IF( GCNT( I ) .LE. MXSPGRP ) THEN
                            GSRC( I,GCNT( I ) ) = S                            
                        END IF

                        MXGCNT = MAX( MXGCNT,GCNT( I ) )

                        LGROUP  = .TRUE.
                        LGRPALL = .TRUE.
                        LANYGRP = .TRUE.
                        EXIT

                    END IF  ! If stack parms are meeting group criteria

                END DO  ! loop through temporary groups at the plant

C.................  If a match was not found for this source, add it's stack
C                   parameters to the temporary group list
                IF( .NOT. LGROUP ) THEN
                    NLOCGRP = NLOCGRP + 1

                    IF ( NLOCGRP .LE. MXLOCGRP ) THEN
                        G_HT  ( NLOCGRP ) = HT
                        G_DM  ( NLOCGRP ) = DM
                        G_TK  ( NLOCGRP ) = TK
                        G_VE  ( NLOCGRP ) = VE
                        G_FL  ( NLOCGRP ) = FL
                        GCNT  ( NLOCGRP ) = 1
                        LOCGID( NLOCGRP ) = NLOCGRP

                        IF( GCNT( NLOCGRP ) .LE. MXSPGRP ) THEN
                            GSRC( NLOCGRP,GCNT( NLOCGRP ) ) = S                            
                        END IF

                        MXGCNT = MAX( MXGCNT,GCNT( NLOCGRP ) )

                    END IF

                END IF

            END IF      ! If same facility or not

        END DO          ! loop through sources

C.........  Ensure that we didn't run out of space for the per-facility groups
        IF ( NLOCGRP .GT. MXLOCGRP ) THEN

            WRITE( MESG,94010 ) 'INTERNAL ERROR: The maximum number ' //
     &             'of groups per plant (', MXLOCGRP, ') was ' //
     &             'exceeded with a number of ', NLOCGRP
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Ensure that we had enough space for the sources per group
        IF( MXGCNT .GT. MXSPGRP ) THEN
            WRITE( MESG,94010 ) 'INTERNAL ERROR: The maximum number ' //
     &             'of sources per group (', MXSPGRP, ') was ' //
     &             'exceeded with a number of ', MXGCNT
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Postprocess group IDs to remove gaps in number and to count the 
C           number of inventory groups.  If LANYGRP is false, then GROUPID
C           will be zeros, so do not want to sort.  GINDEX will be maintained
C           as the source number.
        IF( LANYGRP ) THEN
            CALL SORTI1( NSRC, GINDEX, GROUPID )
        ELSE
            MESG = 'NOTE: No stack groups assigned.'
            CALL M3MSG2( MESG )
        END IF

        PRVG = 0
        G    = 0
        DO S = 1, NSRC

            I    = GINDEX( S )
            LOCG = GROUPID( I )   

            IF ( LOCG .NE. 0 ) THEN

C.................  If a new group, increment group ID
                IF ( LOCG .NE. PRVG ) G = G + 1

C.................  Store global group number
                GROUPID( I ) = G
                PRVG = LOCG

            END IF

        END DO

C.........  Store the actual number of groups present
        NINVGRP = G

C.........  Allocate memory for group information and populate group arrays
        ALLOCATE( GRPGID( NINVGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPGID', PROGNAME )
        ALLOCATE( GRPLAT( NINVGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPLAT', PROGNAME )
        ALLOCATE( GRPLON( NINVGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPLON', PROGNAME )
        ALLOCATE( GRPDM( NINVGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPDM', PROGNAME )
        ALLOCATE( GRPHT( NINVGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPHT', PROGNAME )
        ALLOCATE( GRPTK ( NINVGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPTK', PROGNAME )
        ALLOCATE( GRPVE( NINVGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPVE', PROGNAME )
        ALLOCATE( GRPFL( NINVGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPFL', PROGNAME )
        ALLOCATE( GRPCNT( NINVGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPCNT', PROGNAME )

        GRPGID  = 0
        GRPLAT  = BADVAL3
        GRPLON  = BADVAL3
        GRPDM   = BADVAL3
        GRPHT   = BADVAL3
        GRPTK   = BADVAL3
        GRPVE   = BADVAL3
        GRPFL   = BADVAL3
        GRPCNT  = 0

C.........  Store the group information more succinctly
        PRVG = 0
        DO S = 1, NSRC
            I = GINDEX( S )
            G = GROUPID( I )

C.............  Skip sources that aren't in a group
            IF ( G .EQ. 0 ) CYCLE

C.............  Store. Location will be from first source encountered in group
            IF ( G .NE. PRVG ) THEN
                GRPGID( G ) = G
                GRPLAT( G ) = YLOCA  ( I )
                GRPLON( G ) = XLOCA  ( I )
                GRPDM ( G ) = TGRPDM ( I )
                GRPHT ( G ) = TGRPHT ( I )
                GRPTK ( G ) = TGRPTK ( I )
                GRPVE ( G ) = TGRPVE ( I )
                GRPFL ( G ) = TGRPFL ( I )
                GRPCNT( G ) = TGRPCNT( I )

C.............  Otherwise, give warnings if stacks in the same group do not have
C               the same stack locations.
            ELSE 

                CALL FMTCSRC( CSOURC( I ), NCHARS, BUFFER, L2 )

                IF ( FLTERR( GRPLAT( G ), YLOCA( I ) ) ) THEN
                    WRITE( MESG,94078 ) 'WARNING: Source latitude',
     &                     YLOCA( I ), 'is inconsistent with group' // 
     &                     CRLF() // BLANK10 // 'latitude', GRPLAT(G),
     &                     'taken from the first source in the group, '
     &                     // 'for source:' // CRLF() // BLANK10 //
     &                     BUFFER( 1:L2 ) 
                    CALL M3MESG( MESG )
                END IF

                IF ( FLTERR( GRPLON( G ), XLOCA( I ) ) ) THEN
                    WRITE( MESG,94078 ) 'WARNING: Source longitude',
     &                     XLOCA( I ), 'is inconsistent with group' // 
     &                     CRLF() // BLANK10 // 'longitude', GRPLON(G),
     &                     'taken from the first source in the group, '
     &                     // 'for source:' // CRLF() // BLANK10 //
     &                     BUFFER( 1:L2 ) 
                    CALL M3MESG( MESG )
                END IF

            END IF

            PRVG = G
                        
        END DO

C.........  Deallocate local memory
        DEALLOCATE( TGRPCNT, TGRPHT, TGRPDM, TGRPTK, 
     &              TGRPVE, TGRPFL, GPSTAT )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 2X  ) )

94078   FORMAT( A, 1X, F6.2, 1X, A, 1X, F6.2, 1X, A )

        END SUBROUTINE ASGNGRPS

