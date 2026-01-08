
        SUBROUTINE MRGELEV( NSRC, NMAJOR, NPING, 
     &                      KEY1, KEY2, KEY4, CNV, LINIT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine multiplies a source-emissions vector with optionally 
C      a speciation array and multiplicative control array. 
C      The first time this routine is called, a PinG- and elevated-source-
C      specific set of arrays are allocated for storing and processing the 
C      PinG and elevated emissions.  This routine computes the elevated 
C      emissions for UAM-style processing and the PinG grouped emissions
C      for CMAQ PinG files.
C
C  PRECONDITIONS REQUIRED:
C      Assumes any source that will be a PinG source is also a major (elevated)
C      source.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE M3UTILIO

        USE MODMERGE, ONLY: ELEVFLAG, PEMSRC, PSMATX, PRINFO, INLINEFLAG,
     &                      SRCGRPFLAG

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: PCUMATX

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: ELEVSIDX, PINGGIDX, NGROUP, GRPGID, PGRPEMIS,
     &                     ELEVEMIS, LMAJOR, LPING, GROUPID,
     &                     ELEVGRPID, EMELEVGRP

        IMPLICIT NONE

C...........   INCLUDES
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS
C       INTEGER    FIND1
C        EXTERNAL   FIND1

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC        ! number of source
        INTEGER     , INTENT (IN) :: NMAJOR      ! no. elevated sources
        INTEGER     , INTENT (IN) :: NPING       ! no. plume-in-grid sources
        INTEGER     , INTENT (IN) :: KEY1        ! inven emissions index
        INTEGER     , INTENT (IN) :: KEY2        ! mult controls index
        INTEGER     , INTENT (IN) :: KEY4        ! speciation index
        REAL        , INTENT (IN) :: CNV         ! units conversion factor
        LOGICAL     , INTENT (IN) :: LINIT       ! true: initialize ELEVEMIS

C.........  Other local variables
        INTEGER         I, J, K, L, S   ! counters and indicies
        INTEGER         IDX             ! index to list of counties in grid   
        INTEGER         IOS             ! i/o status

        REAL(8)         SUM1            ! sum for GOUT1   
        REAL(8)         SUM2            ! sum for GOUT2 
        REAL(8)         MULT            ! tmp value with multiplictv controls
        REAL(8)         REAC            ! tmp value with reactivity controls
        REAL(8)         VAL             ! tmp value  
        REAL(8)         VMP             ! tmp market penetration value  

        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine called

        CHARACTER(300)  MESG            ! message buffer

        CHARACTER(16) :: PROGNAME = 'MRGELEV' ! program name

C***********************************************************************
C   begin body of subroutine MRGELEV

C.........  For the first time this routine is called, process the plume-in-
C           grid indicator array to allocate and generate the necessary
C           indices
        IF( FIRSTIME ) THEN

C.............  Allocate indices used for processing in this routine
            ALLOCATE( ELEVSIDX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ELEVSIDX', PROGNAME )
            ALLOCATE( PINGGIDX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PINGGIDX', PROGNAME )
            ELEVSIDX = 0   ! array
            PINGGIDX = 0   ! array

C.............  Allocate arrays for storing emissions (arrays from MODMERGE)
            ALLOCATE( ELEVEMIS( NMAJOR ), STAT=IOS )    ! ungrouped Elev
            CALL CHECKMEM( IOS, 'ELEVEMIS', PROGNAME )
            ALLOCATE( PGRPEMIS( NGROUP ), STAT=IOS )    ! grouped PinG
            CALL CHECKMEM( IOS, 'PGRPEMIS', PROGNAME )

C.............  Create indices used for processing in this routine
            I = 0
            J = 0
            DO S = 1, NSRC

C.................  Set elevated sources index
                IF( ELEVFLAG .AND. LMAJOR( S ) ) THEN
                    I = I + 1
                    IF( I .GT. NMAJOR ) CYCLE

C.....................  Store index to major count 
                    ELEVSIDX( S ) = I

                END IF

C.................  Set elevated sources index
                IF( INLINEFLAG .AND. LMAJOR( S ) .AND. (.NOT. LPING(S)) ) THEN

                    K = FIND1( GROUPID( S ), NGROUP, GRPGID )
                    IF( K .GT. 0 ) ELEVSIDX( S ) = K
                END IF
		
C.................  Set PinG indicies
                IF( LPING( S ) ) THEN

C.....................  Find group ID in list
                    K = FIND1( GROUPID( S ), NGROUP, GRPGID )

C.....................  Store index to groups
                    IF( K .GT. 0 ) PINGGIDX( S ) = K

                END IF

            END DO

C.............  Abort if dimensions incorrect
            IF( ELEVFLAG .AND. I .NE. NMAJOR ) THEN
                WRITE( MESG,94010 ) 'INTERNAL ERROR: Number of ' //
     &                  'elevated sources I=', I, 
     &                  'unequal to dimension NMAJOR=', NMAJOR
                CALL M3MSG2( MESG ) 
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF


	    
            FIRSTIME = .FALSE.

        END IF  ! end of firstime section

C.........  Initialize emissions values if calling routine indicates
C           that this is a new species being processed.
        IF ( LINIT ) THEN
            PGRPEMIS = 0  ! array
            ELEVEMIS = 0  ! array
            IF( SRCGRPFLAG ) EMELEVGRP = 0.
        END IF

C.........  Check if this is a valid inventory pollutant for this call
        IF( KEY1 .GT. 0 ) THEN

C............. If multiplicative controls & speciation
            IF( KEY2 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                DO S = 1, NSRC

C.....................  Skip if source is not elevated and not PinG
                    IF( ELEVSIDX( S ) .EQ. 0 .AND.
     &                  PINGGIDX( S ) .EQ. 0       ) CYCLE

                    VAL  = PEMSRC ( S,KEY1 ) * PSMATX( S,KEY4 ) ! Spec value
                    MULT = VAL * PCUMATX( S,KEY2 )              ! w/ control

                    VMP  = PRINFO( S,2 )
                    VAL = ( MULT * (1.-VMP) + PRINFO( S,1 ) * VMP ) ! w/ reactivity control

C......................  NOTE - apply units conversion after application
C                        of reactivity matrix.
                    IDX = PINGGIDX( S )
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL * CNV

                    IDX = ELEVSIDX( S )
                    IF( IDX .GT. 0 )
     &                  ELEVEMIS( IDX ) = ELEVEMIS( IDX ) + VAL * CNV

                    IF( SRCGRPFLAG ) THEN
                        IDX = ELEVGRPID( S )
                        EMELEVGRP( IDX ) = EMELEVGRP( IDX ) + VAL * CNV
                    END IF

                END DO

C............. If multiplicative controls only
            ELSE IF( KEY2 .GT. 0 ) THEN

                DO S = 1, NSRC

                    IF( ELEVSIDX( S ) .EQ. 0 .AND.
     &                  PINGGIDX( S ) .EQ. 0       ) CYCLE

                    VAL = PEMSRC( S,KEY1 ) * PCUMATX( S,KEY2 )

                    IDX = PINGGIDX( S )
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL * CNV

                    IDX = ELEVSIDX( S )
                    IF( IDX .GT. 0 )
     &                  ELEVEMIS( IDX ) = ELEVEMIS( IDX ) + VAL * CNV

                    IF( SRCGRPFLAG ) THEN
                        IDX = ELEVGRPID( S )
                        EMELEVGRP( IDX ) = EMELEVGRP( IDX ) + VAL * CNV
                    END IF

                END DO

C.............  If speciation only
            ELSE IF( KEY4 .GT. 0 ) THEN

                DO S = 1, NSRC

                    IF( ELEVSIDX( S ) .EQ. 0 .AND.
     &                  PINGGIDX( S ) .EQ. 0       ) CYCLE

                    VAL  = PEMSRC ( S,KEY1 ) * PSMATX( S,KEY4 ) 
                    VMP  = PRINFO( S,2 )
                    VAL = ( VAL * (1.-VMP) + PRINFO( S,1 ) * VMP )

                    IDX = PINGGIDX( S )
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL * CNV

                    IDX = ELEVSIDX( S )
                    IF( IDX .GT. 0 )
     &                  ELEVEMIS( IDX ) = ELEVEMIS( IDX ) + VAL * CNV

                    IF( SRCGRPFLAG ) THEN
                        IDX = ELEVGRPID( S )
                        EMELEVGRP( IDX ) = EMELEVGRP( IDX ) + VAL * CNV
                    END IF

                END DO

C.............  If inventory pollutant only
            ELSE 

                DO S = 1, NSRC

                    IF( ELEVSIDX( S ) .EQ. 0 .AND.
     &                  PINGGIDX( S ) .EQ. 0       ) CYCLE

                    VAL = PEMSRC( S,KEY1 )

                    IDX = PINGGIDX( S )
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL * CNV

                    IDX = ELEVSIDX( S )
                    IF( IDX .GT. 0 )
     &                  ELEVEMIS( IDX ) = ELEVEMIS( IDX ) + VAL * CNV

                    IF( SRCGRPFLAG ) THEN
                        IDX = ELEVGRPID( S )
                        EMELEVGRP( IDX ) = EMELEVGRP( IDX ) + VAL * CNV
                    END IF

                END DO

            END IF  ! End which of controls and speciation

        END IF      ! End if no inventory emissions

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

94000   FORMAT( A )

94010   FORMAT( 10 ( A, :, I8, :, 2X ) )

        END SUBROUTINE MRGELEV
