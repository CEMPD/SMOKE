
        SUBROUTINE MRGELEV( NSRC, NMAJOR, NPING, 
     &                      KEY1, KEY2, KEY3, KEY4, CNV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine multiplies a source-emissions vector with optionally 
C      a speciation array and multiplicative control array. An additive 
C      control array can be added to the emissions.  The first time this
C      routine is called, a PinG- and elevated-source-specific set of arrays are
C      allocated for storing and processing the PinG and elevated emissions.
C
C  PRECONDITIONS REQUIRED:
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
C.........  This module contains the major data structure and control flags
        USE MODMERGE

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS
        INTEGER    FIND1
        EXTERNAL   FIND1

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC        ! number of source
        INTEGER     , INTENT (IN) :: NMAJOR      ! no. elevated sources
        INTEGER     , INTENT (IN) :: NPING       ! no. plume-in-grid sources
        INTEGER     , INTENT (IN) :: KEY1        ! inven emissions index
        INTEGER     , INTENT (IN) :: KEY2        ! mult controls index
        INTEGER     , INTENT (IN) :: KEY3        ! additive controls index
        INTEGER     , INTENT (IN) :: KEY4        ! speciation index
        REAL        , INTENT (IN) :: CNV         ! units conversion factor

C.........  Other local variables
        INTEGER         I, J, K, L, S   ! counters and indicies
        INTEGER         IDX             ! index to list of counties in grid   
        INTEGER         IOS             ! i/o status

        REAL*8          SUM1            ! sum for GOUT1   
        REAL*8          SUM2            ! sum for GOUT2 
        REAL*8          ADD             ! tmp value with additive controls
        REAL*8          MULT            ! tmp value with multiplictv controls
        REAL*8          REAC            ! tmp value with reactivity controls
        REAL*8          VAL             ! tmp value  
        REAL*8          VMP             ! tmp market penetration value  

        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine called

        CHARACTER*300   MESG            ! message buffer

        CHARACTER*16 :: PROGNAME = 'MRGELEV' ! program name

C***********************************************************************
C   begin body of subroutine MRGELEV

C.........  For the first time this routine is called, process the plume-in-
C           grid indicator array to allocate and generate the necessary
C           indices
        IF( FIRSTIME ) THEN

C.............  Allocate indices used for processing in this routine
            ALLOCATE( ELEVSIDX( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ELEVSIDX', PROGNAME )
            ALLOCATE( PINGGIDX( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PINGGIDX', PROGNAME )
            ELEVSIDX = 0   ! array
            PINGGIDX = 0   ! array

C.............  Allocate arrays for storing emissions (arrays from MODMERGE)
            ALLOCATE( ELEVEMIS( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ELEVEMIS', PROGNAME )
            ALLOCATE( PGRPEMIS( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PGRPEMIS', PROGNAME )

C.............  Create indices used for processing in this routine
            I = 0
            J = 0
            DO S = 1, NSRC

C.................  Set elevated sources index
                IF( LMAJOR( S ) ) THEN
                    I = I + 1
                    IF( I .GT. NMAJOR ) CYCLE

C.....................  Store index to sources 
                    ELEVSIDX( I ) = S

                END IF

C.................  Set PinG indicies
                IF( LPING( S ) ) THEN
                    J = J + 1
                    IF( J .GT. NMAJOR ) CYCLE

C.....................  Find group ID in list
                    K = FIND1( GROUPID( S ), NGROUP, GRPGID )

C.....................  Store index to sources and index to groups
                    PINGGIDX( J ) = K

                END IF

            END DO

C.............  Abort if dimensions incorrect
            IF( I .NE. NMAJOR ) THEN
                WRITE( MESG,94010 ) 'INTERNAL ERROR: Number of ' //
     &                  'elevated sources I=', I, 
     &                  'unequal to dimension NMAJOR=', NMAJOR
                CALL M3MSG2( MESG ) 
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF

            IF( J .NE. NPING ) THEN
                WRITE( MESG,94010 ) 'INTERNAL ERROR: Number of ' //
     &                  'plume-in-grid sources J=', J, 
     &                  'unequal to dimension NPING=', NPING
                CALL M3MSG2( MESG ) 
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF

            FIRSTIME = .FALSE.

        END IF  ! end of firstime section

C.........  Initialize emissions values
        PGRPEMIS = 0  ! array
        ELEVEMIS = 0  ! array

C.........  Check if this is a valid inventory pollutant for this call
        IF( KEY1 .GT. 0 ) THEN

C............. If multiplicative controls, additive controls, and speciation
            IF( KEY2 .GT. 0 .AND. KEY3 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                DO K = 1, NMAJOR

                    S   = ELEVSIDX( K )   ! index to source arrays

                    VAL  = PEMSRC ( S,KEY1 ) * PSMATX( S,KEY4 ) 
                    MULT = VAL * PCUMATX( S,KEY2 )

                    ADD = PCAMATX( S,KEY3 ) * PSMATX( S,KEY4 )

                    VAL  = ADD + MULT
                    VMP  = PRINFO( S,2 )
                    VAL = ( VAL * (1.-VMP) + PRINFO( S,1 ) * VMP )

                    IDX = PINGGIDX( K )   ! index to group arrays
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL

                    ELEVEMIS( K ) = ELEVEMIS( K ) + VAL * CNV                   

                END DO

C............. If multiplicative controls & additive controls
            ELSE IF( KEY2 .GT. 0 .AND. KEY3 .GT. 0 ) THEN

                DO K = 1, NMAJOR

                    S   = ELEVSIDX( K )   ! index to source arrays

                    MULT = PEMSRC ( S,KEY1 ) * PCUMATX( S,KEY2 )
                    ADD  = PCAMATX( S,KEY3 )
                    VAL  = ADD + MULT

                    IDX = PINGGIDX( K ) 
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL

                    ELEVEMIS( K ) = ELEVEMIS( K ) + VAL * CNV                    

                END DO

C............. If multiplicative controls & speciation
            ELSE IF( KEY2 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                DO K = 1, NMAJOR

                    S   = ELEVSIDX( K )   ! index to source arrays

                    VAL  = PEMSRC ( S,KEY1 ) * PSMATX( S,KEY4 ) * CNV
                    MULT = VAL * PCUMATX( S,KEY2 )

                    VMP  = PRINFO( S,2 )
                    VAL = ( MULT * (1.-VMP) + PRINFO( S,1 ) * VMP )

                    IDX = PINGGIDX( K )
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL

                    ELEVEMIS( K ) = ELEVEMIS( K ) + VAL * CNV

                END DO

C............. If additive controls & speciation
            ELSE IF( KEY3 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

                DO K = 1, NMAJOR

                    S   = ELEVSIDX( K )   ! index to source arrays

                    VAL = PEMSRC ( S,KEY1 ) * PSMATX( S,KEY4 ) 
                    ADD = PCAMATX( S,KEY3 ) * PSMATX( S,KEY4 )

                    VAL = ADD + VAL
                    VMP  = PRINFO( S,2 )
                    VAL = ( VAL * (1.-VMP) + PRINFO( S,1 ) * VMP )

                    IDX = PINGGIDX( K )
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL

                    ELEVEMIS( K ) = ELEVEMIS( K ) + VAL * CNV

                END DO

C............. If multiplicative controls only
            ELSE IF( KEY2 .GT. 0 ) THEN

                DO K = 1, NMAJOR

                    S   = ELEVSIDX( K )   ! index to source arrays

                    VAL = PEMSRC( S,KEY1 ) * PCUMATX( S,KEY2 )

                    IDX = PINGGIDX( K )
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL

                    ELEVEMIS( K ) = ELEVEMIS( K ) + VAL * CNV

                END DO

C............. If additive controls only
            ELSE IF( KEY3 .GT. 0 ) THEN

                DO K = 1, NMAJOR

                    S   = ELEVSIDX( K )   ! index to source arrays

                    VAL = PEMSRC( S,KEY1 ) + PCAMATX( S,KEY3 )

                    IDX = PINGGIDX( K )
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL

                    ELEVEMIS( K ) = ELEVEMIS( K ) + VAL * CNV
                END DO

C.............  If speciation only
            ELSE IF( KEY4 .GT. 0 ) THEN

                DO K = 1, NMAJOR

                    S   = ELEVSIDX( K )   ! index to source arrays

                    VAL  = PEMSRC ( S,KEY1 ) * PSMATX( S,KEY4 ) 
                    VMP  = PRINFO( S,2 )
                    VAL = ( VAL * (1.-VMP) + PRINFO( S,1 ) * VMP )

                    IDX = PINGGIDX( K )
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL

                    ELEVEMIS( K ) = ELEVEMIS( K ) + VAL * CNV

                END DO

C.............  If inventory pollutant only
            ELSE 

                DO K = 1, NMAJOR

                    S   = ELEVSIDX( K )   ! index to source arrays

                    VAL = PEMSRC( S,KEY1 )

                    IDX = PINGGIDX( K )
                    IF( IDX .GT. 0 ) 
     &                  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL

                    ELEVEMIS( K ) = ELEVEMIS( K ) + VAL * CNV

                END DO

            END IF  ! End which of controls and speciation

        END IF      ! End if no inventory emissions

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

94000   FORMAT( A )

94010   FORMAT( 10 ( A, :, I8, :, 2X ) )

        END SUBROUTINE MRGELEV
