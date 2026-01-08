
        SUBROUTINE POSTPLM( EMLAYS, S, ZBOT, ZTOP, PRESF, ZZF, TA, ZH, 
     &                      LBOT, LTOP, LFRAC )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C    Subroutine POSTPLM computes plume fractions given a top and bottom
C    height of the plume.  It assumes a uniform distribution in pressure
C    (mass concentration -- minor hydrostatic assumption) from bottom to top.
C
C  PRECONDITIONS REQUIRED:
C    Top and bottom of plume as input, vertical grid structure defined, 
C    vertical pressure distribution and temperature provided.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       I/O API
C
C  REVISION  HISTORY:
C    Copied from postplm.f v 1.3 in DAQM-V2 Emissions Preprocessor by
C           M. Houyoux 3/99
C    Replaced POLY calls with direct hydrostatic calculations by
C           G. Pouliot 2/05
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
C***********************************************************************
 
        USE M3UTILIO

        IMPLICIT NONE
 
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'
C        INCLUDE 'PARMS3.EXT'
C        INCLUDE 'IODECL3.EXT'
C        INCLUDE 'FDESC3.EXT'
        INCLUDE 'CONST3.EXT'

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)  CRLF

C        EXTERNAL      CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: EMLAYS           ! no. emissions layers
        INTEGER, INTENT (IN) :: S                ! source ID
        REAL   , INTENT (IN) :: ZBOT             ! plume bottom elevation (m)
        REAL   , INTENT (IN) :: ZTOP             ! plume top elevation (m)
        REAL   , INTENT (IN) :: PRESF( 0:EMLAYS )! pressure at full-levels (mb)
        REAL   , INTENT (IN) :: ZZF  ( 0:EMLAYS )! elevation at full-levels (m)
        REAL   , INTENT (IN) :: TA   ( 1:EMLAYS )! temperature at half-levels (K)
        REAL   , INTENT (IN) :: ZH   ( 1:EMLAYS )! layer center  height (m)
        INTEGER, INTENT(OUT) :: LBOT             ! plume bottom layer
        INTEGER, INTENT(OUT) :: LTOP             ! plume top layer
        REAL   , INTENT(OUT) :: LFRAC( EMLAYS )  ! layer fractions for source

C...........   Local variables

        INTEGER       L

        DOUBLE PRECISION    DDP
        DOUBLE PRECISION    PDIFF
        
        REAL          PBOT, PTOP
        REAL          TEMP

        CHARACTER(300) MESG

C***********************************************************************
C   begin body of subroutine POSTPLM

C...........   Compute LBOT, LTOP so that
C...........   ZZF( LBOT-1 ) <= ZBOT < ZZF( LBOT ) and
C...........   ZZF( LTOP-1 ) <= ZTOP < ZZF( LTOP )
 
        DO L = 1, EMLAYS - 1
            IF ( ZBOT <= ZZF( L ) ) THEN
                LBOT = L
                GO TO  122   ! end loop and skip reset of LBOT
            ELSE
                LFRAC( L ) = 0.0      ! fractions below plume
            END IF
        END DO
        LBOT = EMLAYS           !  fallback

122     CONTINUE                !  loop exit:  bottom found at LBOT
 
        IF ( ZTOP <= ZZF( LBOT ) ) THEN  !  plume in this layer
 
            LFRAC( LBOT ) = 1.0
            LTOP = LBOT
 
            DO L = LBOT + 1, EMLAYS  ! fractions above plume
                LFRAC( L ) = 0.0
            END DO
 
C.........  Note- this check not in original algorithm, but w/o it,
C                         can end up with fractions > 1.0
        ELSE IF( LBOT == EMLAYS ) THEN    ! plume above top layer
 
            LFRAC( LBOT ) = 1.0
 
            DO L = 1, EMLAYS-1       ! fractions below plume
                LFRAC( L ) = 0.0
            END DO
 
        ELSE                               ! plume crosses layers
 
            DO L = LBOT + 1, EMLAYS
                IF ( ZTOP <= ZZF( L ) ) THEN
                    LTOP = L
                    GO TO 126  ! end loop and skip reset of LTOP
                END IF
            END DO
            LTOP = EMLAYS

126         CONTINUE
 
C...........   Compute corresponding PBOT,PTOP so that
C...........   PRESF( LBOT-1 ) <= PBOT < PRESF( LBOT ) and
C...........   PRESF( LTOP-1 ) <= PTOP < PRESF( LTOP )

C............. If above full layer and below half layer... 
            IF( ZBOT < ZH( LBOT ) .AND. ZBOT > ZZF( LBOT-1 ) ) THEN

C.................  Special case near ground
                IF( ZBOT < ZH( 1 ) ) THEN
                    TEMP = TA( 1 )
                ELSE
                    TEMP = ( TA( LBOT ) + TA( LBOT-1 ) ) / 2.
                END IF

C.............  Otherwise, above full layer and above half layer
            ELSE
                TEMP = ( TA( LBOT ) + TA( LBOT+1 ) ) / 2.
            END IF

C.............  Calculate bottom using hydrostatic assumption            
            PBOT = PRESF( LBOT ) * 
     &             EXP( GRAV / (RDGAS*TEMP) * (ZZF( LBOT ) - ZBOT ))
            
C.............  If above full layer and below half layer... 
            IF( ZTOP < ZH( LTOP ) .AND. ZTOP > ZZF( LTOP-1 ) ) THEN

C.................  Special case near ground
                IF( ZTOP < ZH( 1 ) ) THEN
                    TEMP = TA( 1 )
                ELSE
                    TEMP = ( TA( LTOP ) + TA( LTOP-1 ) ) / 2.
                END IF

C.............  Otherwise, above full layer and above half layer
            ELSE
                TEMP = ( TA( LTOP ) + TA( LTOP+1 ) ) / 2.
            END IF

C.............  Calculate top using hydrostatic assumption            
            PTOP = PRESF( LTOP-1 ) * 
     &             EXP( -GRAV / (RDGAS*TEMP) * (ZTOP - ZZF( LTOP-1 )) )
            
            PDIFF = DBLE( PBOT ) - DBLE( PTOP )
            
            IF( PDIFF > 0. ) THEN
            
                DDP = 1.0D0 / ( PDIFF )
                LFRAC( LBOT ) = DDP * 
     &                          (DBLE( PBOT ) - DBLE( PRESF( LBOT ) ))
                LFRAC( LTOP ) = DDP * 
     &                          (DBLE( PRESF( LTOP-1 ) ) - DBLE( PTOP ))
 
            ELSE
                WRITE( MESG,94010 )
     &           'Infinitely small plume created for source ', S,
     &           CRLF() // BLANK5 // 
     &           'because of inverted vertical pressure gradient!' //
     &           CRLF() // BLANK5 // 
     &           'All emissions put in first layer.'
                CALL M3WARN( 'POSTPLM', 0,0, MESG )

                LBOT = 1
                LTOP = 1
                LFRAC( LBOT ) = 1.0
 
            ENDIF
 
            DO L = LBOT+1, LTOP-1 !  layers in plume
                LFRAC( L ) = DDP * 
     &                       (DBLE( PRESF( L-1 ) ) - DBLE( PRESF( L ) ))
            END DO
 
            DO L = LTOP+1, EMLAYS !  fractions above plume
                LFRAC( L ) = 0.0
            END DO
 
        END IF          !  if ztop in same layer as zbot, or not

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

        END SUBROUTINE POSTPLM
