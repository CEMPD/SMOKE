
        SUBROUTINE POSTPLM( EMLAYS, S, PSFC, ZBOT, ZTOP, PRES, ZZF, ZH, 
     &                      LTOP, LFRAC )

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
C    vertical pressure distribution provided.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       I/O API
C
C  REVISION  HISTORY:
C    Copied from postplm.f v 1.3 in DAQM-V2 Emissions Preprocessor by
C           M. Houyoux 3/99
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C***********************************************************************
 
        IMPLICIT NONE
 
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'
        INCLUDE 'PARMS3.EXT'
        INCLUDE 'IODECL3.EXT'
        INCLUDE 'FDESC3.EXT'

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2   CRLF
        REAL          POLY

        EXTERNAL      CRLF, POLY

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: EMLAYS           ! no. emissions layers
        INTEGER, INTENT (IN) :: S                ! source ID
        REAL   , INTENT (IN) :: PSFC             ! surface pressure
        REAL   , INTENT (IN) :: ZBOT             ! plume bottom elevation (m)
        REAL   , INTENT (IN) :: ZTOP             ! plume top elevation (m)
        REAL   , INTENT (IN) :: PRES ( EMLAYS )  ! pressure (Pa)
        REAL   , INTENT (IN) :: ZZF( 0:EMLAYS )  ! elevation at full-levels
        REAL   , INTENT (IN) :: ZH   ( EMLAYS )  ! layer center  height (m)
        INTEGER, INTENT(OUT) :: LTOP             ! plume top layer
        REAL   , INTENT(OUT) :: LFRAC( EMLAYS )  ! layer fractions for source

C...........   Local variables

        INTEGER       L, M

        INTEGER       LBOT 

        REAL          DDP
        REAL          PBOT, PTOP
        REAL          PDIFF
        REAL          PRESF( 0:EMLAYS )

        CHARACTER*300 MESG

C***********************************************************************
C   begin body of subroutine POSTPLM

C...........   Compute LBOT, LTOP so that
C...........   ZZF( LBOT-1 ) <= ZBOT < ZZF( LBOT ) and
C...........   ZZF( LTOP-1 ) <= ZTOP < ZZF( LTOP )
 
        DO L = 1, EMLAYS - 1
            IF ( ZBOT .LE. ZZF( L ) ) THEN
                LBOT = L
                GO TO  122   ! end loop and skip reset of LBOT
            ELSE
                LFRAC( L ) = 0.0      ! fractions below plume
            END IF
        END DO
        LBOT = EMLAYS           !  fallback

122     CONTINUE                !  loop exit:  bottom found at LBOT
 
        IF ( ZTOP .LE. ZZF( LBOT ) ) THEN  !  plume in this layer
 
            LFRAC( LBOT ) = 1.0
            LTOP = LBOT
 
            DO L = LBOT + 1, EMLAYS  ! fractions above plume
                LFRAC( L ) = 0.0
            END DO
 
C.........  Note: this check not in original algorithm, but w/o it,
C                         can end up with fractions > 1.0
        ELSEIF( LBOT .EQ. EMLAYS ) THEN    ! plume above top layer
 
            LFRAC( LBOT ) = 1.0
 
            DO L = 1, EMLAYS-1       ! fractions below plume
                LFRAC( L ) = 0.0
            END DO
 
        ELSE                               ! plume crosses layers
 
            DO L = LBOT + 1, EMLAYS
                IF ( ZTOP .LE. ZZF( L ) ) THEN
                    LTOP = L
                    GO TO 126  ! end loop and skip reset of LTOP
                END IF
            END DO
            LTOP = EMLAYS

126         CONTINUE
 
C...........   Compute corresponding PBOT,PTOP so that
C...........   PRESF( LBOT-1 ) <= PBOT < PRESF( LBOT ) and
C...........   PRESF( LTOP-1 ) <= PTOP < PRESF( LTOP )
C...........   (Use 3rd order polynomial via POLY() )
 
            PRESF( 0 ) = PSFC
            PRESF( 1 ) = POLY( ZZF( 1 ), ZH( 1 ), PRES( 1 ), 3 )
            PRESF( 2 ) = POLY( ZZF( 2 ), ZH( 1 ), PRES( 1 ), 3 )
 
            DO L = 3, EMLAYS-2
                M = L - 2
                PRESF( L ) = POLY( ZZF( L ), ZH( M ), PRES( M ), 3 )
            END DO
 
            PRESF( EMLAYS-1 ) = POLY( ZZF ( EMLAYS-1 ),
     &                                ZH  ( EMLAYS-3 ),
     &                                PRES( EMLAYS-3 ), 3 )
            PRESF( EMLAYS   ) = POLY( ZZF ( EMLAYS ),
     &                                ZH  ( EMLAYS-1 ),
     &                                PRES( EMLAYS-1 ), 1 )
 
            M    = MIN( MAX( 0, LBOT-2 ), EMLAYS-3 )
            PBOT = POLY( ZBOT, ZZF( M ), PRESF( M ), 3 )
 
            M    = MIN( MAX( 0, LTOP-2 ), EMLAYS-3 )
            PTOP = POLY( ZTOP, ZZF( M ), PRESF( M ), 3 )
 
            PDIFF = PBOT - PTOP
            IF( PDIFF .GT. 0 ) THEN

               DDP  = 1.0 / ( PDIFF )  !  = d(plumefrac)/d
               LFRAC( LBOT ) = DDP * ( PBOT - PRESF( LBOT ) )
               LFRAC( LTOP ) = DDP * ( PRESF( LTOP-1 ) - PTOP )
 
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
                LFRAC( L ) = DDP*( PRESF( L-1 ) - PRESF( L ) )
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
