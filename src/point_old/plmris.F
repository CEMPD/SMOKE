
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  PLMRIS( ZH, ZF, WSPD, TA, DTHDZ,
     &                      USTAR, HMIX, HFX,
     &                      STKHT, STKDM, STKTK, STKVE, 
     &                      ZSTK, TSTK, WSTK, LSTK, LPBL, 
     &                      TOP, BOT )

C***********************************************************************
C  subroutine body starts at line  138
C
C  FUNCTION:  
C       computes elevation TOP, BOT of plume top and bottom.
C
C  PRECONDITIONS REQUIRED:
C	meteorology and stack parameters
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC, based on Briggs algorithm adapted from
C       RADM 2.6 subroutine PLUMER() (but with completely different 
C       data structuring).
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'CONST3.EXT'    ! physical and mathematical constants
        INCLUDE 'PARMS3.EXT'    ! I/O API constants
        INCLUDE 'FDESC3.EXT'    ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'LADIMS3.EXT'   !  emis temporal modeling parms


C...........   ARGUMENTS and their descriptions:


        REAL    ZH   ( NLAYS )  !  layer surface height (m)
        REAL    ZF   ( NLAYS )  !  layer center  height (m)
        REAL    WSPD ( NLAYS )  !  wind speed (m/s)
        REAL    TA   ( NLAYS )  !  temperature (deg K)
        REAL    DTHDZ( NLAYS )  !  gradient of THETV
        REAL    USTAR           !  friction velocity (m/s)
        REAL    STKHT           !  stack height (m)
        REAL    STKDM           !  stack diameter (m)
        REAL    STKTK           !  exhaust temperature (deg K)
        REAL    STKVE           !  exhaust velocity (m/s)
        REAL    HMIX            !  mixing height (m)
        REAL    HFX             !  sensible heat flux (M K / S )
        REAL    ZSTK ( NLAYS )  !  zh( l )   - stkht   (m)
        REAL    TSTK            !  temperature at top of stack (deg K)
        REAL    WSTK            !  wind speed  at top of stack (m/s)
        INTEGER LSTK            !  first L: ZH(L) > STKHT
        INTEGER LPBL            !  first L: ZH(L) > HMIX
        REAL    TOP             !  plume top    elevation (m)
        REAL    BOT             !  plume bottom elevation (m)


C...........   PARAMETERS and their descriptions:

        REAL        GAMA        !  ?? plume spread parameter
        REAL        HCRIT	!  criterion for HSTAR
        REAL        CRDIST	!  criterion distance between HMIX and HS
        REAL        SMALL	!  Criterion for stable stability
        REAL        D3     	!  1/ 3
        REAL        D30		!  1/30
        REAL        D2664	!  1/ 2.664
        REAL        D59319	!  1/59.319
        REAL        TWOTHD	!  2/ 3
        REAL        FIVETHD	!  5/ 3 = 1 + 2/3

        PARAMETER ( GAMA   = -0.0098 ,
     &              CRDIST = 200.0 ,
     &              HCRIT  =   1.0E-4 * 0.03 ,	! hfx min * tolerance
     &              SMALL  =   1.0E-5 ,
     &              D3     = 1.0 /  3.0 ,
     &              D30    = 1.0 / 30.0 ,
     &              D2664  = 1.0 /  2.664 ,
     &              D59319 = 1.0 / 59.319 ,
     &              TWOTHD = 2.0 /  3.0 ,
     &             FIVETHD = 5.0 /  3.0 )


C...........   EXTERNAL FUNCTIONS and their descriptions:

        REAL        POLY
        EXTERNAL    POLY


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER M               !  layer counter
        INTEGER IQ              !  stab class:  1-unstbl,2-neut,3-stbl
        INTEGER LPLM            !  first L: ZH(L) > STKHT
        INTEGER NITER           !  iteration counter
        REAL    HSTAR           !  convective scale at stack (m?)
        REAL    TPLM            !  temperature at top of plume (m/s)
        REAL    WPLM            !  wind speed  at top of plume (m/s)
        REAL    BFLX            !  buoyancy flux (m**4/s**3)
        REAL    RBFLX           !  residual buoyancy flux (m**4/s**3)
        REAL    ZMIX            !  hmix - hs
        REAL    ZPLM            !  current plume height
        REAL    DH		!  plume rise increment
        REAL    DHN
        REAL    P, R, S         !  scratch coefficients


C...........   STATEMENT FUNCTIONS:

        REAL    B, H, U, US     !  arguments

        REAL    NEUTRL		!  neutral-stability plume rise function
        REAL    STABLE		!  stable            plume rise function
        REAL    UNSTBL		!  unstable          plume rise function

        NEUTRL( H, B, U, US ) =
     &     MIN( 10.0 * H, 
     &          1.2 * ( (              B/( U*US*US )   )**3 *
     &                  (  H + 1.3 * ( B/( U*US*US ) ) )**2 )**0.2  )

        STABLE( B, U, S ) =  2.6 * ( B / ( U * S ) )**D3

        UNSTBL( B, U )    = 30.0 * ( B / U )**0.6

C***********************************************************************
C   begin body of subroutine  PLMRIS

C.......   Compute buoyancy flux, convective scale.

        LPLM  = LSTK
        HSTAR = GRAV * HFX / TA( 1 )
        BFLX  = MAX( 0.0 , 
     &          0.25 * PI * ( STKTK - TSTK ) * STKVE * STKDM**2 / TSTK
     &             )

C.......   Compute initial plume rise from stack top to next level surface:

        IF( HSTAR .GT. HCRIT ) THEN		!  unstable case:
            
            ZMIX = HMIX - STKHT

            IF ( ZMIX .LE. 0.0 ) THEN           !  Stack above mixing height:

                LPLM = MIN( EMLAYS-1, LPBL+1 )
                S    = MAX( GRAV * DTHDZ( LPLM ) / TSTK, SMALL )
                
                IF ( LPLM .NE. LSTK ) THEN
                    M    = MIN( MAX( 1, LPLM - 2 ), EMLAYS - 3 )
                    WSTK = POLY( ZSTK( LPLM ), ZH( M ), WSPD( M ), 3 )
                END IF

                IF ( S .LE. SMALL ) THEN
                    DH = 1.5 * NEUTRL( STKHT, BFLX, WSTK, USTAR )
                    IQ = 2
                ELSE
                    DH = 1.5 * STABLE( BFLX, WSTK, S )
                    IQ = 3
                END IF

                TOP = STKHT + DH
                BOT = STKHT + D3 * DH

                RETURN

            ELSE IF ( ZMIX .LE. CRDIST ) THEN	!  need to compute penetration:

                S   = MAX( GRAV * DTHDZ( LSTK ) / TSTK, SMALL )
                DH  = 1.5 * STABLE( BFLX, WSTK, S )

                IF ( ZMIX .GE. DH ) THEN
                    TOP = HMIX
                    BOT = 0.0
                ELSE
                    TOP = ZF( LSTK )
                    BOT = 0.0
                END IF

                RETURN

            ELSE				!  unstable case:

                S   = MAX( GRAV * DTHDZ( LSTK ) / TSTK, SMALL )
                DH  = 1.5 * UNSTBL( BFLX, WSTK )
                DHN = 1.5 * NEUTRL( STKHT, BFLX, WSTK, USTAR )

                IF ( DHN .LT. DH ) THEN
                    DH = DHN
                    IQ = 2
                ELSE
                    IQ = 1
                END IF

            END IF

            TOP = STKHT +       DH
            BOT = STKHT + 0.5 * DH

        ELSE IF( HSTAR .LT. -HCRIT ) THEN      !  stable case:

            S   = MAX( GRAV * DTHDZ( LSTK ) / TSTK, SMALL )
            DH  = 1.5 * STABLE( BFLX, WSTK, S )
            DHN = 1.5 * NEUTRL( STKHT, RBFLX, WSTK, USTAR )

            IF ( DHN .LT. DH ) THEN
                DH = DHN
                IQ = 2
            ELSE
                IQ = 3
            END IF

        ELSE					!  neutral case:

            S   = MAX( GRAV * DTHDZ( LSTK ) / TSTK, SMALL )
            DH  = 1.5 * NEUTRL( STKHT, RBFLX, WSTK, USTAR )
            IQ  = 2
            
        END IF			!  hstar ==> unstable, stable, or neutral


C.......   Compute further plume rise from between level surfaces:

        NITER = 0
        RBFLX = BFLX
        ZPLM  = DH
        LPLM  = LSTK

101     CONTINUE	!  loop computing further plume rise

            IF( ZPLM .LE. ZSTK( LPLM ) ) THEN
                GO TO  199  ! exit plume rise loop
            ELSE IF ( LPLM .LT. EMLAYS ) THEN
                LPLM = LPLM + 1
            ELSE
                ZPLM = MIN( ZPLM, ZF( EMLAYS ) - STKHT )
                GO TO  199  ! exit plume rise loop
            END IF
            NITER = NITER + 1
            
C...........   compute residual bflx by stability case IQ:
            
            R = ZPLM
            IF ( R .LE. 0.0 ) THEN
                RBFLX = 0.0
            ELSE IF( IQ .EQ. 1 ) THEN	!  now compute resid bflx by stab case:
                R     = D30 * R
                RBFLX = WPLM * R**FIVETHD
            ELSE IF ( IQ .EQ. 2 ) THEN
                P = STKHT + TWOTHD * ZPLM
                RBFLX = D2664 * R * WPLM * USTAR**2 * (R / P )**TWOTHD
            ELSE	!  else iq = 3:
                RBFLX = D59319 * WPLM * S * R**3
            END IF	!  if stability flag iq is 1, 2, or 3
            
C...........   Re-interpolate met data:
            
            M    = MIN( MAX( 1, LPLM - 2 ), EMLAYS - 3 )
            WPLM = POLY( ZPLM,  ZH( M ), WSPD( M ), 3 )
            TPLM = POLY( ZPLM,  ZH( M ),   TA( M ), 3 )
            S    = GRAV * DTHDZ( LPLM ) / TPLM
            
C...........   Process according to stability cases:
            
            IF( S .GT. SMALL ) THEN               ! stable case:
                
                DH  = 1.5 * STABLE( RBFLX, WPLM, S )
                DHN = 1.5 * NEUTRL( STKHT, RBFLX, WPLM, USTAR )
                IF ( DHN .LT. DH ) THEN
                    DH = DHN
                    IQ = 2
                ELSE
                    IQ = 3
                END IF

            ELSE          ! if upper layer is not stable, use neutral formula

                S  = SMALL
                DH = 1.5 * NEUTRL( STKHT, RBFLX, WPLM, USTAR )
                IQ = 2
                
            END IF
            
            ZPLM = ZSTK( LPLM-1 ) + DH

            GO TO  101

199     CONTINUE	!  end loop computing further plume rise


C.......   Compute plume spread:

        TOP = STKHT +      ZPLM
        BOT = STKHT + D3 * ZPLM

        RETURN
        END

