
        SUBROUTINE  PLMRIS( ZH, ZF, WSPD, TA, DTHDZ,
     &                      USTAR, HMIX, HFX,
     &                      STKHT, STKDM, STKTK, STKVE, 
     &                      ZSTK, TSTK, WSTK, LSTK, LPBL, 
     &                      TOP, BOT )

C***********************************************************************
C  subroutine body starts at line  157
C
C  DESCRIPTION:  
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

        INCLUDE 'CONST3.EXT'    ! physical and mathematical constants
        INCLUDE 'PARMS3.EXT'    ! I/O API constants
        INCLUDE 'FDESC3.EXT'    ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
        INCLUDE 'VGRD3.EXT'     !  vertical direction parms


C...........   ARGUMENTS and their descriptions:

        REAL    ZH   ( EMLAYS ) !  layer center height (m) 
        REAL    ZF   ( 0:EMLAYS )!  layer surface height (m)
        REAL    WSPD ( EMLAYS ) !  wind speed (m/s)
        REAL    TA   ( EMLAYS ) !  temperature (deg K)
        REAL    DTHDZ( EMLAYS ) !  gradient of THETV
        REAL    USTAR           !  friction velocity (m/s)
        REAL    STKHT           !  stack height (m)
        REAL    STKDM           !  stack diameter (m)
        REAL    STKTK           !  exhaust temperature (deg K)
        REAL    STKVE           !  exhaust velocity (m/s)
        REAL    HMIX            !  mixing height (m)
        REAL    HFX             !  sensible heat flux (M K / S )
        REAL    ZSTK ( EMLAYS ) !  zf( l )   - stkht   (m)
        REAL    TSTK            !  temperature at top of stack (deg K)
        REAL    WSTK            !  wind speed  at top of stack (m/s)
        INTEGER LSTK            !  Layer of top of stack! same as RADM's KSTK
        INTEGER LPBL            !  Layer of height of PBL! same as RADM's KMIX
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

        INTEGER IQ              !  stab class:  1-unstbl,2-neut,3-stbl
        INTEGER LPLM      !  first L: ZH(L) > Plume height ! same as RADM's KPR

        REAL    BFLX            !  buoyancy flux (m**4/s**3)
        REAL    DH		!  plume rise increment to center of the plume
        REAL    DHM             !  plume rise from momentum
        REAL    DHN             !  plume rise for neutral case
        REAL    DHT             !  plume rise increment to the top of the plume
        REAL    HSTAR           !  convective scale at stack (m**2/s**3)
        REAL    P, R, S         !  scratch coefficients
        REAL    RBFLX           !  residual buoyancy flux (m**4/s**3)
        REAL    TPLM            !  temperature at top of plume (m/s)
        REAL    WPLM            !  wind speed  at top of plume (m/s)
        REAL    ZMIX            !  hmix - hs
        REAL    ZPLM            !  current plume height above stack 
                                !    (can be greater than the height of the 
                                !     top of the EMLAYS layer)

C...........   STATEMENT FUNCTIONS:

        REAL    B, H, U, US     !  arguments

        REAL    NEUTRL		!  neutral-stability plume rise function
        REAL    STABLE		!  stable            plume rise function
        REAL    UNSTBL		!  unstable          plume rise function

        NEUTRL( H, B, U, US ) =
     &     MIN( 10.0 * H, 
     &          1.2 * ( (              B/( U*US*US )   )**0.6 *    ! pwr 3 * 0.2
     &                  (  H + 1.3 * ( B/( U*US*US ) ) )**0.4   )) ! pwr 2 * 0.2

        STABLE( B, U, S ) =  2.6 * ( B / ( U * S ) )**D3

        UNSTBL( B, U )    = 30.0 * ( B / U )**0.6

C***********************************************************************
C   begin body of subroutine  PLMRIS

C.......   Compute buoyancy flux, convective scale.

        HSTAR = GRAV * HFX / TA( 1 )   ! Using surface temperature is correct
        BFLX  = 0.25*GRAV * ( STKTK-TSTK ) * STKVE * STKDM**2 / STKTK

C.......   Initialize layer of plume
        LPLM  = LSTK

C.......   Compute momentum rise
        DHM   = 3.0 * STKDM * STKVE / WSTK

C.......   When BFLX <= zero, use momentum rise only
C.......   NOTE: This part of algorithm added based on Models-3 plume rise

        IF( BFLX .LE. 0.0 ) THEN

            TOP = STKHT + 1.5 * DHM
            BOT = STKHT + 0.5 * DHM
            RETURN

        ENDIF

C.......   Compute initial plume rise from stack top to next level surface:

        IF( HSTAR .GT. HCRIT ) THEN		!  unstable case:

            ZMIX = HMIX - STKHT

            IF ( ZMIX .LE. 0.0 ) THEN           !  Stack above mixing height:

                LPLM = MIN( EMLAYS-1, LPBL+1 )
                S    = MAX( GRAV * DTHDZ( LPLM ) / TSTK, SMALL )

C.................  Reset the wind speed at stack to the wind speed at plume
C.................  when the layer of the plume is not equal to the layer of
C.................  the stack.  This is from Models-3, and we have asked
C.................  EPA 10/8/97 why this is done and haven't gotten a response.
                IF( LPLM .NE. LSTK ) THEN
                    WSTK = WSPD( LPLM )
                ENDIF

                DHN = NEUTRL( STKHT, BFLX, WSTK, USTAR )
                DH  = STABLE( BFLX, WSTK, S )

C.............  NOTE: The Models-3 version of plumer rise recalculates the
C.............        momentum plume rise here with the new WSTK.  We have
C.............        asked EPA on 10/8/97 if this is a bug but have not heard.
C.............        DHM = 3.0 * STKDM * STKVE / WSTK

                IF( DHN .LT. DH ) THEN  ! Take the minimum of neutral and stable
                    IQ = 2
                    DH = DHN
                ELSE 
                    IQ = 3 
                ENDIF

                IF( DHM .GT. DH ) THEN
                    IQ = 4
                    DH = DHM
                ENDIF

                DHT = 1.5 * DH

            ELSE IF ( ZMIX .LE. CRDIST ) THEN	!  need to compute penetration

                S   = MAX( GRAV * DTHDZ( LSTK ) / TSTK, SMALL )
                DHT  = 1.5 * STABLE( BFLX, WSTK, S )

C...............  Set the plume rise based on the partial penetration method
C...............    in Models-3.  There were several problems in the original
C...............    RADM algorithm which have been changed.  This method is
C...............    still not the same as Daewon's paper because that is only
C...............    appropriate for a Gaussian plume rise model.
                IF ( ZMIX .GE. DHT ) THEN
                    TOP = HMIX
                    BOT = TWOTHD * TOP
                    BOT = MAX( BOT, STKHT )
                    BOT = 0.5 * BOT

                ELSE
                    TOP = ZF( LSTK )
                    BOT = TWOTHD * TOP
                    BOT = MAX( BOT, STKHT )
                    BOT = 0.5 * BOT

                END IF

                RETURN

            ELSE				!  unstable case:

                DH  = UNSTBL( BFLX, WSTK )
                DHN = NEUTRL( STKHT, BFLX, WSTK, USTAR )

                IF ( DHN .LT. DH ) THEN
                    DH = DHN
                    IQ = 2
                ELSE
                    IQ = 1
                END IF

                IF( DHM .GT. DH ) THEN
                    DH = DHM
                    IQ = 4 
                ENDIF

                DHT = 1.5 * DH

            END IF

        ELSE IF( HSTAR .LT. -HCRIT ) THEN      !  stable case:

            S   = MAX( GRAV * DTHDZ( LSTK ) / TSTK, SMALL )
            DHT = 1.5 * STABLE( BFLX, WSTK, S )
            DHN = 1.5 * NEUTRL( STKHT, BFLX, WSTK, USTAR )

            IF ( DHN .LT. DHT ) THEN
                DHT = DHN
                IQ = 2
            ELSE
                IQ = 3
            END IF

        ELSE					!  neutral case:

            DHT = 1.5 * NEUTRL( STKHT, BFLX, WSTK, USTAR )
            IQ  = 2
            
        END IF			!  hstar ==> unstable, stable, or neutral
  
C.......   Compute further plume rise from between level surfaces:

        RBFLX = BFLX
        ZPLM  = DHT

C.......   End calculations if the momentum rise was used in the calculation

        IF( IQ .EQ. 4 ) GO TO 199  ! to point past iterative bouyancy loop

C.......   NOTE: LPLM has been initialized at line 161, and may have been
C                reset at line 185
101     CONTINUE	!  loop computing further plume rise

            R = ZPLM - ZSTK( LPLM )
            IF( R .LE. 0.0 ) THEN
                GO TO  199  ! exit plume rise loop
            ELSE IF ( LPLM .LT. EMLAYS ) THEN
                LPLM = LPLM + 1
            ELSE
                ZPLM = MIN( ZPLM, ZF( EMLAYS ) - STKHT )
                GO TO  199  ! exit plume rise loop
            END IF

C...........   Re-set met data. NOTE: the original RADM code submitted the 
C...........   WSPD and TA to a interpolator, but then requested the a height of
C...........   interpolation identical to ZH( LPLM ).

            WPLM = WSPD( LPLM )
            TPLM = TA  ( LPLM )

C...........   compute residual bflx by stability case IQ:

            IF( IQ .EQ. 1 ) THEN	!  now compute resid bflx by stab case:
                R     = D30 * R
                RBFLX = WPLM * R**FIVETHD
            ELSE IF ( IQ .EQ. 2 ) THEN
                P = STKHT + TWOTHD * ZPLM
                RBFLX = D2664 * R * WPLM * USTAR**2 * ( R / P )**TWOTHD
            ELSE	!  else iq = 3:
                RBFLX = D59319 * WPLM * S * R**3
            END IF	!  if stability flag iq is 1, 2, or 3

C...........   Process according to stability cases:

            S    = GRAV * DTHDZ( LPLM ) / TPLM
            IF( S .GT. SMALL ) THEN               ! stable case:

                DHT = 1.5 * STABLE( RBFLX, WPLM, S )
                DHN = 1.5 * NEUTRL( STKHT, RBFLX, WPLM, USTAR )
                IF ( DHN .LT. DHT ) THEN
                    DHT = DHN
                    IQ  = 2
                ELSE
                    IQ  = 3
                END IF

            ELSE          ! if upper layer is not stable, use neutral formula
                            
                DHT = 1.5 * NEUTRL( STKHT, RBFLX, WPLM, USTAR )
                IQ  = 2

            END IF

            ZPLM = ZSTK( LPLM-1 ) + DHT

            GO TO  101

199     CONTINUE	!  end loop computing further plume rise

C.......   Compute plume spread:

        TOP = STKHT +      ZPLM
        BOT = STKHT + D3 * ZPLM

        RETURN
        END

