
        SUBROUTINE  PLMRIS( EMLAYS, LPBL, LSTK, HFX, HMIX, STKDM, STKHT, 
     &                      STKTK, STKVE, TSTK, DTHDZ, TA, WSPD, ZF, ZH,
     &                      ZSTK, WSTK, ZPLM )

C***********************************************************************
C  subroutine body starts at line 141
C
C  DESCRIPTION:  
C       computes plume rise height
C
C  PRECONDITIONS REQUIRED:
C       meteorology and stack parameters
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC, based on Briggs algorithm adapted from
C         RADM 2.6 subroutine PLUMER() (but with completely different 
C         data structuring).
C       Copied from plmris.F 4.4 by M Houyoux 3/99 
C       Updated with code from Jim Godowitch at EPA, 9/03
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
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'    ! I/O API constants
        INCLUDE 'FDESC3.EXT'    ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
        INCLUDE 'CONST3.EXT'    ! physical and mathematical constants

C...........   ARGUMENTS and their descriptions:

        INTEGER, INTENT (IN) :: EMLAYS    ! no. of emission layers
        INTEGER, INTENT (IN) :: LPBL      ! lyr of height of PBL, = RADM's KMIX
        INTEGER, INTENT (IN) :: LSTK      ! lyr of top of stack, = RADM's KSTK
        REAL   , INTENT (IN) :: HFX            ! sensible heat flux (M K / S )
        REAL   , INTENT (IN) :: HMIX           ! mixing height (m)
        REAL   , INTENT (IN) :: STKDM          ! stack diameter (m)
        REAL   , INTENT (IN) :: STKHT          ! stack height (m)
        REAL   , INTENT (IN) :: STKTK          ! exhaust temperature (deg K)
        REAL   , INTENT (IN) :: STKVE          ! exhaust velocity (m/s)
        REAL   , INTENT (IN) :: TSTK           ! tmptr at top of stack (deg K)
        REAL   , INTENT (IN) :: DTHDZ( EMLAYS )! gradient of THETV
        REAL   , INTENT (IN) :: TA   ( EMLAYS )! temperature (deg K)
        REAL   , INTENT (IN) :: WSPD ( EMLAYS )! wind speed (m/s)
        REAL   , INTENT (IN) :: ZF ( 0:EMLAYS )! layer surface height (m)
        REAL   , INTENT (IN) :: ZH   ( EMLAYS )! layer center height (m) 
        REAL   , INTENT (IN) :: ZSTK ( EMLAYS )! zf( l )   - stkht   (m)
        REAL, INTENT(IN OUT) :: WSTK           ! wind speed @ top of stack (m/s)
        REAL   , INTENT(OUT) :: ZPLM           ! current plume height above stack

C...........   PARAMETERS and their descriptions:

        REAL, PARAMETER :: GAMA    = -0.0098         ! ?? plume spread param
        REAL, PARAMETER :: HCRIT   =  1.0E-4 * 0.03  ! hfx min * tolerance
        REAL, PARAMETER :: CRDIST  =200.0            ! crit dstnce b/w HMIX & HS
        REAL, PARAMETER :: SMALL   =  1.0E-5         ! Criterion for stability
        REAL, PARAMETER :: D3      =  1.0 /  3.0     ! 1/ 3
        REAL, PARAMETER :: D30     =  1.0 / 30.0     ! 1/30
        REAL, PARAMETER :: D1355   =  1.0 /  1.355   ! 1/ 1.355
        REAL, PARAMETER :: D17576  =  1.0 / 17.576   ! 1/17.576
        REAL, PARAMETER :: TWOTHD  =  2.0 /  3.0     ! 2/ 3
        REAL, PARAMETER :: FIVETHD =  5.0 /  3.0     ! 5/ 3S
        REAL, PARAMETER :: NODIV0  =  1.0E-10        ! Prevent divide by zero

C...........   EXTERNAL FUNCTIONS and their descriptions:

        REAL        POLY
        EXTERNAL    POLY

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER IQ              !  stab class:  1-unstbl,2-neut,3-stbl
        INTEGER LPLM      !  first L: ZH(L) > Plume height ! same as RADM's KPR

        REAL    BFLX            !  buoyancy flux (m**4/s**3)
        REAL    DH              !  plume rise increment to center of the plume
        REAL    DHM             !  plume rise from momentum
        REAL    DHN             !  plume rise for neutral case
        REAL    DHT             !  plume rise increment to the top of the plume
        REAL    HSTAR           !  convective scale at stack (m**2/s**3)
        REAL    P, R, S         !  scratch coefficients
        REAL    RBFLX           !  residual buoyancy flux (m**4/s**3)
        REAL    TPLM            !  temperature at top of plume (m/s)
        REAL    WPLM            !  wind speed  at top of plume (m/s)
        REAL    ZMIX            !  hmix - hs

C...........   STATEMENT FUNCTIONS:

        REAL    B, H, U         !  arguments

        REAL    NEUTRL          !  neutral-stability plume rise function
        REAL    STABLE          !  stable            plume rise function
        REAL    UNSTBL          !  unstable          plume rise function

        NEUTRL( H, B, U ) =
     &     MIN( 10.0 * H, 
     &          1.2 * ( (    144.*B/( U*WSPD(1)*WSPD(1) )   )**0.6 *    ! pwr 3 * 0.2
     &          ( H + 1.3 * (144.*B/( U*WSPD(1)*WSPD(1) ) ) )**0.4   )) ! pwr 2 * 0.2

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
        
            ZPLM = STKHT + MAX( DHM, +2. )

            RETURN

        ENDIF

C.......   Compute initial plume rise from stack top to next level surface:

        IF( HSTAR .GT. HCRIT ) THEN             !  unstable case:

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
                    IF( WSTK .EQ. 0. ) WSTK = NODIV0
                ENDIF

                DHN = NEUTRL( STKHT, BFLX, WSTK )
                DH  = STABLE( BFLX, WSTK, S )

C.............  NOTE- The Models-3 version of plume rise recalculates the
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

                DHT = DH

            ELSE                        !  unstable case:

                DH  = UNSTBL( BFLX, WSTK )
                DHN = NEUTRL( STKHT, BFLX, WSTK )

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

                DHT = DH
     
            END IF

        ELSE IF( HSTAR .LT. -HCRIT ) THEN      !  stable case:

            S   = MAX( GRAV * DTHDZ( LSTK ) / TSTK, SMALL )
            DHT =  STABLE( BFLX, WSTK, S )
            DHN =  NEUTRL( STKHT, BFLX, WSTK )

            IF ( DHN .LT. DHT ) THEN
                DHT = DHN
                IQ = 2
            ELSE
                IQ = 3
            END IF

        ELSE                                   !  neutral case:

            DHT =  NEUTRL( STKHT, BFLX, WSTK )
            IQ  = 2
            
        END IF                 !  hstar ==> unstable, stable, or neutral
  
C.......   Compute further plume rise from between level surfaces:

        RBFLX = BFLX
        ZPLM  = DHT

C.......   End calculations if the momentum rise was used in the calculation

        IF( IQ .EQ. 4 ) GO TO 199  ! to point past iterative bouyancy loop

C.......   NOTE- LPLM has been initialized at line 145, and may have been
C                reset at line 169
        DO       !  loop computing further plume rise

            R = ZPLM - ZSTK( LPLM )
            IF( R .LE. 0.0 ) THEN
                EXIT  ! exit plume rise loop

            ELSE IF ( LPLM .LT. EMLAYS ) THEN
                LPLM = LPLM + 1

            ELSE
                ZPLM = MIN( ZPLM, ZF( EMLAYS ) - STKHT )
                EXIT  ! exit plume rise loop

            END IF

C...........   Re-set met data. NOTE- the original RADM code submitted the 
C...........   WSPD and TA to a interpolator, but then requested the a height of
C...........   interpolation identical to ZH( LPLM ).

            WPLM = WSPD( LPLM )
            TPLM = TA  ( LPLM )

C...........   Compute residual bflx by stability case IQ:

            IF( IQ .EQ. 1 ) THEN           !  now compute resid bflx by stab case:
                R     = D30 * R
                RBFLX = WPLM * R**FIVETHD
            ELSE IF ( IQ .EQ. 2 ) THEN
                P = STKHT + TWOTHD * ZPLM
                RBFLX = D1355 * R * WPLM * ( ( WSPD( 1 )**2. )/144. ) 
     &               * ( R / P )**TWOTHD
            ELSE        !  else iq = 3:
                RBFLX = D17576 * WPLM * S * R**3
            END IF      !  if stability flag iq is 1, 2, or 3

C...........   Prevent divide-by-zero by WPLM

            IF( WPLM .EQ. 0. ) WPLM = NODIV0

C...........   Process according to stability cases:

            S    = GRAV * DTHDZ( LPLM ) / TPLM
            IF( S .GT. SMALL ) THEN               ! stable case:

               DHT =  STABLE( RBFLX, WPLM, S )
               DHN =  NEUTRL( STKHT, RBFLX, WPLM )
               IF ( DHN .LT. DHT ) THEN
                    DHT = DHN
                    IQ  = 2
                ELSE
                    IQ  = 3
                END IF

            ELSE          ! if upper layer is not stable, use neutral formula

                DHN =  NEUTRL( STKHT, RBFLX, WPLM )
                DH  = UNSTBL( BFLX, WSTK )
                IQ = 1
                IF ( DHN .LT. DH ) THEN
                    DH = DHN
                    IQ  = 2
                ENDIF
                DHT = DH

            END IF
 
            ZPLM = ZSTK( LPLM-1 ) + DHT

        END DO    !  end loop computing further plume rise

199     CONTINUE

C.........  Adjust for layer 1 combustion pt. source stacks with plume
C           rise limited to layer 1; put plume height in middle of layer 2  
        IF ( ZPLM .LE. ZF( 1 ) .AND. STKTK .GT. TA( 1 ) ) THEN
            ZPLM = ZH( 2 )
        ENDIF

C.........   Determine actual height of plume centerline after rise
        ZPLM = ZPLM + STKHT
      
        RETURN

        END SUBROUTINE PLMRIS

