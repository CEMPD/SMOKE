
        SUBROUTINE FIRE_PLMRIS( NS, EMLAYS, LPBL, LSTK, HFX, HMIX, BFLX,
     &                          TSTK, USTAR, DTHDZ, TA, WSPD, ZF, ZH,
     &                          ZSTK, WSTK, TOP, BOT, ZPLM )

C***********************************************************************
C  subroutine body starts at line 141
C
C  DESCRIPTION:  
C       computes plume top and bottom for fires
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
C       Updated with code from J. Godowitch and G. Pouliot, 2/05
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

        INTEGER, INTENT (IN) :: NS
        INTEGER, INTENT (IN) :: EMLAYS    ! no. of emission layers
        INTEGER, INTENT (IN) :: LPBL      ! lyr of height of PBL, = RADM's KMIX
        INTEGER, INTENT (IN) :: LSTK      ! lyr of top of stack, = RADM's KSTK
        REAL   , INTENT (IN) :: HFX            ! sensible heat flux (M K / S )
        REAL   , INTENT (IN) :: HMIX           ! mixing height (m)
        REAL   , INTENT (IN) :: BFLX           ! bouyancy flux (m^4/s^3)
        REAL   , INTENT (IN) :: TSTK           ! tmptr at top of stack (deg K)
        REAL   , INTENT (IN) :: USTAR          ! friction velocity (m/s)
        REAL   , INTENT (IN) :: DTHDZ( EMLAYS )! gradient of THETV
        REAL   , INTENT (IN) :: TA   ( EMLAYS )! temperature (deg K)
        REAL   , INTENT (IN) :: WSPD ( EMLAYS )! wind speed (m/s)
        REAL   , INTENT (IN) :: ZF ( 0:EMLAYS )! layer surface height (m)
        REAL   , INTENT (IN) :: ZH   ( EMLAYS )! layer center height (m) 
        REAL   , INTENT (IN) :: ZSTK ( EMLAYS )! zf( l )   - stkht   (m)
        REAL, INTENT(IN OUT) :: WSTK           ! wind speed @ top of stack (m/s)
        REAL   , INTENT(OUT) :: TOP            ! plume top elevation (m)
        REAL   , INTENT(OUT) :: BOT            ! plume bottom elevation (m)
        REAL   , INTENT(OUT) :: ZPLM           ! plume centerline height
                                               ! (can be greater than the height of
                                               !  the top of the EMLAYS layer)

C...........   PARAMETERS and their descriptions:

        REAL, PARAMETER :: GAMA    = -0.0098         ! ?? plume spread param
        REAL, PARAMETER :: HCRIT   =  1.0E-4 * 0.03  ! hfx min * tolerance
        REAL, PARAMETER :: CRDIST  =200.0            ! crit dstnce b/w HMIX & HS
        REAL, PARAMETER :: SMALL   =  3.0E-5         ! Criterion for stability
        REAL, PARAMETER :: D3      =  1.0 /  3.0
        REAL, PARAMETER :: D6      =  1.0 /  6.0
        REAL, PARAMETER :: D30     =  1.0 / 30.0
        REAL, PARAMETER :: D45     =  1.0 / 45.0
        REAL, PARAMETER :: D2664   =  1.0 /  2.664
        REAL, PARAMETER :: D59319  =  1.0 / 59.319
        REAL, PARAMETER :: TWOTHD  =  2.0 /  3.0
        REAL, PARAMETER :: FIVETHD =  5.0 /  3.0
        REAL, PARAMETER :: NODIV0  =  1.0            ! Prevent divide by zero
        REAL, PARAMETER :: ZERO    =  0.0

C...........   EXTERNAL FUNCTIONS and their descriptions:

        REAL        POLY
        EXTERNAL    POLY

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER IQ              !  stab class:  1-unstbl,2-neut,3-stbl
        INTEGER LPLM            !  first L: ZH(L) > Plume height ! same as RADM's KPR
        INTEGER NN              !  layer counter

        REAL    DH              !  plume rise increment to center of the plume
        REAL    DHM             !  plume rise from momentum
        REAL    DHSM            !  stable momentum plume rise
        REAL    DHN             !  plume rise for neutral case
        REAL    DHT             !  plume rise increment to the top of the plume
        REAL    HSTAR           !  convective scale at stack (m**2/s**3)
        REAL    P, R, S         !  scratch coefficients
        REAL    RBFLX           !  residual buoyancy flux (m**4/s**3)
        REAL    TPLM            !  temperature at top of plume (m/s)
        REAL    WPLM            !  wind speed  at top of plume (m/s)
        REAL    ZMIX            !  hmix - hs
        REAL    ZB              !  height of bottom of layer

C...........   STATEMENT FUNCTIONS:

        REAL    B, H, U, US     !  arguments

        REAL    NEUTRL          !  neutral-stability plume rise function
        REAL    STABLE          !  stable            plume rise function
        REAL    UNSTBL          !  unstable          plume rise function

        NEUTRL( H, B, U, US ) =
     &     1.2 * ((B/(U*US*US))**0.6 * (H+1.3*(B/(U*US*US)))**0.4)

        STABLE( B, U, S ) =  2.6 * ( B / ( U * S ) )**D3

        UNSTBL( B, U )    = 30.0 * ( B / U )**0.6

C***********************************************************************
C   begin body of subroutine  PLMRIS

C.........  Compute convective scale
        HSTAR = GRAV * HFX / TA( 1 )   ! Using surface temperature is correct

C.........  Initialize layer of plume
        LPLM  = LSTK

C.........  Set minimum wind speed to 1 m/s
        WSTK  = MAX( WSTK, 1.0 )

C.........  Set momentum rise to zero
        DHM   = ZERO

C.........  Compute initial plume rise from stack top to next level surface:
        IF( HSTAR > HCRIT ) THEN             !  unstable case:

            ZMIX = HMIX

            IF ( ZMIX <= 0.0 ) THEN           !  Stack above mixing height:

                S = MAX( GRAV * DTHDZ( LPLM ) / TSTK, SMALL )

C.................  Reset the wind speed at stack to the wind speed at plume
C.................  when the layer of the plume is not equal to the layer of
C.................  the stack.  This is from Models-3, and we have asked
C.................  EPA 10/8/97 why this is done and haven't gotten a response.
                IF( LPLM .NE. LSTK ) THEN
                    WSTK = WSPD( LPLM )
                    IF( WSTK == 0. ) WSTK = NODIV0
                ENDIF

C.................  Set stable momentum rise
                IF( DTHDZ( LPLM ) > 0.001 ) THEN
                    DHSM = ZERO
                ELSE
                    DHSM = DHM
                END IF
                
                DHM = MIN( DHSM, DHM )

C.................  Compute neutral and stable plume rises
                DHN = NEUTRL( ZERO, BFLX, WSTK, USTAR )
                DH  = STABLE( BFLX, WSTK, S )

                IF( DHN < DH ) THEN  ! Take the minimum of neutral and stable
                    IQ = 2
                    DH = DHN
                ELSE 
                    IQ = 3 
                ENDIF

                IF( DHM > DH .AND. WSTK > 1. ) THEN
                    IQ = 4
                    DH = DHM
                ENDIF

                DHT = 1.5 * DH

            ELSE                        !  unstable case:

                DH  = UNSTBL( BFLX, WSTK )
                DHN = NEUTRL( ZERO, BFLX, WSTK, USTAR )

                IF ( DHN < DH ) THEN
                    DH = DHN
                    IQ = 2
                ELSE
                    IQ = 1
                END IF

                IF( DHM > DH .AND. WSTK > 1. ) THEN
                    DH = DHM
                    IQ = 4 
                ENDIF

                DHT = 1.5 * DH
     
            END IF

        ELSE IF( HSTAR < -HCRIT .OR. DTHDZ( LSTK ) > 0.001 ) THEN      !  stable case:

            S   = MAX( GRAV * DTHDZ( LSTK ) / TSTK, SMALL )
            
            DHT = 1.5 * STABLE( BFLX, WSTK, S )
            DHN = 1.5 * NEUTRL( ZERO, BFLX, WSTK, USTAR )

            IF ( DHN < DHT ) THEN
                DHT = DHN
                IQ = 2
            ELSE
                IQ = 3
            END IF

        ELSE                                   !  neutral case:

            DHT = 1.5 * NEUTRL( ZERO, BFLX, WSTK, USTAR )
            IQ  = 2
            
        END IF                 !  hstar ==> unstable, stable, or neutral
  
C.........  Compute further plume rise from between level surfaces:
        NN = 0
        RBFLX = BFLX
        ZPLM  = DHT

C.........  End calculations if the momentum rise was used in the calculation
        IF( IQ == 4 ) GO TO 199  ! to point past iterative bouyancy loop


        DO       !  loop computing further plume rise

            R = ZPLM - ZSTK( LPLM )
            IF( R <= 0.0 ) THEN
                EXIT  ! exit plume rise loop
            END IF
            
            IF( LPLM == EMLAYS ) THEN
                ZPLM = MIN( ZPLM, ZSTK( EMLAYS ) )
                EXIT  ! exit plume rise loop
            END IF

C.............  Re-set met data. NOTE- the original RADM code submitted the 
C.............  WSPD and TA to a interpolator, but then requested the a height of
C.............  interpolation identical to ZH( LPLM ).
            NN = NN + 1
            IF( NN > 1 ) THEN
                WPLM = WSPD( LPLM )
                TPLM = TA  ( LPLM )
            ELSE     ! 1st time, use stack values
                WPLM = WSTK
                TPLM = TSTK
            END IF

C.............  Compute residual bflx by stability case IQ:

            IF( IQ == 1 ) THEN
                R     = D45 * R     ! includes 1.5 factor for plume top
                RBFLX = WPLM * ( R**FIVETHD )
                
            ELSE IF ( IQ == 2 ) THEN
                P     = TWOTHD * ZPLM
                RBFLX = D2664 * ( R**FIVETHD ) * WPLM * ( USTAR**2. ) / 
     &                  P**TWOTHD
     
            ELSE        !  else iq = 3:
                RBFLX = D59319 * WPLM * S * R**3
     
            END IF      !  if stability flag iq is 1, 2, or 3

            IF( LPLM < EMLAYS ) LPLM = LPLM + 1
            WPLM = WSPD( LPLM )
            TPLM = TA  ( LPLM )

C.............  Prevent divide-by-zero by WPLM
            IF( WPLM == 0. ) WPLM = NODIV0

C.............  Process according to stability cases:
            S    = GRAV * DTHDZ( LPLM ) / TPLM
            
            IF( S > SMALL ) THEN               ! stable case:

                DHT = 1.5 * STABLE( RBFLX, WPLM, S )
                DHN = 1.5 * NEUTRL( ZERO, RBFLX, WPLM, USTAR )
                IF ( DHN < DHT ) THEN
                    DHT = DHN
                    IQ  = 2
                ELSE
                    IQ  = 3
                END IF
                DH = DHT / 1.5

            ELSE          ! if upper layer is not stable, use neutral formula

                DHN = NEUTRL( ZERO, RBFLX, WPLM, USTAR )
                DH  = UNSTBL( RBFLX, WPLM )
                IQ = 1
                IF ( DHN < DH ) THEN
                    DH = DHN
                    IQ  = 2
                END IF
                DHT = 1.5 * DH

            END IF
 
            ZPLM = ZSTK( LPLM-1 ) + DHT
            DH   = ZSTK( LPLM-1 ) + DH

        END DO    !  end loop computing further plume rise

199     CONTINUE

C.........  Compute plume top and bottom and determine actual height of
C           plume centerline after rise
        TOP  = ZPLM
        DH   = TWOTHD * ZPLM
        BOT  = 0.5 * DH
        ZPLM = DH
      
        RETURN

        END SUBROUTINE FIRE_PLMRIS

