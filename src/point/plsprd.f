
       SUBROUTINE PLSPRD( DTHDZ, ZF, KZ, CEFSTK, HTMIX,
     &                    PLTOP, PLBOT )
     
C***********************************************************************
C  subroutine body starts at line 70
C
C  DESCRIPTION:  
C       Calculates the initial vertical spread of a plume; modified
C       from Gillani's model.
C
C  PRECONDITIONS REQUIRED:
C
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C
C  REVISION  HISTORY:
C       Created from code by Jim Godowitch at EPA, 9/03
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
       
C...........   ARGUMENTS and their descriptions:

       REAL,    INTENT  (IN) :: DTHDZ( KZ )    ! potential temperature lapse rate (K/m)
       REAL,    INTENT  (IN) :: ZF( 0:KZ )     ! full-layer heights (m)
       INTEGER, INTENT  (IN) :: KZ             ! number of emissions layers
       REAL,    INTENT  (IN) :: CEFSTK         ! effective stack height (m)
       REAL,    INTENT  (IN) :: HTMIX          ! mixing height (m)
       REAL,    INTENT (OUT) :: PLTOP          ! plume top (m)
       REAL,    INTENT (OUT) :: PLBOT          ! plume bottom (m)
       
C...........   PARAMETERS and their descriptions:
       REAL, PARAMETER :: SZ0FAC = 3.545    ! factor used to derive plume depth
       REAL, PARAMETER :: SPRFAC = 15.      ! empirical coefficient for vertical spread
       REAL, PARAMETER :: GAMA   = -0.0098  ! adiabatic lapse rate (K/m)
       
C...........   Local variables
       INTEGER    K
       
       REAL       SIGZ0
       REAL       DTDZ
       REAL       DPTH
       
C***********************************************************************
C   begin body of subroutine  PLSPRD

C........  Get ambient temperature above plume rise height (effective stack height)
       K = 0
       DO 
           K = K + 1
           IF( K == KZ .OR. CEFSTK <= ZF( K ) ) EXIT
       END DO
       DTDZ  = DTHDZ( K ) + GAMA
       
C........  Compute initial vertical spread
       SIGZ0 = MAX( 10.0, SPRFAC * EXP( -117. * DTDZ ) )
       DPTH  = SZ0FAC * SIGZ0
       
C........  Compute plume top and bottom heights; plume is either completely
C          within or outside mixing layer
       PLTOP = CEFSTK + DPTH/2.
       PLBOT = CEFSTK - DPTH/2.
       
C........  Make sure plume bottom is at least zero
       PLBOT = MAX( 0.0, PLBOT )

C........  Make sure that plume top and bottom heights are less than 
C          the top layer's top and bottom heights
       PLTOP = MIN( ZF( KZ ), PLTOP )
       PLBOT = MIN( ZF( KZ ) - 1., PLBOT )
       
       RETURN
       
       END SUBROUTINE PLSPRD
       
