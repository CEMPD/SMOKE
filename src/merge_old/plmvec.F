
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  PLMVEC( FL, P, NP, IP, V )

C***********************************************************************
C  subroutine body starts at line  49
C
C  FUNCTION:  multiply vector P by sparse incidence matrix <NP,IP> 
C       and layer fractions FL, and add to the result V for CSGLMERGE 
C       program.
C
C  PRECONDITIONS REQUIRED:  none
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 2/96 by CJC
C
C***********************************************************************

        IMPLICIT NONE
        
        INCLUDE 'GRDIMS3.EXT'
        INCLUDE 'LADIMS3.EXT'
        INCLUDE 'PTDIMS3.EXT'
        INCLUDE 'ARDIMS3.EXT'

C...........   ARGUMENTS and their descriptions:
        
        REAL            P ( NPSRC )     ! point source input vector
        INTEGER         NP( NGRID )     ! # of entries per row in point matrix
        INTEGER         IP( NPSRC )     ! columns list         in point matrix
        REAL            FL( EMLAYS, NPSRC )     ! pt-srce layering fractions
        REAL            V ( NGRID, EMLAYS )     ! output vector


C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER         R, C, L, I, KP
        REAL            ACC( EMLAYS )


C***********************************************************************
C   begin body of subroutine  SMRGVEC

C.......   Layer 1 calculations use both area and point sources

        KP = 0
        DO  66  R = 1, NGRID

            DO  22  L = 1, EMLAYS
                ACC( L ) = V( R,L )
22          CONTINUE            

            DO  44  C = 1, NP( R )              !  loop on point sources C

                KP  = KP + 1
                I   = IP( KP )

                DO  33 L = 1, EMLAYS
                    ACC( L ) = ACC( L )  +  FL( L,R ) * P( I )
33              CONTINUE            

44          CONTINUE                            !  end loop on point sources C

            DO  55  L = 1, EMLAYS
                V( R,L ) = ACC( L )
55          CONTINUE            

66      CONTINUE        !  end loop on R (grid cells)

        RETURN
        END

