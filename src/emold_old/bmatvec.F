
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  BMATVEC( M, N, P, IX, AX, V, C )

C***********************************************************************
C  subroutine body starts at line  43
C
C  FUNCTION:  apply a 4-band sparse matrix to an array ("layered vector")
C	      as used for bilinear interpolation of meteorology in LAYPOINT
C
C  NOTE:  Transposes from I/O subscript order V( M,P )
C         to computational subscript order    C( P,N )
C
C       For bilinear interpolation of gridded data having dimension NC,NR
C       to a list of locations having grid-normal coordinates <X(S),Y(S)>:
C       let C(S) and R(S) be INT( X(S) ) and INT( Y(S) ), P(S) and Q(S)
C       be AMOD( X(S), 1.0 ) and AMOD( Y(S), 1.0 ).  Then IX has the 
C	single-indexing subscripts into the grid for the cell=corners
C       surrounding the <X(S),Y(S)>
C           IX(1,S) = C + NC * ( R - 1 ) + 1
C           IX(2,S) = C + NC * ( R - 1 ) 
C           IX(3,S) = C + NC *   R
C           IX(4,S) = C + NC *   R       + 1
C       and AX has the bilinear-interpolation coefficients:
C           AX(1,S) = ( 1.0 - P( S ) )*( 1.0 - Q( S ) )
C           AX(2,S) =         P( S )  *( 1.0 - Q( S ) )
C           AX(3,S) = ( 1.0 - P( S ) )*        Q( S )       
C           AX(4,S) =         P( S )  *        Q( S )
C
C  PRECONDITIONS REQUIRED:
C       Index and coefficients set up so that the equation for
C       the multiplication (at each row R and layer L) is
C
C       C(L,R) = SUM_{ J=1...4 }  A( J,R ) * V( I( J,R ),L )
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 12/95 by CJC
C
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:
        
        INTEGER         M               ! length of input  vector
        INTEGER         N               ! length of output vector
        INTEGER         P               ! length of output vector
        INTEGER         IX( 4,N )       ! index array
        REAL            AX( 4,N )       ! 4-band coeff matrix
        REAL            V( M,P )        ! P-layered input  vector
        REAL            C( P,N )        ! P-layered output vector


C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        
        INTEGER         S, L, J1, J2, J3, J4


C***********************************************************************
C   begin body of subroutine  BMATVEC

        DO  22  S = 1, N
            
            J1 = IX( 1,S )
            J2 = IX( 2,S )
            J3 = IX( 3,S )
            J4 = IX( 4,S )

            DO  11  L = 1, P
                C( L,S ) = AX( 1,S ) * V( J1,L )  +
     &                     AX( 2,S ) * V( J2,L )  +
     &                     AX( 3,S ) * V( J3,L )  +
     &                     AX( 4,S ) * V( J4,L )
11          CONTINUE
        
22      CONTINUE

        RETURN
        END

