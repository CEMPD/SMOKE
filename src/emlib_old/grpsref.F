
        SUBROUTINE GRPSREF( KEY2, KEY3, PROFID, FAC1, FAC2, FAC3, FAC4,
     &                      PROF, ARR1, ARR2, ARR3, ARR4, 
     &                      MAX1, MAX2, MAX3, N, CNT2, CNT3 )

C***********************************************************************
C  subroutine body starts at line  79
C
C  DESCRIPTION:  
C       Groups x-ref entries based on keys for use in processing various x-ref
C        inputs.  MAX3 and CNT3 should be 1 if PROF, ARR1, and ARR2 are 2-d.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 10/96 by M Houyoux  -  copied from grpxref.F
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
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'

C...........   ARGUMENTS and their descriptions:

        INTEGER      CNT2         ! Loop size for dimension 2
        INTEGER      CNT3         ! Loop size for dimension 3
        INTEGER      KEY2         ! Second key (0 is default)
        INTEGER      KEY3         ! Third key  (0 is default)
        INTEGER      MAX1         ! Dimension 1
        INTEGER      MAX2         ! Dimension 2
        INTEGER      MAX3         ! Dimension 3
        INTEGER      PROFID       ! Value for current PROF
        INTEGER      N            ! Location in PROF, ARR2, and ARR1 
                                  !   for dimension 1
        REAL         FAC1         ! Value of current ARR1
        REAL         FAC2         ! Value of current ARR2
        REAL         FAC3         ! Value of current ARR3
        REAL         FAC4         ! Value of current ARR4

        INTEGER      PROF( MAX1, MAX2, MAX3 )  ! Profile codes
        REAL         ARR1( MAX1, MAX2, MAX3 )  ! 1st array
        REAL         ARR2( MAX1, MAX2, MAX3 )  ! 2nd array
        REAL         ARR3( MAX1, MAX2, MAX3 )  ! 3rd array
        REAL         ARR4( MAX1, MAX2, MAX3 )  ! 4th array

C...........   LOCAL VARIABLES their descriptions:
       
        INTEGER      I, J         ! Loop counters

C***********************************************************************
C   begin body of subroutine GRPSREF

        IF ( KEY2 .EQ. 0 ) THEN
            IF ( KEY3 .EQ. 0 ) THEN

                DO 111 I = 1, CNT2
                    DO 110 J = 1, CNT3
                        PROF( N,I,J ) = PROFID
                        ARR1( N,I,J ) = FAC1
                        ARR2( N,I,J ) = FAC2
                        ARR3( N,I,J ) = FAC3
                        ARR4( N,I,J ) = FAC4
  110               CONTINUE
  111           CONTINUE

            ELSE  ! KEY3 .NE. 0

                DO 113  I = 1, CNT2  
                    
                    PROF( N,I,KEY3 ) = PROFID
                    ARR1( N,I,KEY3 ) = FAC1
                    ARR2( N,I,KEY3 ) = FAC2
                    ARR3( N,I,KEY3 ) = FAC3
                    ARR4( N,I,KEY3 ) = FAC4

                    DO 112  J = 1, CNT3
                     IF( PROF(N,I,J) .LE. IMISS3 ) PROF(N,I,J)= PROFID
                     IF( ARR1(N,I,J) .LE. AMISS3 ) ARR1(N,I,J)= FAC1
                     IF( ARR2(N,I,J) .LE. AMISS3 ) ARR2(N,I,J)= FAC2
                     IF( ARR3(N,I,J) .LE. AMISS3 ) ARR3(N,I,J)= FAC3
                     IF( ARR4(N,I,J) .LE. AMISS3 ) ARR4(N,I,J)= FAC4
  112               CONTINUE

  113           CONTINUE

            END IF

        ELSE  ! KEY2 .NE. 0
 
            IF ( KEY3 .EQ. 0 ) THEN

                DO  115  J = 1, CNT3
                    
                    PROF( N,KEY2,J ) = PROFID
                    ARR1( N,KEY2,J ) = FAC1
                    ARR2( N,KEY2,J ) = FAC2
                    ARR3( N,KEY2,J ) = FAC3
                    ARR4( N,KEY2,J ) = FAC4

                    DO  114  I = 1, CNT2
                     IF( PROF(N,I,J) .LE. IMISS3 ) PROF(N,I,J)= PROFID
                     IF( ARR1(N,I,J) .LE. AMISS3 ) ARR1(N,I,J)= FAC1
                     IF( ARR2(N,I,J) .LE. AMISS3 ) ARR2(N,I,J)= FAC2
                     IF( ARR3(N,I,J) .LE. AMISS3 ) ARR3(N,I,J)= FAC3
                     IF( ARR4(N,I,J) .LE. AMISS3 ) ARR4(N,I,J)= FAC4
  114               CONTINUE

  115           CONTINUE

            ELSE  ! KEY3 .NE. 0
                
                PROF( N,KEY2,KEY3 ) = PROFID
                ARR1( N,KEY2,KEY3 ) = FAC1
                ARR2( N,KEY2,KEY3 ) = FAC2
                ARR3( N,KEY2,KEY3 ) = FAC3
                ARR4( N,KEY2,KEY3 ) = FAC4

                DO 117  I = 1, CNT2
                    DO 116  J = 1, CNT3
                     IF( PROF(N,I,J) .LE. IMISS3 ) PROF(N,I,J)= PROFID
                     IF( ARR1(N,I,J) .LE. AMISS3 ) ARR1(N,I,J)= FAC1
                     IF( ARR2(N,I,J) .LE. AMISS3 ) ARR2(N,I,J)= FAC2
                     IF( ARR3(N,I,J) .LE. AMISS3 ) ARR3(N,I,J)= FAC3
                     IF( ARR4(N,I,J) .LE. AMISS3 ) ARR4(N,I,J)= FAC4
  116               CONTINUE
  117           CONTINUE
            END IF
        END IF
 
        RETURN

        END
