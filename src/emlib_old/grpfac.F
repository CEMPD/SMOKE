
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE GRPFAC( KEY2, KEY3, FAC,
     &                      FTR, MAX1, MAX2, MAX3, 
     &                      N, CNT2, CNT3 )

C***********************************************************************
C  subroutine body starts at line  48
C
C  FUNCTION:  Groups factor entries based on keys for use in processing
C             various profile inputs.  MAX3 and CNT3 should be 1 if
C             FAC is 2-d.
C
C  PRECONDITIONS REQUIRED: 
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 4/96 by M Houyoux
C
C***********************************************************************

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
        INTEGER      N            ! Location in FTR for dimension 1

        REAL         FAC          ! Current factors
        REAL         FTR( MAX1, MAX2, MAX3 ) ! List of all factors
       
C...........   LOCAL VARIABLES their descriptions:

        INTEGER      I, J         ! Loop counters

C***********************************************************************
C   begin body of subroutine GRPFAC

        IF ( KEY2 .EQ. 0 ) THEN
            IF ( KEY3 .EQ. 0 ) THEN

                DO 111 I = 1, CNT2
                    DO 110 J = 1, CNT3
                        FTR( N,I,J ) = FAC
  110               CONTINUE
  111           CONTINUE

            ELSE  ! KEY3 .NE. 0

                DO 113  I = 1, CNT2  
                    
                    FTR( N,I,KEY3 ) = FAC

                    DO 112  J = 1, CNT3
                        IF( FTR(N,I,J) .EQ. AMISS3 ) FTR(N,I,J)= FAC
  112               CONTINUE

  113           CONTINUE

            END IF

        ELSE  ! KEY2 .NE. 0
 
            IF ( KEY3 .EQ. 0 ) THEN

                DO  115  J = 1, CNT3
                    
                    FTR( N,KEY2,J ) = FAC

                    DO  114  I = 1, CNT2
                        IF( FTR(N,I,J) .EQ. AMISS3 ) FTR(N,I,J)= FAC
  114               CONTINUE

  115           CONTINUE

            ELSE  ! KEY3 .NE. 0
                
                FTR( N,KEY2,KEY3 ) = FAC

                DO 117  I = 1, CNT2
                    DO 116  J = 1, CNT3
                        IF( FTR(N,I,J) .EQ. AMISS3 ) FTR(N,I,J)= FAC
  116               CONTINUE
  117           CONTINUE
            END IF
        END IF
 
        RETURN

        END
