
        SUBROUTINE GRPXREF( KEY2, KEY3, MONCODE, WEKCODE, DIUCODE, 
     &                      MPR, WPR, DPR, MAX1, MAX2, MAX3, 
     &                      N, CNT2, CNT3 )

C***********************************************************************
C  subroutine body starts at line  74
C
C  DESCRIPTION:  
C       Groups x-ref entries based on keys for use in processing various x-ref
C       inputs.  MAX3 and CNT3 should be 1 if MPR, WPR, and DPR are 2-d.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 4/96 by M Houyoux
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
        INTEGER      DIUCODE      ! Value of current DPR
        INTEGER      KEY2         ! Second key (0 is default)
        INTEGER      KEY3         ! Third key  (0 is default)
        INTEGER      MAX1         ! Dimension 1
        INTEGER      MAX2         ! Dimension 2
        INTEGER      MAX3         ! Dimension 3
        INTEGER      MONCODE      ! Value for current MPR
        INTEGER      N            ! Location in MPR, DPR, and WPR for dimension 1
        INTEGER      WEKCODE      ! Value for current WPR

        INTEGER      MPR( MAX1, MAX2, MAX3 )  ! Monthly profile codes
        INTEGER      WPR( MAX1, MAX2, MAX3 )  ! Weekly profile codes
        INTEGER      DPR( MAX1, MAX2, MAX3 )  ! Diurnal profile codes

C...........   LOCAL VARIABLES their descriptions:
       
        INTEGER      I, J         ! Loop counters

C***********************************************************************
C   begin body of subroutine GRPXREF

        IF ( KEY2 .EQ. 0 ) THEN
            IF ( KEY3 .EQ. 0 ) THEN

                DO 111 I = 1, CNT2
                    DO 110 J = 1, CNT3
                        MPR( N,I,J ) = MONCODE
                        WPR( N,I,J ) = WEKCODE
                        DPR( N,I,J ) = DIUCODE
  110               CONTINUE
  111           CONTINUE

            ELSE  ! KEY3 .NE. 0

                DO 113  I = 1, CNT2  
                    
                    MPR( N,I,KEY3 ) = MONCODE
                    WPR( N,I,KEY3 ) = WEKCODE
                    DPR( N,I,KEY3 ) = DIUCODE

                    DO 112  J = 1, CNT3
                        IF( MPR(N,I,J) .EQ. IMISS3 ) MPR(N,I,J)= MONCODE
                        IF( WPR(N,I,J) .EQ. IMISS3 ) WPR(N,I,J)= WEKCODE
                        IF( DPR(N,I,J) .EQ. IMISS3 ) DPR(N,I,J)= DIUCODE
  112               CONTINUE

  113           CONTINUE

            END IF

        ELSE  ! KEY2 .NE. 0
 
            IF ( KEY3 .EQ. 0 ) THEN

                DO  115  J = 1, CNT3
                    
                    MPR( N,KEY2,J ) = MONCODE
                    WPR( N,KEY2,J ) = WEKCODE
                    DPR( N,KEY2,J ) = DIUCODE

                    DO  114  I = 1, CNT2
                        IF( MPR(N,I,J) .EQ. IMISS3 ) MPR(N,I,J)= MONCODE
                        IF( WPR(N,I,J) .EQ. IMISS3 ) WPR(N,I,J)= WEKCODE
                        IF( DPR(N,I,J) .EQ. IMISS3 ) DPR(N,I,J)= DIUCODE
  114               CONTINUE

  115           CONTINUE

            ELSE  ! KEY3 .NE. 0
                
                MPR( N,KEY2,KEY3 ) = MONCODE
                WPR( N,KEY2,KEY3 ) = WEKCODE
                DPR( N,KEY2,KEY3 ) = DIUCODE

                DO 117  I = 1, CNT2
                    DO 116  J = 1, CNT3
                        IF( MPR(N,I,J) .EQ. IMISS3 ) MPR(N,I,J)= MONCODE
                        IF( WPR(N,I,J) .EQ. IMISS3 ) WPR(N,I,J)= WEKCODE
                        IF( DPR(N,I,J) .EQ. IMISS3 ) DPR(N,I,J)= DIUCODE
  116               CONTINUE
  117           CONTINUE
            END IF
        END IF
 
        RETURN

        END
