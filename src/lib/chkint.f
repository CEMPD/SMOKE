
        LOGICAL FUNCTION CHKINT( STRING )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function returns true if the string appears to be an integer
C      and false otherwise
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

        IMPLICIT NONE


C...........   SUBROUTINE ARGUMENTS
        CHARACTER*(*)   STRING   ! input character string

C...........   Other local variables
        INTEGER         K, L  ! counters and indices

        LOGICAL         SPACFLAG  ! true if already encountered a space in string
        LOGICAL         NEGVFLAG  ! true if already encountered a '-' in string

        CHARACTER*1     CBUF 
        CHARACTER*256   BUFFER 

        CHARACTER*16 :: PROGNAME = 'CHKINT' ! program name

C***********************************************************************
C   begin body of function CHKINT

        CHKINT = .TRUE.

        BUFFER = ADJUSTL( STRING )
        L = LEN_TRIM( BUFFER )

        SPACFLAG = .FALSE.
        NEGVFLAG = .FALSE.
        DO K = 1, L

            CBUF = BUFFER( K:K )

            IF( CBUF .GT. '9' .OR.
     &        ( CBUF .LT. '0' .AND. 
     &          CBUF .NE. ' ' .AND. 
     &          CBUF .NE. '-'       ) .OR.  
     &        ( SPACFLAG .AND. CBUF .EQ. ' ' ) .OR.
     &        ( NEGVFLAG .AND. CBUF .EQ. '-' )           ) THEN

                CHKINT = .FALSE.
                RETURN

            ENDIF

            IF( CBUF .EQ. ' ' ) SPACFLAG = .TRUE.
            IF( CBUF .EQ. '-' ) NEGVFLAG = .TRUE.

        ENDDO    


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END

