
        LOGICAL FUNCTION CHKREAL( STRING )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function returns true if the string appears to be a real
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
C Project Title: EDSS Tools Library
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
C***************************************************************************

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
C       CHARACTER(*), INTENT (IN OUT) :: STRING   ! input character string
        CHARACTER(*), INTENT (IN) :: STRING   ! character string used as input only

C...........   Other local variables
        INTEGER         K, L  ! counters and indices

        LOGICAL         EXPOFLAG  ! true if exponential
        LOGICAL         SPACFLAG  ! true if already encountered a space in string
        LOGICAL         PERDFLAG  ! true if already encountered a period in string
        LOGICAL         NEGVFLAG  ! true if already encountered a '-' in string

        CHARACTER       CBUF 
        CHARACTER(256)  BUFFER 

        CHARACTER(16) :: PROGNAME = 'CHKREAL' ! program name

C***********************************************************************
C   begin body of function CHKREAL

        CHKREAL = .TRUE.

        BUFFER = ADJUSTL( STRING )
        L = LEN_TRIM( BUFFER )

        CALL UPCASE( BUFFER )

        EXPOFLAG = .FALSE.
        SPACFLAG = .FALSE.
        PERDFLAG = .FALSE.
        NEGVFLAG = .FALSE.
        DO K = 1, L

            CBUF = BUFFER( K:K )

            IF( 
     &        ( CBUF .GT. '9' .AND. 
     &          CBUF .NE. 'E'       ) .OR.
     &        ( CBUF .LT. '0' .AND. 
     &          CBUF .NE. ' ' .AND. 
     &          CBUF .NE. '.' .AND.
     &          CBUF .NE. '-' .AND.
     &          CBUF .NE. '+'       ) .OR.  
     &        ( ( SPACFLAG .OR. PERDFLAG .OR. NEGVFLAG ) .AND.
     &            CBUF .EQ. ' '                                ) .OR.
     &        ( NEGVFLAG .AND. CBUF .EQ. '-' ) .OR.
     &        ( PERDFLAG .AND. CBUF .EQ. '.' ) .OR.
     &        ( EXPOFLAG .AND. CBUF .EQ. 'E' )      ) THEN

                CHKREAL = .FALSE.
                RETURN

            END IF

            IF( CBUF .EQ. ' ' ) SPACFLAG = .TRUE.
            IF( CBUF .EQ. '.' ) PERDFLAG = .TRUE.
            IF( CBUF .EQ. '-' ) NEGVFLAG = .TRUE.
            IF( CBUF .EQ. 'E' ) THEN
		EXPOFLAG = .TRUE.
                NEGVFLAG = .FALSE.  ! Reinitialize after exponent
            END IF
        END DO    

C       IF( BUFFER .EQ. '.' ) STRING = '0'  ! UNC-IE Jan2025: commented out to not modify STRING value

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END

