
        SUBROUTINE FMTCSRC( INSTRING, NCHARS, OUTBUFF, LENOUT )

C***********************************************************************
C  subroutine body starts at line 78
C
C  DESCRIPTION:
C      This subroutine formats the INSTRING assuming that is a string of 
C      source characteristics that has proper segment lengths.  It adds
C      labels and newlines as necessary and returns the formatted OUTBUFF
C      and its length.
C
C  PRECONDITIONS REQUIRED:
C      INSTRING with correct segement lengths
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Functions: I/O API functions
C
C  REVISION  HISTORY:
C      Created 11/98 by M. Houyoux
C
C**************************************************************************
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2     CRLF
        INTEGER         TRIMLEN

        EXTERNAL        CRLF, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER*(*) INSTRING     ! Input string
        INTEGER       NCHARS       ! No. of characteristics to format (1 to 8)
        CHARACTER*(*) OUTBUFF      ! Formatted output string
        INTEGER       LENOUT       ! Length of output string

C...........   Other local variables
        INTEGER         I, L, K, L1, L2, L3  ! counters and indices

        CHARACTER*5  :: LABEL( 8 ) =  !  message buffer
     &                ( / 'FIPS ', 'PLANT', 'CHAR1', 'CHAR2',
     &                    'CHAR3', 'CHAR4', 'CHAR5', 'POL  ' / )
 
        CHARACTER*300   BUFFER  !  string buffer

        CHARACTER*16 :: PROGNAME = 'FMTCSRC'  ! program name

C***********************************************************************
C   begin body of subroutine FMTCSRC

C.........  Make sure not to exceed NCHARS legal value
        NCHARS = MIN( NCHARS, 8 )

C.........  Initialize output buffer
        IF( NCHARS .GE. 1 ) THEN
            L1 = LENARR3( 1 )
            L2 = LENARR3( 2 ) - 1
  
            WRITE( OUTBUFF, 94900 ) LABEL( 1 ), INSTRING( L1:L2 )

        ENDIF

C.........  Loop through the remaining source chars, writing the output string
C           for each populated source characterstic
        K = TRIMLEN( OUTBUFF )
        DO I = 2, NCHARS

            L = TRIMLEN( OUTBUFF )

            L1 = LENARR3( I )    ! retrieve stored field lengths
            L2 = LENARR3( I+1 ) - 1
 
            BUFFER = ADJUSTL( INSTRING( L1:L2 ) )   ! Left-justifying
            L3 = TRIMLEN( BUFFER )
            K = K + 8 + L3

C.............  Continue to contribute to buffer if not blank
            IF( BUFFER .NE. ' ' ) THEN

C.................  Include a return and indent if line gets too long
                IF( K .GT. EMOUTLN3 ) THEN
                    K = 18 + L3
                    WRITE( OUTBUFF, 94900 ) OUTBUFF( 1:L ) //
     &                                      CRLF() // BLANK10 //
     &                                      LABEL( I ), BUFFER( 1:L3 )
                ELSE
                    WRITE( OUTBUFF, 94900 ) OUTBUFF( 1:L ) // ' ' //
     &                                      LABEL( I ), BUFFER( 1:L3 )

                ENDIF
            ENDIF

        ENDDO

        LENOUT = TRIMLEN( OUTBUFF )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94900   FORMAT( A, ': ', A )

        END





