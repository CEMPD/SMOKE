
        INTEGER FUNCTION GETINVYR( LINE )

C***********************************************************************
C  function body starts at line 
C
C  DESCRIPTION:
C      This function returns the inventory year from a string if the
C      string contains the INVYEAR packet.
C
C  PRECONDITIONS REQUIRED:

C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 12/98
C
C****************************************************************************/
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

C...........   INCLUDES

        INCLUDE 'PARMS3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         STR2INT
        INTEGER         TRIMLEN

        EXTERNAL        STR2INT, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER*(*) LINE    !  description of source category

        INTEGER         INY
        INTEGER         L1, L2

        CHARACTER*300   MESG        !  message buffer

        CHARACTER*16 :: PROGNAME = 'GETINVYR' ! program name

C***********************************************************************
C   begin body of function GETINVYR

        L1 = INDEX( LINE, 'INVYEAR' )
        L2 = TRIMLEN( LINE )

C.........  Process INVYEAR packet 
        IF( L1 .GT. 0 ) THEN
 
            INY = STR2INT( LINE( L1+7:L2 ) )
 
            IF( INY .LE. 0 ) THEN
 
                MESG = 'Incorrectly set year using ' //
     &                 'INVYEAR packet.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSEIF( INY .LT. 1970 ) THEN

                MESG = 'INVYEAR packet has set 4-digit year ' //
     &                 'below 1970 minimum in PTINV file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 

            ENDIF

            GETINVYR = INY

        ELSE

            GETINVYR = IMISS3

        ENDIF


        RETURN
 

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
