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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES

C        INCLUDE 'PARMS3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
C       INTEGER         STR2INT
        INTEGER         TRIMLEN

C        EXTERNAL        STR2INT, TRIMLEN
        EXTERNAL     TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*) LINE    !  description of source category

        INTEGER         INY
        INTEGER         L1, L2

        CHARACTER(300)  MESG        !  message buffer

        CHARACTER(16) :: PROGNAME = 'GETINVYR' ! program name

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
