C copied by: mhouyoux
C origin: initem.F 1.2

        SUBROUTINE INITEM( LDEV )

***********************************************************************
C  program body starts at line 60
C
C  DESCRIPTION:
C       Writes out abridged copyright information
C
C  PRECONDITIONS REQUIRED:
C       None
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Prototype  02/98 by M Houyoux
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

C.........  Subroutine arguments
        INTEGER       LDEV  ! Log file unit number

C.........  Parameters
        INTEGER       YEAR
        PARAMETER   ( YEAR = 1998 )

C.........  External functions
        INTEGER       TRIMLEN
        EXTERNAL      TRIMLEN

C.........  Local variables
        CHARACTER*256 LINE0, LINE1, LINE2, LINE3

C***********************************************************************
C   begin body of program INITEM

        LINE0 = 'SMOKE ---------------' 
        WRITE( LINE1,94020 ) 'Copyright (c)', YEAR, 
     &                      'MCNC--North Carolina Supercomputing Center'

        LINE2 = 'All right reserved'
        LINE3 = 'See file COPYRIGHT for conditions of use.'

C.........  Write to log file

        IF( LDEV .NE. 6 ) THEN
            WRITE( LDEV,92000 ) LINE0( 1:TRIMLEN( LINE0 ) )
            WRITE( LDEV,92000 ) LINE1( 1:TRIMLEN( LINE1 ) )
            WRITE( LDEV,92000 ) LINE2( 1:TRIMLEN( LINE2 ) )
            WRITE( LDEV,92000 ) 
            WRITE( LDEV,92000 ) LINE3( 1:TRIMLEN( LINE3 ) )
            WRITE( LDEV,92000 ) 

        ENDIF

C.........  Write to screen

        WRITE( *,92000 ) LINE0( 1:TRIMLEN( LINE0 ) )
        WRITE( *,92000 ) LINE1( 1:TRIMLEN( LINE1 ) )
        WRITE( *,92000 ) LINE2( 1:TRIMLEN( LINE2 ) )
        WRITE( *,92000 ) 
        WRITE( *,92000 ) LINE3( 1:TRIMLEN( LINE3 ) )
        WRITE( *,92000 ) 

        RETURN 

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx
 
92000   FORMAT( 5X, A )

94020   FORMAT( 10( A, 1X, :, I4, :, 1X ) )

        END
