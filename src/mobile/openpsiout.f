
        SUBROUTINE OPENPSIOUT( ENAME, TNAME, FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine opens the intermediate files for processing 
C      parameter scheme indices (PSIs). One file contains the min/max
C      temperature combinations used for each PSI and activity.
C 
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 6/99 by M. Houyoux
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

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER                PROMPTFFILE

        EXTERNAL       PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT    (IN) :: ENAME  ! name of inventory file
        CHARACTER(*), INTENT(IN OUT) :: TNAME  ! name of output tmpr/PSIs file
        INTEGER     , INTENT   (OUT) :: FDEV   ! Unit no. of output tmpr/PSIs

C...........   Other local variables
        CHARACTER*300   MESG    ! message buffer

        CHARACTER*16 :: PROGNAME = 'OPENPSIOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENPSIOUT

C.........  NOTE- This is so simple because source/PSIs file was removed

C.........  Prompt for file with min-max temperatures per PSI
        FDEV  = PROMPTFFILE(
     &          'Enter logical name for PSI/MIN-MAX-TEMPERATURE file',
     &          .FALSE., .TRUE., TNAME, PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94030   FORMAT( A, F10.3, 1X, A )

        END SUBROUTINE OPENPSIOUT


