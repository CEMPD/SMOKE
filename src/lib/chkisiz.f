
        SUBROUTINE CHKISIZ( FILNAM, FILDESC, COMPDESC, NSRC, STATUS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Checks the number of sources and sets error status (0 okay, 1 bad)
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

C...........   INCLUDE FILES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
C       CHARACTER(2)    CRLF 
        INTEGER         TRIMLEN

C        EXTERNAL        CRLF, TRIMLEN
        EXTERNAL     TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*)   FILNAM         ! Name of file being checked
        CHARACTER(*)   FILDESC        ! Description of file being checked
        CHARACTER(*)   COMPDESC       ! Description of comparison value
        INTEGER        NSRC           ! Number of sources comparing against
        INTEGER        STATUS         ! Exit status

C...........   Other local variables

        CHARACTER(256)  MESG 

        CHARACTER(16) :: PROGNAME = 'CHKISIZ' ! program name

C***********************************************************************
C   begin body of subroutine CHKISIZ

        STATUS = 0

        IF( .NOT. DESC3( FILNAM ) ) THEN
            STATUS = 1
            MESG = 'Could not read description for "' //
     &             FILNAM( 1:TRIMLEN( FILNAM ) ) // '"'
            CALL M3MSG2( MESG )
        ENDIF

C.............  Compare the number of sources to the NSRC value
        IF( NCOLS3D .NE. NSRC ) THEN
            STATUS = 1
            WRITE( MESG,94010 )
     &             'Number of inventory sources mismatch.' //
     &             CRLF() // BLANK16 //
     &             COMPDESC // ':', NSRC, 
     &             CRLF() // BLANK16 //
     &             FILDESC // ':', NCOLS3D

            CALL M3MSG2( MESG )
        ENDIF
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END

