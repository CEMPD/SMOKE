
        INTEGER FUNCTION GETIFDSC( FILEINFO, KEY )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION: 
C     Retreives an integer from the FDESC array
C
C  PRECONDITIONS REQUIRED:
C     
C
C  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
C
C  REVISION  HISTORY:
C       prototype 12/98 by M Houyoux
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

C...........   Include files
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   ARGUMENTS and their descriptions:

        CHARACTER*(*) FILEINFO( * )! Array of file information
        CHARACTER*(*) KEY          ! Key to search for in FILEINFO

C...........   EXTERNAL FUNCTIONS:
        INTEGER    STR2INT
        INTEGER    TRIMLEN

        EXTERNAL   STR2INT, TRIMLEN

C...........   LOCAL VARIABLES their descriptions:

        INTEGER       L1, L2       ! length of strings
        INTEGER       I            ! loop index
        INTEGER       IVAL         ! temporary integer value
        INTEGER       K            ! description string position of key

        CHARACTER*256 MESG         ! Message buffer

        CHARACTER*16 :: PROGNAME = 'GETIFDSC'    ! Program name

C***********************************************************************
C   begin body of subroutine GETIFDSC

        L1 = TRIMLEN( KEY )

        DO I = 1, MXDESC3

            K = INDEX( FILEINFO( I ), KEY( 1:L1 ) )

            IF( K .GT. 0 ) THEN

                L2 = TRIMLEN( FILEINFO( I ) )
                IVAL = STR2INT( FILEINFO( I )( K+L1-1:L2 ) )

                IF( IVAL .EQ. IMISS3 ) THEN
                    MESG = 'ERROR: non-integer result found at FDESC '//
     &                     'entry "'// KEY( 1:L1 )// '" in NetCDF file'
                    CALL M3EXIT( MESG, 0, 0, PROGNAME, 2 )

                ELSE
                    GETIFDSC = IVAL
                    RETURN

                ENDIF
            ENDIF

        ENDDO

C.........  If we get here, then key was not found in FDESC!

        MESG = 'ERROR: key "' // KEY( 1:L1 ) // 
     &         '" not found in NetCDF file!'
        CALL M3EXIT( MESG, 0, 0, PROGNAME, 2 )
    
C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

        END
