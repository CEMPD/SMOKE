        INTEGER FUNCTION GETIFDSC( FILEINFO, KEY, REQUIRED )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION: 
C     Retreives an integer from the FDESC array or returns -1 if key is not
C     required to be found and key is not present
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C****************************************************************************
        USE M3UTILIO

        IMPLICIT NONE

C...........   Include files
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   ARGUMENTS and their descriptions:

        CHARACTER(*), INTENT (IN) :: FILEINFO( * ) ! Array of file information
        CHARACTER(*), INTENT (IN) :: KEY           ! Key to find in FILEINFO
        LOGICAL     , INTENT (IN) :: REQUIRED      ! true: key must be found


C...........   EXTERNAL FUNCTIONS:
C       INTEGER    STR2INT
        INTEGER    TRIMLEN

C        EXTERNAL   STR2INT, TRIMLEN
        EXTERNAL     TRIMLEN

C...........   LOCAL VARIABLES their descriptions:

        INTEGER       L1, L2       ! length of strings
        INTEGER       I            ! loop index
        INTEGER       IVAL         ! temporary integer value
        INTEGER       K            ! description string position of key

        CHARACTER(300) BUFFER      ! Key buffer
        CHARACTER(300) MESG        ! Message buffer

        CHARACTER(16) :: PROGNAME = 'GETIFDSC'    ! Program name

C***********************************************************************
C   begin body of function GETIFDSC

        L1 = TRIMLEN( KEY )
        BUFFER = ADJUSTL( KEY( 1:L1 ) )

        L1 = TRIMLEN( BUFFER )

        DO I = 1, MXDESC3

            K = INDEX( FILEINFO( I ), BUFFER( 1:L1 ) )

            IF( K .GT. 0 ) THEN

                L2 = TRIMLEN( FILEINFO( I ) )
                IVAL = STR2INT( FILEINFO( I )( K+L1:L2 ) )

                IF( IVAL .EQ. IMISS3 ) THEN
                    MESG = 'ERROR: non-integer result found at FDESC '//
     &                     'entry "'// KEY( 1:L1 )// '" in NetCDF file'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ELSE
                    GETIFDSC = IVAL
                    RETURN

                ENDIF
            ENDIF

        ENDDO

C.........  If we get here, then key was not found in FDESC, so if it was
C           required, then abort.

        IF( REQUIRED ) THEN
            MESG = 'FDESC3D packet "' // KEY( 1:L1 ) // 
     &             '" was not found in NetCDF file!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE
            GETIFDSC = -1
            RETURN

        END IF
    
C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

        END FUNCTION GETIFDSC
