
        CHARACTER*(*) FUNCTION GETCFDSC( FILEINFO, KEY )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION: 
C     Retreives a character string from the FDESC array
C
C  PRECONDITIONS REQUIRED:
C     
C
C  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
C
C  REVISION  HISTORY:
C       prototype 1/99 by M Houyoux
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
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   ARGUMENTS and their descriptions:

        CHARACTER*(*) FILEINFO( * )! Array of file information
        CHARACTER*(*) KEY          ! Key to search for in FILEINFO

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2  CRLF
        INTEGER      STR2INT
        INTEGER      TRIMLEN

        EXTERNAL     CRLF, STR2INT, TRIMLEN

C...........   LOCAL VARIABLES their descriptions:

        INTEGER       L1, L2       ! length of strings
        INTEGER       I            ! loop index
        INTEGER       K            ! description string position of key
        INTEGER       LENGTH       ! length of function

        CHARACTER*300 CVAL         ! temporary output buffer
        CHARACTER*300 BUFFER       ! key buffer
        CHARACTER*300 MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'GETCFDSC'    ! Program name

C***********************************************************************
C   begin body of subroutine GETCFDSC

C.........  Ensure left-justified keyword
        L1 = TRIMLEN( KEY )
        BUFFER = ADJUSTL( KEY( 1:L1 ) )

        L1 = TRIMLEN( BUFFER )

C.........  Get length of function
        LENGTH = LEN( GETCFDSC )

        DO I = 1, MXDESC3

            K = INDEX( FILEINFO( I ), BUFFER( 1:L1 ) )

            IF( K .GT. 0 ) THEN
                
                L2 = TRIMLEN( FILEINFO( I ) )
                CVAL = FILEINFO( I )( K+L1:L2 )

                L2 = TRIMLEN( CVAL )

                IF( L2 .GT. LENGTH ) THEN
                    WRITE( MESG,94010 )
     &                     'INTERNAL WARNING: Length of string used ' //
     &                     'to call function "' //
     &                     PROGNAME( 1:TRIMLEN( PROGNAME ) ) // '"' //
     &                     CRLF() // BLANK16 //
     &                     'is not long enough. Needs to be at least',
     &                     LENGTH, '. Trimming FDESC3D entry.'
                    CALL M3MSG2( MESG )

                ENDIF

                GETCFDSC = ADJUSTL( CVAL( 1:LENGTH ) )
                RETURN

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
