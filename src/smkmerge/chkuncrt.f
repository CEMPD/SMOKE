
        SUBROUTINE GETUNCRT( UINAME )

C***********************************************************************
C  subroutine CHKUNCRT body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to indetify a pollutants as 
C      an uncertainty.  If a pollutant is an uncertainty it will be 
C      labelled as such allowing the proper data to be loaded for 
C      merging
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 6/02 by G. Cano
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$ $
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $ $
C Last updated: $ $ 
C
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

C.........  This module contains the global variables for uncertainty
        USE MODUNCERT

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
c        CHARACTER*2     CRLF
c        LOGICAL         DSCM3GRD
c        LOGICAL         ENVYN
c        CHARACTER*50    GETCFDSC  
c        INTEGER         GETIFDSC  
c        INTEGER         PROMPTFFILE  
c        CHARACTER*16    PROMPTMFILE  
c        INTEGER         SECSDIFF  

        INTEGER         INDEX1
        INTEGER         GETIFDSC  
        CHARACTER*16    PROMPTMFILE  

        EXTERNAL  INDEX1, GETIFDSC, PROMPTMFILE

c        EXTERNAL  CRLF, ENVYN, GETCFDSC, GETIFDSC, PROMPTFFILE, 
c     &            PROMPTMFILE, SECSDIFF

C...........   Subroutine arguments
        CHARACTER(*), INTENT(IN)      :: UINAME ! I/O API uncertainty inven input

C.........  Other local variables

        INTEGER         I, J, U, L    ! counters and indices

        INTEGER         IOS           ! tmp I/O status
        INTEGER         ISECS         ! tmp duration in seconds
        INTEGER         NPACT         ! no. variables per activity
        INTEGER         NPPOL         ! no. variables per pollutant
        INTEGER         NDIM          ! tmp dimensioning variable 
        INTEGER         NUVAR         ! tmp no. of uncertainty variables 

        CHARACTER*16 :: VARBUF        ! temporary buffer
        CHARACTER*300   MESG          ! message buffer

        CHARACTER*16 :: PROGNAME = 'GETUNCRT' ! program name

C***********************************************************************
C   begin body of subroutine GETUNCRT

        CALL RDUSTATI( .FALSE., UINAME, ' ', ' ' )

        ALLOCATE( UEANAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UEANAM', PROGNAME )
        UEANAM = .FALSE.     ! Array initialization

        U = 0
C.........  Load emission factors pollutants and activities
        DO I = 1, NIPPA

            L = LEN_TRIM( EANAM( I ) )
            VARBUF = 'MTH_EF_' // EANAM( I )( 1:L )
            J = INDEX1( VARBUF, NVARS3D, VNAME3D)
            IF ( J .GT. 0 ) THEN

                U = U + 1
                UEANAM( J ) = .TRUE.

             END IF

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS

            SUBROUTINE LASTCALL
            END SUBROUTINE LASTCALL

        END SUBROUTINE GETUNCRT
