
        SUBROUTINE GETUNCRT( RLZN, MXGRP, MXVARPGP, NGRP, UINAME )

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
        LOGICAL         SETENVVAR

        EXTERNAL  INDEX1, GETIFDSC, PROMPTMFILE, SETENVVAR

c        EXTERNAL  CRLF, ENVYN, GETCFDSC, GETIFDSC, PROMPTFFILE, 
c     &            PROMPTMFILE, SECSDIFF

C...........   Subroutine arguments
        INTEGER,      INTENT(IN)        :: RLZN     ! realization to process
        INTEGER,      INTENT(IN)        :: NGRP     ! number of groups
        INTEGER,      INTENT(IN)        :: MXVARPGP ! max variables in group
        INTEGER,      INTENT(IN)        :: MXGRP    ! max number of groups
        CHARACTER(*), INTENT(IN OUT)    :: UINAME   ! I/O API uncert inven input

C.........  Other local variables

        INTEGER         I, J, L, N, U, V  ! counters and indices
        INTEGER         L1, L2            ! length counts

        INTEGER         IOS           ! tmp I/O status
        INTEGER         ISECS         ! tmp duration in seconds
        INTEGER         NPACT         ! no. variables per activity
        INTEGER         NPPOL         ! no. variables per pollutant
        INTEGER         NDIM          ! tmp dimensioning variable 
        INTEGER         NUVAR         ! tmp no. of uncertainty variables 

        CHARACTER*16 :: FBUF          ! file name buffer
        CHARACTER*16    LNAME         ! name buffer
        CHARACTER*16 :: TMPBUF        ! temporary buffer
        CHARACTER*16 :: VARBUF        ! temporary buffer
        CHARACTER*300   DIRBUF        ! directory buffer
        CHARACTER*300   MESG          ! message buffer

        CHARACTER*16 :: PROGNAME = 'GETUNCRT' ! program name

C***********************************************************************
C   begin body of subroutine GETUNCRT
 
        IF( RLZN .EQ. 1 ) THEN

            LNAME = ATNAME( 1 )
            IF( UCAT( 1:1 ) .EQ. 'A' ) THEN
                CALL RDIUSRC( NASRC, UINAME )
                DO I = 1, 7
                    L = LEN_TRIM( ATNAME( I ) )
c                    UTNAME( I ) = ATNAME( I )( 1:L ) // 'U'
                    UTNAME( I ) = ATNAME( I )( 1:L ) 
                END DO
            ELSE IF( UCAT( 1:1 ) .EQ. 'P' ) THEN
                CALL RDIUSRC( NPSRC, UINAME )
                DO I = 1, 7
                    L = LEN_TRIM( PTNAME( I ) )
                    UTNAME( I ) = PTNAME( I )( 1:L ) // 'U'
                END DO
            END IF

            ALLOCATE( UEANAM( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'UEANAM', PROGNAME )
            UEANAM = .FALSE.     ! Array initialization
    
            ALLOCATE( NU_EXIST ( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NU_EXIST', PROGNAME )
            NU_EXIST = A_EXIST   ! Array initialization
    
            ALLOCATE( UI_EXIST ( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'UI_EXIST', PROGNAME )
            UI_EXIST = A_EXIST   ! Array initialization

C.............  Load emission factors pollutants and activities
            DO I = 1, NIPPA
                L = LEN_TRIM( EANAM( I ) )
                VARBUF = 'MTH_' // EANAM( I )( 1:L )
                J = INDEX1( VARBUF, NVARS3D, VNAME3D)
                IF ( J .GT. 0 ) UEANAM( I ) = .TRUE.
            END DO

C.........  Load emission factors pollutants and activities
            DO U = 1, NIPPA
                DO N = 1, MXGRP
                    DO V = 1, MXVARPGP
                        IF ( UEANAM( U ) .AND. 
     &                       NU_EXIST( V,N ) .EQ. U ) THEN
C.............................  Uncertainty pollutants to be skipped
                            NU_EXIST( V,N ) = 0
                        ELSE IF ( .NOT.( UEANAM( U ) ) .AND. 
     &                      UI_EXIST( V,N ) .EQ. U ) THEN
C.............................  Uncertainty pollutants to be read
                            UI_EXIST( V,N ) = 0
                        END IF
                    END DO
                END DO
            END DO

        END IF

        L = LEN_TRIM( LNAME )
C.........  For processing an uncertainty realization
        IF( .NOT.( CLOSE3( LNAME ) ) ) THEN
            WRITE( MESG, 94010 ) 
     &             'WARNING: Unable to close file ' 
     &             // LNAME( 1:L ), RLZN
            CALL M3MSG2( MESG )
        END IF

        LNAME = UTNAME( 1 )
        L = LEN_TRIM( LNAME )
C.........  Set up name for temporal input file
        CALL SETTMPU
     
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS

            SUBROUTINE SETTMPU

            INTEGER,  SAVE  ::  L1, L2
            INTEGER                 L3 

            CHARACTER*16,  SAVE ::  FNAM      ! TMPU file name buffer 
            CHARACTER*16            TMPBUF    ! temporary string buffer   
            CHARACTER*300, SAVE ::  DIRBUF    ! directory buffer 
            CHARACTER*300,          PHYSNAME  ! directory buffer 

C**********************************************************************

            IF( RLZN .EQ. 1 ) THEN

C.................  Set name of temporal uncertainty file
                CALL NAMEVAL( UCAT( 1:1 ) // 'TMPU_ODIR', DIRBUF )
                L1 = LEN_TRIM( DIRBUF )

                CALL NAMEVAL( UCAT( 1:1 ) // 'TMPU_FNAM', FNAM )
                L2 = LEN_TRIM( FNAM )

            END IF

            WRITE( TMPBUF,94010 ) '', RLZN
            TMPBUF = ADJUSTL( TMPBUF )
            CALL PADNZERO( RMXLEN, TMPBUF )
            L3 = LEN_TRIM( TMPBUF )

            PHYSNAME = DIRBUF( 1:L1 )  // '/'    // 
     &                 FNAM( 1: L2 )   // '_'    //
     &                 TMPBUF( 1: L3 ) // '.ncf' 

            IF( .NOT.( SETENVVAR( LNAME, PHYSNAME ) ) )
     &           CALL M3EXIT( PROGNAME, 0, 0, 
     &                        'Unable to assign setenv for' // 
     &                         LNAME( 1:L ), 2 )

            LNAME = PROMPTMFILE( 
     &               'Enter logical name for temporal ' //
     &               'unceratainty input file',
     &               FSREAD3, LNAME, PROGNAME )

C...........   Formatted file I/O formats............ 93xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

            END SUBROUTINE SETTMPU

        END SUBROUTINE GETUNCRT
