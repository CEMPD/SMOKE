
        SUBROUTINE GETMRGEV

C***********************************************************************
C  subroutine GETMRGEV body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
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
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        LOGICAL         ENVYN  

        EXTERNAL CRLF, ENVINT, ENVYN

C...........   Other local variables

        INTEGER         IOS      ! tmp I/O status
        INTEGER         INVIOS   ! i/o status for MEG_REPINV_YN
        INTEGER         SPCIOS   ! i/o status for MEG_REPSPC_YN
        INTEGER         I, J     ! counter

        CHARACTER*5     MRGSRC   ! value of MRG_SOURCE E.V.
        CHARACTER*5     CTLMULT  ! value of MRG_CTLMAT_MULT E.V.
        CHARACTER*5     CTLADD   ! value of MRG_CTLMAT_ADD  E.V.
        CHARACTER*5     CTLREAC  ! value of MRG_CTLMAT_REAC E.V.
        CHARACTER*5     TMPBYDAY ! value of MRG_BYDAY E.V.
        CHARACTER*100   BUFFER   ! text buffer
        CHARACTER*300   MESG     ! message buffer

        CHARACTER*16 :: PROGNAME = 'GETMRGEV' ! program name

C***********************************************************************
C   begin body of subroutine GETMRGEV
        
C.........  Retrieve the environment variables that can be set to contain any
C           combination of "A" "B" "M" "P".
C.........  The default value of these will be blank, so that by not setting 
C           the environment variable(s), nothing will happen with the program

        BUFFER = 'Control for which source categories are in merge'
        CALL ENVSTR( 'MRG_SOURCE', BUFFER, ' ', MRGSRC, IOS  )

        IF( MRGSRC .EQ. ' ' ) THEN
            MESG = 'No source categories selected! You must define ' //
     &             'the MRG_SOURCE ' // CRLF() // BLANK10 // 
     &             'environment variable using the letters A, B, ' //
     &             'M, and P'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        BUFFER = 'Setting for which source categories get ' //
     &           'multiplicative control matrix'
        CALL ENVSTR( 'MRG_CTLMAT_MULT', BUFFER, ' ', CTLMULT, IOS  )

        BUFFER = 'Setting for which source categories get ' //
     &           'additive control matrix'
        CALL ENVSTR( 'MRG_CTLMAT_ADD', BUFFER, ' ', CTLADD, IOS  )

        BUFFER = 'Setting for which source categories get ' //
     &           'reactivity control matrix'
        CALL ENVSTR( 'MRG_CTLMAT_REAC', BUFFER, ' ', CTLREAC, IOS  )

        BUFFER = 'Setting for which source categories get ' //
     &           'by-day processing'
        CALL ENVSTR( 'MRG_BYDAY', BUFFER, ' ', TMPBYDAY, IOS  )

C.........  Process the character strings to set the various control flags 
        AFLAG = ( INDEX( MRGSRC, 'A' ) .GT. 0 )
        BFLAG = ( INDEX( MRGSRC, 'B' ) .GT. 0 )
        MFLAG = ( INDEX( MRGSRC, 'M' ) .GT. 0 )
        PFLAG = ( INDEX( MRGSRC, 'P' ) .GT. 0 )        

        AUFLAG = ( INDEX( CTLMULT, 'A' ) .GT. 0 .AND. AFLAG )
        MUFLAG = ( INDEX( CTLMULT, 'M' ) .GT. 0 .AND. MFLAG ) 
        PUFLAG = ( INDEX( CTLMULT, 'P' ) .GT. 0 .AND. PFLAG )

        AAFLAG = ( INDEX( CTLADD, 'A' ) .GT. 0 .AND. AFLAG )
        MAFLAG = ( INDEX( CTLADD, 'M' ) .GT. 0 .AND. MFLAG )
        PAFLAG = ( INDEX( CTLADD, 'P' ) .GT. 0 .AND. PFLAG )

        ARFLAG = ( INDEX( CTLREAC, 'A' ) .GT. 0 .AND. AFLAG )
        MRFLAG = ( INDEX( CTLREAC, 'M' ) .GT. 0 .AND. MFLAG )
        PRFLAG = ( INDEX( CTLREAC, 'P' ) .GT. 0 .AND. PFLAG )

C.........  Determine if there are multiple source categories
        J = 0
        IF( AFLAG ) J = J + 1
        IF( BFLAG ) J = J + 1
        IF( MFLAG ) J = J + 1
        IF( PFLAG ) J = J + 1
        XFLAG = ( J .GT. 1 )

C.........  Retrieve the on/off environment variables 

        TFLAG   = ENVYN( 'MRG_TEMPORAL_YN', 
     &                   'Use hourly emissions or not', .FALSE., IOS )

        SFLAG   = ENVYN( 'MRG_SPCMAT_YN', 
     &                   'Use speciation matrices or not', .FALSE., IOS)

        LMETCHK = ENVYN( 'MRG_METCHK_YN', 'Check consistency ' //
     &                   'of met headers or not', .TRUE., IOS )

        LMKTPON = ENVYN( 'MRG_MARKETPEN_YN', 'Use market penetration '//
     &                   'from reactivity matrices', .TRUE., IOS )

        LGRDOUT = ENVYN( 'MRG_GRDOUT_YN', 'Output gridded I/O API ' //
     &                   'NetCDF files or not', .FALSE., IOS )

        LREPSTA = ENVYN( 'MRG_REPSTA_YN', 'Output state total ' //
     &                   'ASCII reports or not', .FALSE., IOS )

        LREPCNY = ENVYN( 'MRG_REPCNY_YN', 'Output county total ' //
     &                   'ASCII reports or not', .FALSE., IOS )

        LREPINV = ENVYN( 'MRG_REPINV_YN', 'Report inventory ' //
     &                   'emissions or not', .TRUE., IOS )
        INVIOS = IOS

        LREPSPC = ENVYN( 'MRG_REPSPC_YN', 'Report speciated ' //
     &                   'emissions or not', .TRUE., IOS )
        SPCIOS = IOS

        LREPCTL = ENVYN( 'MRG_REPCTL_YN', 'Report controlled ' //
     &                   'emissions separately or not', .TRUE., IOS )

        LREPANY = ( LREPSTA .OR. LREPCNY )

C.........  Point-source specific environment variables
        IF ( PFLAG ) THEN

            MESG = 'Use layer fractions or not'
            LFLAG   = ENVYN( 'MRG_LAYERS_YN', MESG, .FALSE., IOS )

            MESG = 'Create CMAQ plume-in-grid outputs or not'
            I = ENVINT( 'SMK_PING_METHOD', MESG, .FALSE., IOS )
            PINGFLAG = ( I .EQ. 1 )

            MESG = 'Create ASCII elevated sources file or not'
            ELEVFLAG = ENVYN( 'SMK_ASCIIELEV_YN', MESG, .FALSE., IOS )

            IF ( ELEVFLAG ) PINGFLAG = .FALSE.

            MESG = 'Indicator for including explicit plume ' //
     &             'rise sources'
            EXPLFLAG = ENVYN( 'EXPLICIT_PLUMES_YN', MESG, .FALSE., IOS )

C.............  Must be running for UAM-style processing to use explicit...
            IF( EXPLFLAG .AND. .NOT. ELEVFLAG ) THEN
                ELEVFLAG = .TRUE.
                MESG = 'NOTE: ASCII elevated output switched on to be'//
     &                 'consitent with ' // CRLF() // BLANK10//
     &                 'EXPLICIT_PLUMES_YN = Y setting.'
                CALL M3MSG2( MESG )
            END IF

        END IF

C.........  If temporal processing, set which source categories get by-day 
C           processing
        IF( TFLAG ) THEN
            AFLAG_BD = ( INDEX( TMPBYDAY, 'A' ) .GT. 0 .AND. AFLAG )
            MFLAG_BD = ( INDEX( TMPBYDAY, 'M' ) .GT. 0 .AND. MFLAG )
            PFLAG_BD = ( INDEX( TMPBYDAY, 'P' ) .GT. 0 .AND. PFLAG )
        END IF

C.........  Retrieve variable to indicate whether to use annual or ozone 
C           season data
        MESG = 'Use annual or ozone season emissions'
        LO3SEAS = ENVYN( 'SMK_O3SEASON_YN', MESG, .FALSE., IOS )

C.........  Set index for extracting pollutant data
        INVPIDX = 1
        IF( ( AFLAG .OR. PFLAG ) .AND. 
     &      .NOT. TFLAG          .AND. 
     &      LO3SEAS                    ) INVPIDX = 2

C.........  Check output flags to ensure at least some output
        IF( .NOT. LGRDOUT .AND.
     &      .NOT. LREPSTA .AND.
     &      .NOT. LREPCNY       ) THEN
            MESG = 'No output flags are set!  You must set at least ' //
     &             'one of the' // CRLF() // BLANK10 //
     &             'following to "Y": '//
     &             'MRG_GRDOUT_YN, MRG_REPSTA_YN, or MRG_REPCNY_YN.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Make biogenics flag consitent with speciation flag (must have
C           speciation to be able to include biogenic emissions)
        IF( BFLAG .AND. XFLAG .AND. .NOT. SFLAG ) THEN

            MESG = 'Speciation control environment variable ' //
     &             '"MRG_SPCMAT_YN" indicates no'//
     &             CRLF()// BLANK10// 'speciation, so biogenics will '//
     &             'not be included despite setting of environment ' //
     &             CRLF() // BLANK10 // 'variable "MRG_SOURCE"'
            CALL M3WARN( PROGNAME, 0, 0, MESG )

            BFLAG = .FALSE.

        ELSE IF ( BFLAG ) THEN
            SFLAG = .TRUE.
               
        ENDIF

C.........  Do not output ASCII elevated unless speciation is being done
        IF( ELEVFLAG .AND. .NOT. SFLAG ) THEN
            MESG = 'NOTE: Turning off ASCII elevated point source ' //
     &             'outputs because speciation is turned off.'
            CALL M3MSG2( MESG )
            ELEVFLAG = .FALSE.

        END IF

C.........  Strange to have elevated ASCII and layer merge at the same time 
        IF( ELEVFLAG .AND. LFLAG ) THEN
            MESG = 'WARNING: Elevated ASCII output and 3-d ' //
     &             'output at the same time!'
            CALL M3MSG2( MESG ) 

        END IF

C.........  Cannot have layer merge at the same time as explicit plume rise
        IF( LFLAG .AND. EXPLFLAG ) THEN

            LFLAG = .FALSE.
            MESG = 'WARNING: Turning off layered merge (MRG_LAYERS_YN'//
     &             ' = Y) because' // CRLF() // BLANK10 //
     &             'explicit plume rise being used ' //
     &             '(EXPLICIT_PLUMES_YN = Y).'
            CALL M3MSG2( MESG )

        END IF

C.........  Cannot have BFLAG without TFLAG
        IF( BFLAG ) TFLAG = .TRUE.

C.........  Report that flags for reporting inventory emissions, speciated
C           emissions or not, and controlled emissions or not do not work yet
        IF( INVIOS .GE. 0 ) THEN   ! If it was set to something
            MESG = 'NOTE: MRG_REPINV_YN control is not yet functional.'
            CALL M3MSG2( MESG )
        END IF

        IF( SPCIOS .GE. 0 ) THEN
            MESG = 'NOTE: MRG_REPSPC_YN control is not yet functional.'
            CALL M3MSG2( MESG )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GETMRGEV
