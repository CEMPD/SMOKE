
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        LOGICAL         ENVYN  

        EXTERNAL CRLF, ENVYN

C...........   Other local variables

        INTEGER         IOS      ! tmp I/O status
        INTEGER         INVIOS   ! i/o status for MEG_REPINV_YN
        INTEGER         SPCIOS   ! i/o status for MEG_REPSPC_YN
        INTEGER         CTLIOS   ! i/o status for MEG_REPCTL_YN
        INTEGER         J        ! counter

        CHARACTER*4     MRGSRC   ! value of MRG_SOURCE E.V.
        CHARACTER*4     CTLMULT  ! value of MRG_CTLMAT_MULT E.V.
        CHARACTER*4     CTLADD   ! value of MRG_CTLMAT_ADD  E.V.
        CHARACTER*4     CTLREAC  ! value of MRG_CTLMAT_REAC E.V.
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

        BUFFER = 'Control for which source categories get ' //
     &           'multiplicative control matrix'
        CALL ENVSTR( 'MRG_CTLMAT_MULT', BUFFER, ' ', CTLMULT, IOS  )

        BUFFER = 'Control for which source categories get ' //
     &           'additive control matrix'
        CALL ENVSTR( 'MRG_CTLMAT_ADD', BUFFER, ' ', CTLADD, IOS  )

        BUFFER = 'Control for which source categories get ' //
     &           'reactivity control matrix'
        CALL ENVSTR( 'MRG_CTLMAT_REAC', BUFFER, ' ', CTLREAC, IOS  )

C.........  Process the character strings to set the various control flags 
        AFLAG = ( INDEX( MRGSRC, 'A' ) .GT. 0 )
        BFLAG = ( INDEX( MRGSRC, 'B' ) .GT. 0 )
        MFLAG = ( INDEX( MRGSRC, 'M' ) .GT. 0 )
        PFLAG = ( INDEX( MRGSRC, 'P' ) .GT. 0 )        

        AUFLAG = ( INDEX( CTLMULT, 'A' ) .GT. 0 )
        MUFLAG = ( INDEX( CTLMULT, 'M' ) .GT. 0 )
        PUFLAG = ( INDEX( CTLMULT, 'P' ) .GT. 0 )

        AAFLAG = ( INDEX( CTLADD, 'A' ) .GT. 0 )
        MAFLAG = ( INDEX( CTLADD, 'M' ) .GT. 0 )
        PAFLAG = ( INDEX( CTLADD, 'P' ) .GT. 0 )

        ARFLAG = ( INDEX( CTLREAC, 'A' ) .GT. 0 )
        MRFLAG = ( INDEX( CTLREAC, 'M' ) .GT. 0 )
        PRFLAG = ( INDEX( CTLREAC, 'P' ) .GT. 0 )

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

        LFLAG   = ENVYN( 'MRG_LAYERS_YN', 
     &                   'Use layer fractions or not', .FALSE., IOS )

        PINGFLAG = ENVYN( 'SMK_PING_YN',  'Create plume-in-grid ' //
     &                    'outputs or not', .FALSE., IOS )

        VFLAG   = ENVYN( 'MRG_VMT_YN', 
     &                   'Use VMT or not', .FALSE., IOS )

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
        CTLIOS = IOS

        LREPANY = ( LREPSTA .OR. LREPCNY )

C.........  Retrieve variable to indicate whether to use annual or ozone 
C           season data
        MESG = 'Use annual or ozone season emissions'
        LO3SEAS = ENVYN( 'SMK_O3SEASON_YN', MESG, .FALSE., IOS )

C.........  Set index for extracting pollutant data
        INVPIDX = 1
        IF( .NOT. TFLAG .AND. LO3SEAS ) INVPIDX = 2

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

C.........  Don't output gridded if gridded biogenics is only input
        IF( BFLAG .AND. .NOT. XFLAG .AND. LGRDOUT ) THEN

            MESG = 'NOTE: Turning off gridded outputs because no ' //
     &             'merging is taking place'
            CALL M3MSG2( MESG )
            LGRDOUT = .FALSE.

        END IF

C.........  Cannot have BFLAG without TFLAG
        IF( BFLAG ) TFLAG = .TRUE.

C.........  Make VMT usage flag consistent with speciation and source merging
C           such that if the merge VMT flag is set to true, speciation and other
C           source categories will be overridden.
        IF( MFLAG .AND. VFLAG .AND. 
     &    ( AFLAG .OR. BFLAG .OR. PFLAG .OR. SFLAG ) ) THEN

            MESG = 'VMT control environment variable "MRG_VMT_YN"' //
     &             'indicates VMT should be used.  This will'//
     &             CRLF()// BLANK10// 'override settings of ' //
     &             'environment variables "MRG_SOURCE" and ' //
     &             '"MRG_SPCMAT_YN"'
            CALL M3WARN( PROGNAME, 0, 0, MESG )

            AFLAG = .FALSE.
            BFLAG = .FALSE.
            PFLAG = .FALSE.
            SFLAG = .FALSE.

        ENDIF

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

        IF( CTLIOS .GE. 0 ) THEN
            MESG = 'NOTE: MRG_REPCTL_YN control is not yet functional.'
            CALL M3MSG2( MESG )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GETMRGEV
