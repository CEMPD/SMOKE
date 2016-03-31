
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY

C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: AFLAG, BFLAG, MFLAG, PFLAG,
     &                      MFLAG_BD, VARFLAG, LREPANY,
     &                      LMETCHK, LGRDOUT, LREPCNY,
     &                      LREPSTA, LREPSCC, LREPSRC,
     &                      SRCGRPFLAG, SMATCHK

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: RPDFLAG, RPHFLAG, RPVFLAG, RPPFLAG, MVFILDIR, TVARNAME,
     &                       SPDFLAG, CFFLAG, EXPCFFLAG, REFCFFLAG, TEMPBIN,
     &                       MOPTIMIZE, GRDENV, TOTENV, MTMP_OUT

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)    CRLF
        REAL            ENVREAL
        LOGICAL         ENVYN  

        EXTERNAL CRLF, ENVREAL, ENVYN

C...........   Other local variables

        INTEGER         IOS      ! tmp I/O status
        INTEGER         I, J, L  ! counter

        CHARACTER(5)    TMPBYDAY ! value of MRG_BYDAY E.V.
        CHARACTER(100)  BUFFER   ! text buffer
        CHARACTER(300)  MESG     ! message buffer

        CHARACTER(16) :: PROGNAME = 'GETMRGEV' ! program name

C***********************************************************************
C   begin body of subroutine GETMRGEV

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  End program if source category is not mobile sources
        IF( CATEGORY /= 'MOBILE' ) THEN
            L = LEN_TRIM( PROGNAME )
            MESG = 'Program ' // PROGNAME( 1:L ) // ' does not ' //
     &             'support ' // TRIM( CATEGORY ) // ' sources.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
C.........  Retrieve the environment variables that can be set to contain any
C           combination of "A" "B" "M" "P".
C.........  The default value of these will be blank, so that by not setting 
C           the environment variable(s), nothing will happen with the program

        BUFFER = 'Setting for which source categories get ' //
     &           'by-day processing'
        CALL ENVSTR( 'MRG_BYDAY', BUFFER, ' ', TMPBYDAY, IOS  )

C.........  Set merge source category flags (impacts library routines shared
C           between Movesmrg and Smkmerge)
        AFLAG = .FALSE.
        BFLAG = .FALSE.
        MFLAG = .TRUE.
        PFLAG = .FALSE.

C.........  Retrieve the on/off environment variables 
        SMATCHK = ENVYN( 'USE_SPCMAT_SPC_YN', 'Use Speciation Matrix '//
     &                   'output file to calculate model species',
     &                   .FALSE., IOS )

        LMETCHK = ENVYN( 'MRG_METCHK_YN', 'Check consistency ' //
     &                   'of met headers or not', .TRUE., IOS )

        LGRDOUT = ENVYN( 'MRG_GRDOUT_YN', 'Output gridded I/O API ' //
     &                   'NetCDF files or not', .FALSE., IOS )

        LREPSTA = ENVYN( 'MRG_REPSTA_YN', 'Output state total ' //
     &                   'ASCII reports or not', .FALSE., IOS )

        LREPCNY = ENVYN( 'MRG_REPCNY_YN', 'Output county total ' //
     &                   'ASCII reports or not', .FALSE., IOS )

        LREPSCC = ENVYN( 'MRG_REPSCC_YN', 'Output SCC total ' //
     &                   'ASCII reports or not', .FALSE., IOS )

        LREPSRC = ENVYN( 'MRG_REPSRC_YN', 'Output source total ' //
     &                   'ASCII reports or not', .FALSE., IOS )

        TEMPBIN = ENVREAL( 'TEMP_BUFFER_BIN', 'Buffer for ' //
     &                   ' min/max temperature ', 0., IOS )

        LREPANY = ( LREPSTA .OR. LREPCNY .OR. LREPSCC .OR. LREPSRC )
                                            
C.........  Retrieve variables for setting the output units for gridded and
C           country/state/county total data
        BUFFER = 'Units for output gridded emissions'
        CALL ENVSTR( 'MRG_GRDOUT_UNIT', BUFFER, ' ', GRDENV, IOS)

        BUFFER = 'Units for output state/county total emissions'
        CALL ENVSTR( 'MRG_TOTOUT_UNIT', BUFFER, ' ', TOTENV, IOS)

        IF( LREPSRC ) THEN
            LREPCNY = .FALSE.
            LREPSTA = .FALSE.
            LREPSCC = .FALSE.
            MESG = 'WARNING: Reset MRG_REPSCC_YN, MRG_REPCNY, and ' //
     &           'MRG_REPSTA_YN to N since MRG_REPSRC_YN is set to Y '
            CALL M3MESG( MESG )
        END IF

C.........  Check for variable grid
        VARFLAG = ENVYN( 'USE_VARIABLE_GRID', 'Use variable grid ' //
     &                 'definition', .FALSE., IOS )

C.........  Check if source grouping should be used
        SRCGRPFLAG = ENVYN( 'SMK_SRCGROUP_OUTPUT_YN', 'Use source ' //
     &                      'grouping', .FALSE., IOS )

C.........  Check for rate-per-distance, rate-per-vehicle, or rate-per-profile processing
        RPDFLAG = ENVYN( 'RPD_MODE', 'Calculate rate-per-distance ' //
     &                   'emissions', .FALSE., IOS )

        RPHFLAG = ENVYN( 'RPH_MODE', 'Calculate rate-per-hour ' //
     &                   'emissions', .FALSE., IOS )

        RPVFLAG = ENVYN( 'RPV_MODE', 'Calculate rate-per-vehicle ' //
     &                   'emissions', .FALSE., IOS )
     
        RPPFLAG = ENVYN( 'RPP_MODE', 'Calculate rate-per-profile ' //
     &                   'emissions', .FALSE., IOS )

        MOPTIMIZE = ENVYN( 'MEMORY_OPTIMIZE_YN', 'Optimize Memory usage' //
     &                    " ", .FALSE., IOS )

        MTMP_OUT = ENVYN( 'MTMP_OUTPUT_YN', 'Output mobile hourly emissions' //
     &                    " ", .FALSE., IOS )

        MESG = 'Apply control factors to emissions'
        CFFLAG = ENVYN( 'USE_CONTROL_FACTORS', MESG, .FALSE., IOS )

        IF( CFFLAG ) THEN
            MESG = 'Use pollutant/species-specific control factor'
            EXPCFFLAG = ENVYN( 'USE_EXP_CONTROL_FAC_YN', MESG, .FALSE., IOS )

            MESG = 'Use reference county-specific control factor'
            REFCFFLAG = ENVYN( 'USE_REF_CONTROL_FAC_YN', MESG, .FALSE., IOS )
        END IF

        IF( .NOT. RPDFLAG .AND.
     &      .NOT. RPHFLAG .AND.
     &      .NOT. RPVFLAG .AND.
     &      .NOT. RPPFLAG ) THEN
            MESG = 'No mode selected!  You must set either RPD_MODE, ' //
     &             'RPH_MODE, RPV_MODE, or RPP_MODE to "Y".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Select default processing mode
        IF( RPDFLAG ) THEN
            RPHFLAG = .FALSE.
            RPVFLAG = .FALSE.
            RPPFLAG = .FALSE.
        ELSE IF( RPVFLAG ) THEN
            RPHFLAG = .FALSE.
            RPPFLAG = .FALSE.
        ELSE IF( RPPFLAG ) THEN
            RPVFLAG = .FALSE.
        END IF

C.........  Check if hourly speeds should be used
        IF( RPDFLAG ) THEN
            SPDFLAG = ENVYN( 'USE_HOURLY_SPEEDS', 'Use hourly ' //
     &                       'speed data', .FALSE., IOS )
        END IF

C.........  Get directory where MOVES output files are stored
        BUFFER = 'Location of MOVES output files'
        CALL ENVSTR( 'SMK_MVSPATH', BUFFER, '.', MVFILDIR, IOS )
        
        IF( IOS /= 0 ) THEN
            MESG = 'ERROR: Could not read value for SMK_MVSPATH.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  If temporal processing, set which source categories get by-day 
C           processing
        MFLAG_BD = ( INDEX( TMPBYDAY, 'M' ) .GT. 0 )

C.........  Get name of temperature variable to read from meteorology files
        BUFFER = 'Name of temperature variable'
        CALL ENVSTR( 'TVARNAME', BUFFER, 'TEMP2', TVARNAME, IOS )

C.........  Check output flags to ensure at least some output
        IF( .NOT. LGRDOUT .AND.
     &      .NOT. LREPSTA .AND.
     &      .NOT. LREPCNY .AND.
     &      .NOT. LREPSRC .AND.
     &      .NOT. LREPSCC       ) THEN
            MESG = 'No output flags are set!  You must set at least ' //
     &             'one of the' // CRLF() // BLANK10 //
     &             'following to "Y": '//
     &             'MRG_GRDOUT_YN, MRG_REPSTA_YN, MRG_REPCNY_YN, ' //
     &             'MRG_REPSCC_YN, or MRG_REPSRC_YN.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GETMRGEV
