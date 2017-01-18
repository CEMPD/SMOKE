
        SUBROUTINE RDSRCORLFR( LINE, CFIP, FIREID, LOCID, SKID,
     &                         SGID, TSCC, NDATPERLN, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an ORL FIRE format wildfire
C      point-source inventory file and returns the unique source characteristics.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by B.H. Baek (02/06) based on rdsrcntifpt.f
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
C***************************************************************************

C...........   MODULES for public variables
C.........  This module contains the lists of unique inventory information
        USE MODINFO, ONLY: TMPNAM

        USE MODLISTS, ONLY: UCASNKEP, NUNIQCAS, UNIQCAS, NINVTBL,
     &                      ITNAMA, ITCASA

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         FINDC
        INTEGER         INDEX1

        EXTERNAL    CRLF, FINDC, INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT (IN) :: LINE      ! input line
        CHARACTER(FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(PLTLEN3), INTENT(OUT) :: FIREID    ! fire ID
        CHARACTER(CHRLEN3), INTENT(OUT) :: LOCID     ! location ID
        CHARACTER(CHRLEN3), INTENT(OUT) :: SKID      ! dummy stack ID
        CHARACTER(CHRLEN3), INTENT(OUT) :: SGID      ! dummy segment ID
        CHARACTER(SCCLEN3), INTENT(OUT) :: TSCC      ! dummy scc code
        INTEGER,            INTENT(OUT) :: NDATPERLN ! no. pollutants per line
        LOGICAL,            INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,            INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 1000  ! arbitrary maximum pollutants in file
        INTEGER, PARAMETER :: NSEG = 10      ! number of segments in line
        INTEGER, PARAMETER :: NEXTRA  = 4    ! number of extra non-data fields that need

C...........   Other local variables
        INTEGER         I, J       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NDAT = -1 !  number of pollutants in file
 
        CHARACTER(CHRLEN3) SEGMENT( NSEG ) ! segments of line
        CHARACTER(300)     MESG            ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDSRCORLFR' ! Program name

C***********************************************************************
C   begin body of subroutine RDSRCORLFR

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .TRUE., 
     &               LINE, ICC, INY, NDAT, IOS )

C.........  Interpret error status
        IF( IOS == 4 ) THEN
            WRITE( MESG,94010 ) 
     &             'Maximum allowed data variables ' //
     &             '(MXPOLFIL=', MXPOLFIL, CRLF() // BLANK10 //
     &             ') exceeded in input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE IF( IOS > 0 ) THEN
            EFLAG = .TRUE.

        END IF

C.........  If a header line was encountered, set flag and return
        IF( IOS >= 0 ) THEN
            HDRFLAG = .TRUE.
            IF( NDAT > 0 ) THEN
                NDATPERLN = NDAT + NEXTRA  
                DO I = 1, NDAT
                    J = FINDC( TMPNAM( I ), NUNIQCAS, UNIQCAS )
                    IF( J > 1 ) NDATPERLN = NDATPERLN + UCASNKEP( J ) - 1
                END DO
            END IF
            RETURN
        ELSE
            HDRFLAG = .FALSE.
        END IF
        
C.........  Give error if #DATA line has not defined the pollutants that
C           are contained in the day-specific data file. NOTE: This code
C           will not be reached until after the last header line)
        IF ( NDAT < 0 ) THEN
            WRITE( MESG, 94010 ) 'First data line reached in ORL '//
     &             'FIRE file with required #DATA header found.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Separate line into segments
        CALL PARSLINE( LINE, NSEG, SEGMENT )

C.........  Use the file format definition to parse the line into
C           the various data fields
        CFIP = REPEAT( '0', FIPLEN3 )
        WRITE( CFIP( FIPEXPLEN3+1:FIPEXPLEN3+1 ), '(I1)' ) ICC  ! country code of FIPS
        CFIP( FIPEXPLEN3+2:FIPEXPLEN3+6 ) = ADJUSTR( SEGMENT( 1 )( 1:5 ) )  ! state/county code

        FIREID = SEGMENT( 2 )   ! fire ID
        LOCID  = SEGMENT( 3 )   ! location ID
        TSCC   = SEGMENT( 4 )   ! scc code
        SKID   = ' '            ! dummy stack ID
        SGID   = ' '            ! dummy segment ID
        
C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        END SUBROUTINE RDSRCORLFR
