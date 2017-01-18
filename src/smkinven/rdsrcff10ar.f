
        SUBROUTINE RDSRCFF10AR( LINE, CFIP, TSCC, NPOLPERLN,
     &                         HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from a FF10 format nonpoint-source inventory
C      file and returns the unique source characteristics.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by B.H. Baek   Aug 2011
C       Version 10/2016 by C. Coats:  USE M3UTILIO
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
        USE M3UTILIO

C...........   MODULES for public variables
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: UCASNKEP, NUNIQCAS, UNIQCAS

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: DAYINVFLAG, HRLINVFLAG, FF10INVFLAG

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL, EXTERNAL :: CHKINT, USEEXPGEO

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT (IN) :: LINE      ! input line
        CHARACTER(FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(SCCLEN3), INTENT(OUT) :: TSCC      ! scc code
        INTEGER,            INTENT(OUT) :: NPOLPERLN ! no. pollutants per line
        LOGICAL,            INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,            INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters
        INTEGER, PARAMETER :: MXDATFIL = 60  ! arbitrary max no. data variables
        INTEGER, PARAMETER :: NSEG = 60      ! number of segments in line

C...........   Other local variables
        INTEGER         I, L1, L2        ! counters and indices

        INTEGER, SAVE:: ICC      !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY      !  inventory year
        INTEGER         IOS      !  i/o status
        INTEGER, SAVE:: NPOA     !  number of pollutants in file

        CHARACTER(50)      SEGMENT( NSEG ) ! segments of line
        CHARACTER(CASLEN3) TCAS            ! tmp cas number
        CHARACTER(300)     MESG            ! message buffer
        CHARACTER(500)     LBUF            ! line  buffer

        CHARACTER(16) :: PROGNAME = 'RDSRCFF10AR' ! Program name

C***********************************************************************
C   begin body of subroutine RDSRCFF10AR

C.........  Scan for header lines and check to ensure all are set
C           properly (country and year required)
        CALL GETHDR( MXDATFIL, .NOT. USEEXPGEO(), .TRUE., .FALSE.,
     &               LINE, ICC, INY, NPOA, IOS )

C.........  Interpret error status
        IF( IOS == 4 ) THEN
            WRITE( MESG,94010 )
     &             'Maximum allowed data variables ' //
     &             '(MXDATFIL=', MXDATFIL, CRLF() // BLANK10 //
     &             ') exceeded in input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE IF( IOS > 0 ) THEN
            EFLAG = .TRUE.

        END IF

C.........  If a header line was encountered, set flag and return
        IF( IOS >= 0 ) THEN

C.............  Determine whether processing daily/hourly inventories or not
            LBUF = LINE
            CALL UPCASE( LBUF )

            L1 = INDEX( LBUF, 'FF10_DAILY_' )
            L2 = INDEX( LBUF, 'FF10_HOURLY_' )

            IF( INDEX( LBUF, '_POINT'  ) > 0 ) THEN
                MESG = 'ERROR: Can not process POINT inventory '//
     &                 'as AREA inventory'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( FF10INVFLAG ) THEN
              IF( .NOT. DAYINVFLAG .AND. L1  > 0 ) THEN
                MESG = 'ERROR: MUST set DAY_SPECIFIC_YN to Y '//
     &               'to process daily FF10_DAILY inventory'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              END IF

              IF( .NOT. HRLINVFLAG .AND. L2 > 0 ) THEN
                MESG = 'ERROR: MUST set HOUR_SPECIFIC_YN to Y '//
     &               'to process hourly FF10_HOURLY inventory'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

              END IF
            END IF

            HDRFLAG = .TRUE.
            RETURN
        ELSE
            HDRFLAG = .FALSE.
        END IF

C.........  Separate line into segments
        CALL PARSLINE( LINE, NSEG, SEGMENT )

C......... Return if the first line is a header line
        IF( .NOT. CHKINT( SEGMENT( 2 ) ) ) THEN
            HDRFLAG = .TRUE.
            RETURN
        END IF

C.........  Use the file format definition to parse the line into
C           the various data fields
        CFIP = REPEAT( '0', FIPLEN3 )
        IF( USEEXPGEO() ) THEN
            CFIP(  1: 3 ) = ADJUSTR( SEGMENT( 1 )( 1:3 ) )
            CFIP(  4: 9 ) = ADJUSTR( SEGMENT( 2 )( 1:6 ) )
            CFIP( 10:12 ) = ADJUSTR( SEGMENT( 3 )( 1:3 ) )
        ELSE
            WRITE( CFIP( FIPEXPLEN3+1:FIPEXPLEN3+1 ), '(I1)' ) ICC  ! country code of FIPS
            CFIP( FIPEXPLEN3+2:FIPEXPLEN3+6 ) = ADJUSTR( SEGMENT( 2 )( 1:5 ) )  ! state/county code
        END IF

C.........  Replace blanks with zeros
        DO I = 1,FIPLEN3
            IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
        END DO

C.........  Determine number of pollutants for this line based on CAS number
        IF( FF10INVFLAG ) THEN
            TSCC = SEGMENT( 8 )                           ! SCC code
            TCAS = ADJUSTL( SEGMENT( 9 ) )
        ELSE
            TSCC = SEGMENT( 6 )                           ! SCC code
            TCAS = ADJUSTL( SEGMENT( 8 ) )
        END IF

        I = FINDC( TCAS, NUNIQCAS, UNIQCAS )
        IF( I < 1 ) THEN
            NPOLPERLN = 0
        ELSE
            NPOLPERLN = UCASNKEP( I )
        END IF

C.........  Return from subroutine
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        END SUBROUTINE RDSRCFF10AR
