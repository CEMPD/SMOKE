
        SUBROUTINE RDSRCFF10MB( LINE, CFIP, CLNK, TSCC, 
     &                         NVARPERLN, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an ORL format mobile-source inventory
C      file and returns the unique source characteristics.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by  B.H. Baek (Aug 2011)
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
        USE MODLISTS, ONLY: UCASNKEP, NUNIQCAS, UNIQCAS, MXIDAT, INVSTAT,
     &                      INVDNAM, NINVTBL, ITCASA, ITNAMA

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         FINDC, INDEX1
        LOGICAL         CHKINT
        LOGICAL         USEEXPGEO
        
        EXTERNAL   CRLF, FINDC, INDEX1, CHKINT, USEEXPGEO

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT (IN) :: LINE      ! input line
        CHARACTER(FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(LNKLEN3), INTENT(OUT) :: CLNK      ! link ID
        CHARACTER(SCCLEN3), INTENT(OUT) :: TSCC      ! scc code
        INTEGER,            INTENT(OUT) :: NVARPERLN ! no. variables per line
        LOGICAL,            INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,            INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXDATFIL = 60  ! arbitrary max data variables in file
        INTEGER, PARAMETER :: NSEG = 30      ! number of segments in line

C...........   Other local variables
        INTEGER         I       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NVAR    !  number of variables in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine is called
 
        CHARACTER(128) SEGMENT( NSEG ) ! segments of line
        CHARACTER(CASLEN3) TCAS            ! tmp cas number
        CHARACTER(300)     MESG            ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDSRCFF10MB' ! Program name

C***********************************************************************
C   begin body of subroutine RDSRCFF10MB

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXDATFIL, .NOT. USEEXPGEO(), .TRUE., .FALSE., 
     &               LINE, ICC, INY, NVAR, IOS )

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

C.........  Processing activity data
        CLNK = ' '                        ! link ID
        TSCC = SEGMENT( 6 )               ! scc code
        TCAS = ADJUSTL( SEGMENT( 9 ) )

C.........  Determine whether activity data or not
        I = INDEX1( TCAS, NINVTBL, ITCASA )
        I = INDEX1( ITNAMA( I ), MXIDAT, INVDNAM )
        IF( INVSTAT( I ) > 0 ) THEN
            MESG='ERROR: MUST process activity data for FF10_ACTIVITY '
     &          // 'format.' // CRLF() // BLANK10 // TRIM( TCAS ) //
     &          ' is not an activity data based on INVTABLE file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        I = FINDC( TCAS, NUNIQCAS, UNIQCAS )
        IF( I < 1 ) THEN
            NVARPERLN = 0
        ELSE
            NVARPERLN = UCASNKEP( I )
        END IF

C.........  Make sure routine knows it's been called already
        FIRSTIME = .FALSE.

C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        END SUBROUTINE RDSRCFF10MB
