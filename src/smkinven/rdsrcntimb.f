
        SUBROUTINE RDSRCNTIMB( LINE, CFIP, CLNK, TSCC, 
     &                         NVARPERLN, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an NTI format mobile-source inventory
C      file and returns the unique source characteristics.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen (01/03) based on rdntimb.f
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
        USE MODLISTS, ONLY: UCASNKEP, NUNIQCAS, UNIQCAS

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         FINDC
        
        EXTERNAL   CRLF, FINDC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*),       INTENT (IN) :: LINE      ! input line
        CHARACTER(LEN=FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(LEN=LNKLEN3), INTENT(OUT) :: CLNK      ! link ID
        CHARACTER(LEN=SCCLEN3), INTENT(OUT) :: TSCC      ! scc code
        INTEGER,                INTENT(OUT) :: NVARPERLN ! no. variables per line
        LOGICAL,                INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,                INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXDATFIL = 60  ! arbitrary max data variables in file
        INTEGER, PARAMETER :: NSEG = 6       ! number of segments in line

C...........   Other local variables
        INTEGER         I       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NVAR    !  number of variables in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
 
        CHARACTER(LEN=SCCLEN3) SEGMENT( NSEG ) ! segments of line
        CHARACTER(LEN=CASLEN3) TCAS            ! tmp cas number
        CHARACTER(LEN=300)     MESG            ! message buffer

        CHARACTER*16 :: PROGNAME = 'RDSRCNTIMB' ! Program name

C***********************************************************************
C   begin body of subroutine RDSRCNTIMB

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXDATFIL, .TRUE., .TRUE., .FALSE., 
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
        
C.........  Use the file format definition to parse the line into
C           the various data fields
        WRITE( CFIP( 1:1 ), '(I1)' ) ICC  ! country code of FIPS        
        CFIP( 2:3 ) = ADJUSTR( SEGMENT( 1 )( 1:2 ) )  ! state code
        CFIP( 4:6 ) = ADJUSTR( SEGMENT( 2 )( 1:3 ) )  ! county code
        CLNK = ' '                        ! link ID
        TSCC = SEGMENT( 3 )               ! scc code

C.........  Replace blanks with zeros        
        DO I = 1,FIPLEN3
            IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
        END DO
        
C.........  Determine number of pollutants for this line based on CAS number
        TCAS = ADJUSTL( SEGMENT( 4 ) )
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

        END SUBROUTINE RDSRCNTIMB
