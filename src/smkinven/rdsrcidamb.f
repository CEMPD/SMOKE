
        SUBROUTINE RDSRCIDAMB( LINE, CFIP, CLNK, TSCC, 
     &                         NVARPERLN, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an IDA format mobile-source inventory
C      file and returns the unique source characteristics.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen (01/03) based on rdidamb.f
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
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NPACT

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)           CRLF
        
        EXTERNAL   CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT (IN) :: LINE      ! input line
        CHARACTER(FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(LNKLEN3), INTENT(OUT) :: CLNK      ! link ID
        CHARACTER(SCCLEN3), INTENT(OUT) :: TSCC      ! scc code
        INTEGER,            INTENT(OUT) :: NVARPERLN ! no. variables per line
        LOGICAL,            INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,            INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXVARFIL = 112 ! maximum data variables in file

C...........   Local allocatable arrays
        CHARACTER(25), ALLOCATABLE, SAVE :: SEGMENT( : )  ! list-formatted strings
        
C...........   Other local variables
        INTEGER         I, J, K, L, N, V  ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NVAR    !  number of variables in file
        INTEGER         NSEG    ! number of input segments

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
        LOGICAL, SAVE:: FIXED    = .TRUE. ! true: input file is fixed-format
 
        CHARACTER(300)     MESG            ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDSRCIDAMB' ! Program name

C***********************************************************************
C   begin body of subroutine RDSRCIDAMB

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXVARFIL, .TRUE., .TRUE., .TRUE., 
     &               LINE, ICC, INY, NVAR, IOS )

C.........  Interpret error status
        IF( IOS == 4 ) THEN
            WRITE( MESG,94010 ) 
     &             'Maximum allowed data variables ' //
     &             '(MXVARFIL=', MXVARFIL, CRLF() // BLANK10 //
     &             ') exceeded in input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE IF( IOS > 0 ) THEN
            EFLAG = .TRUE.

        END IF

C.........  If a header line was encountered, check for free or fixed format
        IF( IOS >= 0 ) THEN
            
            IF( LINE( 2:5 ) == 'TYPE' ) THEN

C.................  Try to find activity in name of data, otherwise, assume
C                   emissions
                MESG = LINE
                CALL UPCASE( MESG )
            
                IF( INDEX( MESG, 'ACTIVITY' ) > 0 ) THEN
                    FIXED = .FALSE.
                ELSE
                    FIXED = .TRUE.
                END IF
            END IF
            
            HDRFLAG = .TRUE.
            RETURN
        ELSE
            HDRFLAG = .FALSE.
            NVARPERLN = NVAR
        END IF

C.........  If not fixed format, allocate memory for number of segments
        IF( .NOT. FIXED .AND. .NOT. ALLOCATED( SEGMENT ) ) THEN
            NSEG = 4 + NVAR * NPACT
            ALLOCATE( SEGMENT( NSEG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )
            SEGMENT = ' '   ! array
    	ELSE IF ( .NOT. FIXED ) THEN
            SEGMENT = ' '   ! array
        END IF
        
C.........  Use the file format definition to parse the LINE into
C           the various data fields
        WRITE( CFIP( 1:1 ), '(I1)' ) ICC  ! country code of FIPS
        
        IF( FIXED ) THEN
            CFIP( 2:3 ) = ADJUSTR( LINE(  1:2  ) ) ! state
            CFIP( 4:6 ) = ADJUSTR( LINE(  3:5  ) ) ! county
            CLNK        =          LINE(  6:15 )   ! link
            TSCC        =          LINE( 16:25 )   ! scc
        ELSE
            CALL PARSLINE( LINE, NSEG, SEGMENT )        
            CFIP( 2:3 ) = ADJUSTR( SEGMENT( 1 )( 1:2 ) )  ! state code
            CFIP( 4:6 ) = ADJUSTR( SEGMENT( 2 )( 1:3 ) )  ! county code
            CLNK = SEGMENT( 3 )               ! link ID
            TSCC = SEGMENT( 4 )               ! scc code
        END IF

C.........  Replace blanks with zeros        
        DO I = 1,FIPLEN3
            IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
        END DO

        CLNK = ADJUSTL( CLNK )
            
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

        END SUBROUTINE RDSRCIDAMB
