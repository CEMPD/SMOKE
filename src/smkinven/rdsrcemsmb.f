
        SUBROUTINE RDSRCEMSMB( LINE, CFIP, CRWT, CLNK, 
     &                         NVARPERLN, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an EMS format mobile-source inventory
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
C***************************************************************************

C...........   MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NPACT

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        LOGICAL                ENVYN
        
        EXTERNAL   CRLF, ENVYN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*),       INTENT (IN) :: LINE      ! input line
        CHARACTER(LEN=FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(LEN=RWTLEN3), INTENT(OUT) :: CRWT      ! roadway type
        CHARACTER(LEN=LNKLEN3), INTENT(OUT) :: CLNK      ! link ID
        INTEGER,                INTENT(OUT) :: NVARPERLN ! no. variables per line
        LOGICAL,                INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,                INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXDATFIL = 60 ! arbitrary max data variables in file

C...........   Local allocatable arrays
        CHARACTER(LEN=25), ALLOCATABLE :: SEGMENT( : )  ! list-formatted strings
        
C...........   Other local variables
        INTEGER         I       ! counters and indices

        INTEGER, SAVE:: FMTCASE !  code for format case
        INTEGER, SAVE:: ICC = 0 !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPRECOL !  no. src char columns for list-directed
        INTEGER, SAVE:: NVAR    !  number of variables in file
        INTEGER         NSEG    ! number of input segments

        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine is called
        LOGICAL, SAVE:: FIXED    = .FALSE. ! true: input file is fixed-format
        LOGICAL, SAVE:: LFLAG    = .FALSE. ! true: link file
 
        CHARACTER(LEN=300)     MESG            ! message buffer

        CHARACTER*16 :: PROGNAME = 'RDSRCEMSMB' ! Program name

C***********************************************************************
C   begin body of subroutine RDSRCEMSMB

        IF( FIRSTIME ) THEN
            MESG = 'Indicator for fixed-column EMS-95 format'
            FIXED = ENVYN( 'SMK_EMS95_FIXFMT', MESG, FIXED, IOS )
            
            FMTCASE = 0
            
            FIRSTIME = .FALSE.
        END IF

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXDATFIL, .FALSE., .FALSE., .TRUE., 
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

C.........  If a header line was encountered, check if file is non-link or link
        IF( IOS >= 0 ) THEN
            
            IF( LINE( 2:LEN_TRIM( LINE ) ) == 'LINK' ) THEN
                LFLAG = .TRUE.
                
                IF( FIXED ) THEN
                    FMTCASE = 3      ! link fixed
                ELSE
                    FMTCASE = 4      ! link list
                    NPRECOL = 8
                END IF
            ELSE IF( LINE( 2:LEN_TRIM( LINE ) ) == 'NONLINK' ) THEN
                LFLAG = .FALSE.
                
                IF( FIXED ) THEN
                    FMTCASE = 1       ! non-link fixed
                ELSE
                    FMTCASE = 2       ! non-link list
                    NPRECOL = 2
                END IF
            END IF
            
            HDRFLAG = .TRUE.
            RETURN
        ELSE
            HDRFLAG = .FALSE.
            NVARPERLN = NVAR
        END IF

C.........  Make sure file format has been set
        IF( FMTCASE == 0 ) THEN
            MESG = 'INTERNAL ERROR: Mobile EMS-95 file format ' //
     &             'is not recognized'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  If not fixed format, allocate memory for number of segments
        IF( .NOT. FIXED .AND. .NOT. ALLOCATED( SEGMENT ) ) THEN
            NSEG = NPRECOL + NVAR
            ALLOCATE( SEGMENT( NSEG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )
            SEGMENT = ' '   ! array
        END IF
        
C.........  Use the file format definition to parse the LINE into
C           the various data fields        
        IF( FIXED ) THEN
            WRITE( CFIP( 1:1 ), '(I1)' ) ICC  ! country code of FIPS
            CFIP( 2:3 ) = LINE( 1:2 )  ! state
            CFIP( 4:6 ) = LINE( 3:5 )  ! county
            CRWT( 1:1 ) = LINE( 6:6 )  ! area type
            CRWT( 2:3 ) = ADJUSTL( LINE( 7:10 ) ) ! facility type  <-- PROBLEM!!!!
            CLNK        = ' '          ! link
        ELSE
            CALL PARSLINE( LINE, NSEG, SEGMENT )
            
            CFIP = ADJUSTR( SEGMENT( 1 )( 1:FIPLEN3 ) )  ! fips code
            CRWT = ADJUSTL( SEGMENT( 2 ) )               ! roadway type
            
            IF( LFLAG ) THEN
                CLNK = ADJUSTL( SEGMENT( 3 ) )   ! link ID
            ELSE
                CLNK = ' '
            END IF
        END IF

C.........  Replace blanks with zeros        
        DO I = 1,FIPLEN3
            IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
        END DO
            
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

        END SUBROUTINE RDSRCEMSMB
