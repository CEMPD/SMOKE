
        SUBROUTINE RDSRCORLFR( LINE, NN, CFIP, FIREID, LOCID, SKID,
     &                         SGID, TSCC, NPOLPERLN, HDRFLAG, EFLAG,
     &                         BKSPFLAG )

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
        USE MODLISTS, ONLY: UCASNKEP, NUNIQCAS, UNIQCAS, NINVTBL,
     &                      ITNAMA, ITCASA

C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: FIREPOL

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
        INTEGER,            INTENT (IN) :: NN        ! no of current FIREPOL
        CHARACTER(FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(PLTLEN3), INTENT(OUT) :: FIREID    ! fire ID
        CHARACTER(CHRLEN3), INTENT(OUT) :: LOCID     ! location ID
        CHARACTER(CHRLEN3), INTENT(OUT) :: SKID      ! dummy stack ID
        CHARACTER(CHRLEN3), INTENT(OUT) :: SGID      ! dummy segment ID
        CHARACTER(SCCLEN3), INTENT(OUT) :: TSCC      ! dummy scc code
        INTEGER,            INTENT(OUT) :: NPOLPERLN ! no. pollutants per line
        LOGICAL,            INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,            INTENT(OUT) :: EFLAG     ! error flag
        LOGICAL,            INTENT (IN) :: BKSPFLAG  ! backspace flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 60  ! arbitrary maximum pollutants in file
        INTEGER, PARAMETER :: NSEG = 63      ! number of segments in line

C...........   Other local variables
        INTEGER         I, II       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
 
        CHARACTER(CHRLEN3) SEGMENT( NSEG ) ! segments of line
        CHARACTER(CASLEN3) TCAS            ! tmp cas number
        CHARACTER( 2 )     ID              ! fake SCC ID for extra pol
        CHARACTER(300)     MESG            ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDSRCORLFR' ! Program name

C***********************************************************************
C   begin body of subroutine RDSRCORLFR

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .FALSE., 
     &               LINE, ICC, INY, NPOL, IOS )

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
            RETURN
        ELSE
            HDRFLAG = .FALSE.
        END IF
        
C.........  Separate line into segments
        CALL PARSLINE( LINE, NSEG, SEGMENT )
        SEGMENT( 11 ) = 'HEATCONTENT'

C.........  Use the file format definition to parse the line into
C           the various data fields
        WRITE( CFIP( 1:1 ), '(I1)' ) ICC  ! country code of FIPS     
        CFIP( 2:6 ) = ADJUSTR( SEGMENT( 1 )( 1:5 ) )  ! state/county code

C.........  Replace blanks with zeros        
        DO I = 1,FIPLEN3
            IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
        END DO

        FIREID = ADJUSTL( SEGMENT( 2 ) )  ! fire ID
        LOCID  = ADJUSTL( SEGMENT( 3 ) )  ! location ID
        TSCC   = ADJUSTL( SEGMENT( 4 ) )  ! scc code
        SKID   = '               '        ! dummy stack ID
        SGID   = '               '        ! dummy segment ID

C.........  Re-define line buffer with other pollutants (CO,NOX,SO2, and others)
        IF( BKSPFLAG ) THEN
            SEGMENT( 11 ) = ADJUSTL( FIREPOL( NN ) )   ! replace with additional pol names
        END IF
        
C.........  Determine number of pollutants for this line based on CAS number
        TCAS = ADJUSTL( SEGMENT( 11 ) )
        I = FINDC( TCAS, NUNIQCAS, UNIQCAS )
        
        IF( I < 1 ) THEN
            II = INDEX1( TCAS, NINVTBL, ITNAMA )
            IF( II > 0 ) THEN
                MESG = 'FATAL ERROR: Pollutant ' // TRIM( TCAS )  //
     &             ' is not available in a list of CAS pollutants'//
     &             ' from $INVDIR/other/INVTABLE.'// CRLF()//BLANK10//
     &             ' Please update a list of pollutant names ' //
     &             '( #DATA ) in a master wildfire inventory file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            NPOLPERLN = 0
        ELSE
            NPOLPERLN = UCASNKEP( I )
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

        END SUBROUTINE RDSRCORLFR
