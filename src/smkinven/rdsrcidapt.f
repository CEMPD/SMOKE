
        SUBROUTINE RDSRCIDAPT( LINE, CFIP, FCID, PTID, SKID, SGID, TSCC,
     &                         NPOLPERLN, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an IDA format point-source inventory
C      file and returns the unique source characteristics.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen (01/03) based on rdidapt.f
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

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        
        EXTERNAL   CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*),       INTENT (IN) :: LINE      ! input line
        CHARACTER(LEN=FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(LEN=PLTLEN3), INTENT(OUT) :: FCID      ! facility ID
        CHARACTER(LEN=CHRLEN3), INTENT(OUT) :: PTID      ! point ID
        CHARACTER(LEN=CHRLEN3), INTENT(OUT) :: SKID      ! stack ID
        CHARACTER(LEN=CHRLEN3), INTENT(OUT) :: SGID      ! segment ID
        CHARACTER(LEN=SCCLEN3), INTENT(OUT) :: TSCC      ! scc code
        INTEGER,                INTENT(OUT) :: NPOLPERLN ! no. pollutants per line
        LOGICAL,                INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,                INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 53  ! maximum pollutants in file

C...........   Other local variables
        INTEGER         I       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
 
        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDSRCIDAPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDSRCIDAPT

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .TRUE., 
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
            NPOLPERLN = NPOL
        END IF

C.........  Use the file format definition to parse the line into
C           the various data fields
        WRITE( CFIP( 1:1 ), '(I1)' ) ICC  ! country code of FIPS
        CFIP( 2:3 ) = LINE( 1:2 )         ! state code
        CFIP( 4:6 ) = LINE( 3:5 )         ! county code

C.........  Replace blanks with zeros        
        DO I = 2,6
            IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
        END DO

        FCID = ADJUSTL( LINE(   6:20  ) ) ! facility ID
        PTID = ADJUSTL( LINE(  21:35  ) ) ! point ID
        SKID = ADJUSTL( LINE(  36:47  ) ) ! stack ID
        SGID = ADJUSTL( LINE(  60:61  ) ) ! segment ID
        TSCC = ADJUSTL( LINE( 102:111 ) ) ! SCC code
            
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

        END SUBROUTINE RDSRCIDAPT
