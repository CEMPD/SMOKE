
        SUBROUTINE RDSRCEMSPT( LINE, CFIP, FCID, SKID, DVID, PRID,
     &                         NPOLPERLN, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an EMS-95 format point-source emission
C      file and returns the unique source characteristics. 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen (01/03) based on rdemspt.f
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
        USE MODLISTS, ONLY: MXIDAT, INVDNAM

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         INDEX1

        EXTERNAL    CRLF, INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*),       INTENT (IN) :: LINE      ! input line
        CHARACTER(LEN=FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(LEN=PLTLEN3), INTENT(OUT) :: FCID      ! facility ID
        CHARACTER(LEN=CHRLEN3), INTENT(OUT) :: SKID      ! stack ID
        CHARACTER(LEN=CHRLEN3), INTENT(OUT) :: DVID      ! device ID
        CHARACTER(LEN=CHRLEN3), INTENT(OUT) :: PRID      ! process ID
        INTEGER,                INTENT(OUT) :: NPOLPERLN ! no. pollutants per line
        LOGICAL,                INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,                INTENT(OUT) :: EFLAG     ! error flag

C...........   Other local variables
        INTEGER         I       ! counters and indices
        INTEGER         IOS     !  i/o status

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
 
        CHARACTER(LEN=IOVLEN3) CPOL            !  pollutant name
        CHARACTER*300          MESG            !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDSRCEMSPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDSRCEMSPT

C.........  Use the file format definition to parse the line into
C           the various data fields
        CFIP( 1:1 ) = '0'                    ! country code of FIPS
        CFIP( 2:3 ) = ADJUSTR( LINE( 1:2 ) ) ! state code
        CFIP( 4:6 ) = ADJUSTR( LINE( 3:5 ) ) ! county code

C.........  Replace blanks with zeros        
        DO I = 1,FIPLEN3
            IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
        END DO

        FCID = ADJUSTL( LINE(  6:20 ) )  ! facility ID
        SKID = ADJUSTL( LINE( 21:32 ) )  ! stack ID
        DVID = ADJUSTL( LINE( 33:44 ) )  ! device ID
        PRID = ADJUSTL( LINE( 45:56 ) )  ! process ID

C.........  Find pollutant name in master list and set number of pollutants per line
        CPOL = ADJUSTL( LINE( 57:61 ) )
        CALL UPCASE( CPOL )
        I = INDEX1( CPOL, MXIDAT, INVDNAM )
        
        IF( I < 1 ) THEN
            NPOLPERLN = 0
        ELSE
            NPOLPERLN = 1
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

        END SUBROUTINE RDSRCEMSPT
