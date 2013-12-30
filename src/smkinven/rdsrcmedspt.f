
        SUBROUTINE RDSRCMEDSPT( LINE, CFIP, FCID, PTID, SKID, SGID, TSCC,
     &                         NPOLPERLN, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an MEDS format point-source inventory
C      file and returns the unique source characteristics.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created on 12/2013 by B.H. Baek  
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
        USE MODSOURC, ONLY: NMEDGAI, COABDST
        
        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)           CRLF
        INTEGER                INDEX1
        INTEGER                STR2INT
 
        EXTERNAL   CRLF, INDEX1, STR2INT

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT (IN) :: LINE      ! input line
        CHARACTER(FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(PLTLEN3), INTENT(OUT) :: FCID      ! facility ID
        CHARACTER(CHRLEN3), INTENT(OUT) :: PTID      ! point ID
        CHARACTER(CHRLEN3), INTENT(OUT) :: SKID      ! stack ID
        CHARACTER(CHRLEN3), INTENT(OUT) :: SGID      ! segment ID
        CHARACTER(SCCLEN3), INTENT(OUT) :: TSCC      ! scc code
        INTEGER,            INTENT(OUT) :: NPOLPERLN ! no. pollutants per line
        LOGICAL,            INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,            INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 53  ! maximum pollutants in file

C...........   Other local variables
        INTEGER         I       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file
        INTEGER         ROW, COL  ! tmp row and col

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
 
        CHARACTER(300)       MESG    !  message buffer
        CHARACTER(CHRLEN3)   GAI     !  GAI lookup code

        CHARACTER(16) :: PROGNAME = 'RDSRCMEDSPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDSRCMEDSPT

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXPOLFIL, .FALSE., .TRUE., .TRUE., 
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
        GAI = ADJUSTL( LINE( 71:73 ) )  ! GAI lookup code

        IF( .NOT. ALLOCATED( COABDST ) ) THEN
           MESG='ERROR: MUST set IMPORT_MEDS_YN to Y to process'
     &        //' pregridded MEDS-formatted inventory'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        I = INDEX1( GAI, NMEDGAI, COABDST( :,1 ) )
        IF( I < 1 ) THEN
            MESG = 'ERROR: Can not find GAI code ' // GAI //
     &             ' from GAI_LOOKUP_TABLE input file'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CFIP = COABDST( I,2 )( 4:9 )     ! FIPS code
        FCID = ADJUSTL( LINE( 43:51 ) )  ! platn/facility ID
        SKID = ADJUSTL( LINE( 52:56 ) )  ! stack ID
        PTID = ADJUSTL( LINE( 37:39 ) )  ! Column ID
        SGID = ADJUSTL( LINE( 40:42 ) )  ! Row ID
        TSCC = ADJUSTL( LINE(  9:22 ) )  ! EIC code

C.........  Make sure routine knows it's been called already
        FIRSTIME = .FALSE.

C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDSRCMEDSPT
