
        SUBROUTINE RDDATAEMSPT( LINE, READDATA, READPOL, 
     &                          NPOLPERLN, TIMEPERIOD, 
     &                          HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an EMS format point-source emission
C      file and returns the inventory data values.
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
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NPPOL, NEM, NDY, NEF, NCE, NRE, NC1, NC2
        
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: MXIDAT, INVDNAM
        
        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         INDEX1

        EXTERNAL    INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*),       INTENT (IN) :: LINE       ! input line
        CHARACTER(LEN=*),       INTENT(OUT) :: READDATA( 1,NPPOL )  ! array of data values
        CHARACTER(LEN=IOVLEN3), INTENT(OUT) :: READPOL( 1 )         ! pollutant name
        INTEGER,                INTENT(OUT) :: NPOLPERLN  ! no. pollutants per line
        CHARACTER(LEN=*),       INTENT(OUT) :: TIMEPERIOD ! time period type
        LOGICAL,                INTENT(OUT) :: HDRFLAG    ! true: line is a header line
        LOGICAL,                INTENT(OUT) :: EFLAG      ! error flag
        
C...........   Other local variables
        INTEGER         I       ! counters and indices
        INTEGER         IOS     !  i/o status

        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine is called
 
        CHARACTER(LEN=IOVLEN3) CPOL            !  pollutant name
        CHARACTER(LEN=300)     MESG            ! message buffer

        CHARACTER*16 :: PROGNAME = 'RDDATAEMSPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAEMSPT

C.........  Use the file format definition to parse the line into
C           the various data fields
        READPOL ( 1 ) = ADJUSTL( LINE( 57:61 ) )
        CALL UPCASE( READPOL( 1 ) )

        READDATA( 1,NEM ) = LINE(  88:100 )
        READDATA( 1,NDY ) = '0'
        READDATA( 1,NEF ) = '0'
        READDATA( 1,NCE ) = LINE( 126:132 )
        READDATA( 1,NRE ) = '100.0'
        READDATA( 1,NC1 ) = '0'
        READDATA( 1,NC2 ) = '0'
        
        TIMEPERIOD = LINE( 114:115 )
        CALL UPCASE( TIMEPERIOD )

C.........  Find pollutant name in master list and set number of pollutants per line
        I = INDEX1( READPOL( 1 ), MXIDAT, INVDNAM )
        
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

        END SUBROUTINE RDDATAEMSPT
