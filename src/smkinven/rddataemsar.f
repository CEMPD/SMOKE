
        SUBROUTINE RDDATAEMSAR( LINE, READDATA, READPOL, NPOLPERLN, 
     &                          IYEAR, TIMEPERIOD, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an EMS format area-source inventory
C      file and returns the inventory data values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen (01/03) based on rdntiar.f
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
        USE MODINFO, ONLY: NEM, NOZ, NEF, NCE, NRE, NRP
        
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
        CHARACTER(LEN=*),       INTENT  (IN) :: LINE                  ! input line
        CHARACTER(LEN=*),       INTENT (OUT) :: READDATA( 1,NARPPOL3 )! array of data values
        CHARACTER(LEN=IOVLEN3), INTENT (OUT) :: READPOL( 1 )          ! pollutant name
        INTEGER,                INTENT (OUT) :: NPOLPERLN             ! no. pollutants per line
        INTEGER,                INTENT (OUT) :: IYEAR                 ! inventory year
        CHARACTER(LEN=*),       INTENT (OUT) :: TIMEPERIOD            ! time period type
        LOGICAL,                INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,                INTENT (OUT) :: EFLAG                 ! error flag
        
C...........   Local parameters
        INTEGER, PARAMETER :: MXDATFIL = 60  ! arbitrary max no. data variables

C...........   Other local variables
        INTEGER         I        ! counters and indices

        INTEGER, SAVE:: ICC      !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY      !  inventory year
        INTEGER         IOS      !  i/o status
        INTEGER, SAVE:: NPOA     !  number of pollutants in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
 
        CHARACTER(LEN=300)     MESG            !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDDATAEMSAR' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAEMSAR

C.........  Scan for header lines and check to ensure all are set 
C           properly (country and year required)
        CALL GETHDR( MXDATFIL, .FALSE., .FALSE., .FALSE., 
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
            HDRFLAG = .TRUE.
            IYEAR = INY
            RETURN
        ELSE
            HDRFLAG = .FALSE.
        END IF

C.........  Use the file format definition to parse the line into
C           the various data fields
        READPOL ( 1 ) = ADJUSTL( LINE( 21:25 ) )
        CALL UPCASE( READPOL( 1 ) )
         
        READDATA( 1,NEM ) = LINE( 52:65 )
        READDATA( 1,NOZ ) = '0'
        READDATA( 1,NEF ) = '0'
        READDATA( 1,NCE ) = LINE( 88:94 )
        READDATA( 1,NRE ) = '100.0'
        READDATA( 1,NRP ) = '100.0'
        
        TIMEPERIOD = LINE( 95:96 )
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

        END SUBROUTINE RDDATAEMSAR
