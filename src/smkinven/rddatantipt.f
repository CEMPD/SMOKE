
        SUBROUTINE RDDATANTIPT( LINE, READDATA, READPOL, IYEAR, 
     &                          DESC, HT, DM, TK, FL, VL, SIC, 
     &                          MACT, NAICS, LAT, LON, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an NTI format point-source inventory
C      file and returns the inventory data values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen (01/03) based on rddataidapt.f
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
        USE MODINFO, ONLY: NEM, NOZ, NEF, NCE, NRE, NC1, NC2

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        INTEGER                FINDC
        
        EXTERNAL   CRLF, FINDC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*),       INTENT  (IN) :: LINE                  ! input line
        CHARACTER(LEN=*),       INTENT (OUT) :: READDATA( 1,NPTPPOL3 )! array of data values
        CHARACTER(LEN=IOVLEN3), INTENT (OUT) :: READPOL( 1 )          ! array of pollutant names
        INTEGER,                INTENT (OUT) :: IYEAR                 ! inventory year
        CHARACTER(LEN=40),      INTENT (OUT) :: DESC                  ! plant description
        CHARACTER(LEN=4),       INTENT (OUT) :: HT                    ! stack height
        CHARACTER(LEN=6),       INTENT (OUT) :: DM                    ! stack diameter
        CHARACTER(LEN=4),       INTENT (OUT) :: TK                    ! exit temperature
        CHARACTER(LEN=10),      INTENT (OUT) :: FL                    ! flow rate
        CHARACTER(LEN=9),       INTENT (OUT) :: VL                    ! exit velocity
        CHARACTER(LEN=SICLEN3), INTENT (OUT) :: SIC                   ! SIC
        CHARACTER(LEN=MACLEN3), INTENT (OUT) :: MACT                  ! MACT code
        CHARACTER(LEN=NAILEN3), INTENT (OUT) :: NAICS                 ! NAICS code
        CHARACTER(LEN=9),       INTENT (OUT) :: LAT                   ! stack latitude
        CHARACTER(LEN=9),       INTENT (OUT) :: LON                   ! stack longitude
        LOGICAL,                INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,                INTENT (OUT) :: EFLAG                 ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 60  ! arbitrary maximum pollutants in file
        INTEGER, PARAMETER :: NSEG = 25      ! number of segments in line

C...........   Other local variables
        INTEGER         I       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
 
        CHARACTER(LEN=25)      SEGMENT( NSEG ) ! segments of line
        CHARACTER(LEN=CASLEN3) TCAS            ! tmp cas number
        CHARACTER(LEN=300)     MESG            ! message buffer

        CHARACTER*16 :: PROGNAME = 'RDDATANTIPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATANTIPT

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
            IYEAR = INY
            RETURN
        ELSE
            HDRFLAG = .FALSE.
        END IF

C.........  Separate line into segments
        CALL PARSLINE( LINE, NSEG, SEGMENT )

C.........  Use the file format definition to parse the line into
C           the various data fields
        DESC  = ADJUSTL( SEGMENT( 7 ) )   ! plant description        
        HT    = SEGMENT( 9 )              ! stack height
        DM    = SEGMENT( 10 )             ! stack diameter
        TK    = SEGMENT( 11 )             ! exit temperature
        FL    = SEGMENT( 12 )             ! flow rate
        VL    = SEGMENT( 13 )             ! exit velocity
        SIC   = SEGMENT( 14 )             ! SIC
        MACT  = ADJUSTL( SEGMENT( 15 ) )  ! MACT code
        NAICS = ADJUSTL( SEGMENT( 16 ) )  ! NAICS code
        LAT   = SEGMENT( 17 )             ! stack latitude
        LON   = SEGMENT( 18 )             ! stack longitude

        READPOL ( 1     ) = SEGMENT( 19 )
        READDATA( 1,NEM ) = SEGMENT( 20 ) ! annual emissions
        READDATA( 1,NOZ ) = SEGMENT( 21 ) ! average-day emissions
        READDATA( 1,NEF ) = ' '           ! emission factor
        READDATA( 1,NCE ) = SEGMENT( 22 ) ! control efficiency
        READDATA( 1,NRE ) = SEGMENT( 23 ) ! rule effectiveness
        READDATA( 1,NC1 ) = SEGMENT( 24 ) ! primary control equipment code
        READDATA( 1,NC2 ) = SEGMENT( 25 ) ! secondary control equipment code
            
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

        END SUBROUTINE RDDATANTIPT
