
        SUBROUTINE RDDATAIDAPT( LINE, READDATA, READPOL, NPOLPERLN, 
     &                          IYEAR, CORS, BLID, DESC, HT, DM, TK,
     &                          FL, VL, SIC, LAT, LON, HDRFLAG, 
     &                          EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an IDA format point-source inventory
C      file and returns the inventory data values.
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
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: TMPNAM

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        
        EXTERNAL   CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*),       INTENT  (IN) :: LINE                  ! input line
        CHARACTER(LEN=*),       INTENT (OUT) :: 
     &                                READDATA( NPOLPERLN,NPTPPOL3 )  ! array of data values
        CHARACTER(LEN=IOVLEN3), INTENT (OUT) :: READPOL( NPOLPERLN )  ! array of pollutant names
        INTEGER,                INTENT(INOUT):: NPOLPERLN             ! no. pollutants per line
        INTEGER,                INTENT (OUT) :: IYEAR                 ! inventory year
        CHARACTER(LEN=ORSLEN3), INTENT (OUT) :: CORS                  ! DOE plant ID
        CHARACTER(LEN=6),       INTENT (OUT) :: BLID                  ! boiler ID
        CHARACTER(LEN=40),      INTENT (OUT) :: DESC                  ! plant description
        CHARACTER(LEN=4),       INTENT (OUT) :: HT                    ! stack height
        CHARACTER(LEN=6),       INTENT (OUT) :: DM                    ! stack diameter
        CHARACTER(LEN=4),       INTENT (OUT) :: TK                    ! exit temperature
        CHARACTER(LEN=10),      INTENT (OUT) :: FL                    ! flow rate
        CHARACTER(LEN=9),       INTENT (OUT) :: VL                    ! exit velocity
        CHARACTER(LEN=SICLEN3), INTENT (OUT) :: SIC                   ! SIC
        CHARACTER(LEN=9),       INTENT (OUT) :: LAT                   ! stack latitude
        CHARACTER(LEN=9),       INTENT (OUT) :: LON                   ! stack longitude
        LOGICAL,                INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,                INTENT (OUT) :: EFLAG                 ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 53  ! maximum pollutants in file
        INTEGER, PARAMETER :: PTOTWIDE = 52  ! total width of all pol fields

C...........   Local parameter arrays...
C...........   Start and end positions in the file format of the first set
C              of pollutant fields.
        INTEGER, PARAMETER :: ISINIT( NPTPPOL3 ) = 
     &                              ( / 250,263,276,283,286,296,299 / )

        INTEGER, PARAMETER :: IEINIT( NPTPPOL3 ) = 
     &                              ( / 262,275,282,285,295,298,301 / )

C...........   Local arrays
        INTEGER         IS( NPTPPOL3 )  ! start position for each pol char
        INTEGER         IE( NPTPPOL3 )  ! end position for each pol char
        
C...........   Other local variables
        INTEGER         I,J     ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
 
        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDDATAIDAPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAIDAPT

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
            NPOLPERLN = NPOL
            IYEAR = INY
            RETURN
        ELSE
            HDRFLAG = .FALSE.
        END IF

C.........  Set pollutants for this line
        READPOL = TMPNAM

C.........  Read source data
        CORS = ADJUSTL( LINE( 48: 53 ) )  ! DOE plant ID
        BLID = ADJUSTL( LINE( 54: 59 ) )  ! boiler ID
        DESC = ADJUSTL( LINE( 62:101 ) )  ! plant description
        
        HT   = LINE( 120:123 )  ! stack height
        DM   = LINE( 124:129 )  ! stack diameter
        TK   = LINE( 130:133 )  ! exit temperature
        FL   = LINE( 134:143 )  ! flow rate
        VL   = LINE( 144:152 )  ! exit velocity
        SIC  = LINE( 227:230 )  ! SIC
        LAT  = LINE( 231:239 )  ! stack latitude
        LON  = LINE( 240:248 )  ! stack longitude

C.........  Initialize start and end positions
        IS = ISINIT  ! array
        IE = IEINIT  ! array

C.........  Use the file format definition to parse the line into
C           the various data fields
        DO I = 1,NPOL
            DO J = 1,NPTPPOL3

C.................  Make sure data array is large enough for value
                IF( IE(J)-IS(J)+1 > LEN( READDATA( 1,1 ) ) ) THEN
                    WRITE( MESG,94010 ) 'INTERNAL ERROR: Length of ' //
     &                 'data array ', LEN( READDATA( 1,1 ) ), 
     &                 ' not sufficient for field size of ', 
     &                 IE(J)-IS(J)+1
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
                READDATA( I,J ) = LINE( IS( J ):IE( J ) )
            END DO
            
C.............  Update start and end positions
            IS = IS + PTOTWIDE  ! array
            IE = IE + PTOTWIDE  ! array            
        
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

        END SUBROUTINE RDDATAIDAPT
