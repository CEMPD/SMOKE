
        SUBROUTINE RDDATAIDAAR( LINE, READDATA, READPOL, NPOLPERLN, 
     &                          IYEAR, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an IDA format area-source inventory
C      file and returns the inventory data values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen (01/03) based on rdidaar.f
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
     &                                READDATA( NPOLPERLN,NARPPOL3 )  ! array of data values
        CHARACTER(LEN=IOVLEN3), INTENT (OUT) :: READPOL( NPOLPERLN )  ! array of pollutant names
        INTEGER,                INTENT(INOUT):: NPOLPERLN             ! no. pollutants per line
        INTEGER,                INTENT (OUT) :: IYEAR                 ! inventory year
        LOGICAL,                INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,                INTENT (OUT) :: EFLAG                 ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 63  ! maximum pollutants in file
        INTEGER, PARAMETER :: AROTWIDE = 47  ! total width of all pol fields

C...........   Local parameter arrays...
C...........   Start and end positions in the file format of the first set
C              of pollutant fields.
        INTEGER, PARAMETER :: ISINIT( NARPPOL3 ) = 
     &                              ( / 16,26,36,47,54,57 / )

        INTEGER, PARAMETER :: IEINIT( NARPPOL3 ) = 
     &                              ( / 25,35,46,53,56,62 / )

C...........   Local arrays
        INTEGER         IS( NARPPOL3 )  ! start position for each pol char
        INTEGER         IE( NARPPOL3 )  ! end position for each pol char
        
C...........   Other local variables
        INTEGER         I,J     ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
 
        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDDATAIDAAR' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAIDAAR

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

C.........  Initialize start and end positions
        IS = ISINIT  ! array
        IE = IEINIT  ! array

C.........  Use the file format definition to parse the line into
C           the various data fields
        DO I = 1,NPOL
            DO J = 1,NARPPOL3

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
            IS = IS + AROTWIDE  ! array
            IE = IE + AROTWIDE  ! array            
        
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

        END SUBROUTINE RDDATAIDAAR
