
        SUBROUTINE RDDATAIDAMB( LINE, READDATA, READPOL, NPOLPERLN, 
     &                          IYEAR, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an IDA format mobile-source inventory
C      file and returns the inventory data values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen (01/03) based on rddataidaar.f
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
        USE MODINFO, ONLY: TMPNAM, NPPOL, NPACT

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)           CRLF
        
        EXTERNAL   CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT  (IN) :: LINE                  ! input line
        CHARACTER(*),       INTENT (OUT) :: 
     &                                READDATA( NPOLPERLN,NPPOL )     ! array of data values
        CHARACTER(IOVLEN3), INTENT (OUT) :: READPOL( NPOLPERLN )  ! array of pollutant names
        INTEGER,            INTENT(INOUT):: NPOLPERLN             ! no. pollutants per line
        INTEGER,            INTENT (OUT) :: IYEAR                 ! inventory year
        LOGICAL,            INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,            INTENT (OUT) :: EFLAG                 ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 112  ! maximum pollutants in file
        INTEGER, PARAMETER :: MBOTWIDE = 20   ! total width of all pol fields

C...........   Local parameter arrays...
C...........   Start and end positions in the file format of the first set
C              of pollutant fields.
        INTEGER  :: ISINIT( NMBPPOL3 ) = ( / 26,36 / )

        INTEGER  :: IEINIT( NMBPPOL3 ) = ( / 35,45 / )

C...........   Local arrays
        INTEGER         IS( NPPOL )  ! start position for each pol char
        INTEGER         IE( NPPOL )  ! end position for each pol char
        
C...........   Local allocatable arrays
        CHARACTER(25), ALLOCATABLE :: SEGMENT( : )  ! list-formatted strings
        
C...........   Other local variables
        INTEGER         I,J,K   ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file
        INTEGER         NSEG    ! number of input segments

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
        LOGICAL, SAVE:: FIXED    = .TRUE. ! true: input file is fixed-format
 
        CHARACTER(5)   TMPBUF  ! temporary string buffer
        CHARACTER(300) MESG    !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDDATAIDAMB' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAIDAMB

C.........  Scan for header lines and check to ensure all are set 
C           properly
C.........  Skip if line is units to avoid problems in gethdr.f
        TMPBUF = LINE( 2:6 )
        CALL UPCASE( TMPBUF )
        IF( TMPBUF /= 'UNITS' ) THEN
            CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .TRUE., 
     &                   LINE, ICC, INY, NPOL, IOS )
        ELSE
            IOS = 0
        END IF

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

C.........  If a header line was encountered, check for free or fixed format
        IF( IOS >= 0 ) THEN
        
            IF( LINE( 2:5 ) == 'TYPE' ) THEN

C.................  Try to find activity in name of data, otherwise, assume
C                   emissions
                MESG = LINE
                CALL UPCASE( MESG )
            
                IF( INDEX( MESG, 'ACTIVITY' ) > 0 ) THEN
                    FIXED = .FALSE.
                ELSE
                    FIXED = .TRUE.
                END IF
            END IF
            
            HDRFLAG = .TRUE.
            NPOLPERLN = NPOL
            IYEAR = INY
            RETURN
        ELSE
            HDRFLAG = .FALSE.
        END IF

C.........  Set pollutants for this line
        READPOL = TMPNAM

C.........  If not fixed format, allocate memory for number of segments
        IF( .NOT. FIXED .AND. .NOT. ALLOCATED( SEGMENT ) ) THEN
            NSEG = 4 + NPOL * NPACT
            ALLOCATE( SEGMENT( NSEG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )
            SEGMENT = ' '   ! array
        END IF

C.........  Initialize start and end positions
        IS = ISINIT  ! array
        IE = IEINIT  ! array

C.........  Use the file format definition to parse the line into
C           the various data fields
        IF( FIXED ) THEN
            DO I = 1,NPOL
                DO J = 1,NPPOL

C.....................  Make sure data array is large enough for value
                    IF( IE(J)-IS(J)+1 > LEN( READDATA( 1,1 ) ) ) THEN
                        WRITE( MESG,94010 ) 'INTERNAL ERROR: Length ' //
     &                     'of data array ', LEN( READDATA( 1,1 ) ), 
     &                     ' not sufficient for field size of ', 
     &                     IE(J)-IS(J)+1
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                
                    READDATA( I,J ) = LINE( IS( J ):IE( J ) )
                END DO
            
C.................  Update start and end positions
                IS = IS + MBOTWIDE  ! array
                IE = IE + MBOTWIDE  ! array            
        
            END DO
        ELSE
            CALL PARSLINE( LINE, NSEG, SEGMENT )
            
            K = 4
            DO I = 1,NPOL
                DO J = 1,NPACT
                    K = K + 1
                    READDATA( I,J ) = SEGMENT( K )
                END DO
            END DO
            
C.............  Since there are 2 values per pollutant but only 1 per activity,
C               need to fill in second spot to avoid errors later on
            READDATA( :,2 ) = ' '
                
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

        END SUBROUTINE RDDATAIDAMB
