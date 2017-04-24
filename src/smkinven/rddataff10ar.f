
        SUBROUTINE RDDATAFF10AR( LINE, READDATA, READPOL, IYEAR, 
     &                           SRCTYP, SIC, SHAPE, EXTORL, 
     &                           HDRFLAG, AVEFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an ORL format area-source inventory
C      file and returns the inventory data values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by B.H. Baek (Aug, 2011)
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
        USE MODINFO, ONLY: NEM, NDY, NEF, NCE, NRE, NRP, INV_MON

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: FF10INVFLAG 

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         FINDC, STR2INT
        REAL            YR2DAY, STR2REAL
        LOGICAL         CHKINT

        EXTERNAL    CRLF, FINDC, STR2INT, STR2REAL, CHKINT, YR2DAY

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT  (IN) :: LINE                  ! input line
        CHARACTER(*),       INTENT (OUT) :: READDATA( 1,NARPPOL3 )! array of data values
        CHARACTER(IOVLEN3), INTENT (OUT) :: READPOL( 1 )          ! pollutant name
        INTEGER,            INTENT (OUT) :: IYEAR                 ! inventory year
        CHARACTER(STPLEN3), INTENT (OUT) :: SRCTYP                ! source type code
        CHARACTER(SICLEN3), INTENT (OUT) :: SIC                   ! SIC
        CHARACTER(SHPLEN3), INTENT (OUT) :: SHAPE                 ! SHAPE_ID
        CHARACTER(EXTLEN3), INTENT (OUT) :: EXTORL                ! additional ext vars
        LOGICAL,            INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,            INTENT (OUT) :: AVEFLAG               ! true: Aveday inv is processed
        LOGICAL,            INTENT (OUT) :: EFLAG                 ! error flag
        
C...........   Local parameters
        INTEGER, PARAMETER :: MXDATFIL = 60  ! arbitrary max no. data variables
        INTEGER, PARAMETER :: NSEG = 60      ! number of segments in line

C...........   Other local variables
        INTEGER         I        ! counters and indices

        INTEGER, SAVE:: ICC      !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY      !  inventory year
        INTEGER         IOS      !  i/o status
        INTEGER, SAVE:: NPOA     !  number of pollutants in file
        INTEGER      :: MDAYS    !  days of modeling inventory month
        INTEGER      :: LYEAR    !  Leap year (366 days per year)
        REAL         :: AVEINV   !  annual total estimate from monthly total VMT

        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine is called
        LOGICAL      :: BLKFLAG  = .TRUE.  ! true when it is blank
 
        CHARACTER(25)      SEGMENT( NSEG ) ! segments of line
        CHARACTER(25)      TMPSEG          ! tmp segments of line
        CHARACTER(CASLEN3) TCAS            ! tmp cas number
        CHARACTER(300)     MESG            !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDDATAFF10AR' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAFF10AR

C.........  Scan for header lines and check to ensure all are set 
C           properly (country and year required)
        CALL GETHDR( MXDATFIL, .TRUE., .TRUE., .FALSE., 
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

C.........  Separate line into segments
        CALL PARSLINE( LINE, NSEG, SEGMENT )

C......... Return if the first line is a header line
        IF( .NOT. CHKINT( SEGMENT( 2 ) ) ) THEN
            HDRFLAG = .TRUE.
            RETURN
        END IF 

C.........  Use the file format definition to parse the line into
C           the various data fields
        READPOL ( 1     ) = ADJUSTL( SEGMENT( 8 ) )
        READDATA( 1,NEM ) = SEGMENT( 9 )
        READDATA( 1,NDY ) = ''
        READDATA( 1,NEF ) = '-9'
        READDATA( 1,NCE ) = '-9'
        READDATA( 1,NRE ) = '-9'
        READDATA( 1,NRP ) = '-9'

        SRCTYP = ' '   ! source type code
        SIC    = ' '
        SHAPE  = ADJUSTL( SEGMENT(  5 ) )
        EXTORL = ' '   ! extended orl (N/A)

C.........  Compute annual total based on monthly total
        AVEINV = 0.0
        AVEFLAG = .FALSE.
        IF( INV_MON > 0 ) THEN

            DO I = 1, 12
                IF( LEN_TRIM( SEGMENT( 20+I ) ) < 1 ) THEN
                    SEGMENT( 20+I ) = '0.0'
                END IF
                AVEINV = AVEINV + STR2REAL( SEGMENT( 20+I ) )
            END DO

            IF( AVEINV > 0.0 ) AVEFLAG = .TRUE.

            IF( AVEFLAG ) THEN

                READDATA( 1,NEM ) = '0.0'
                READDATA( 1,NDY ) = SEGMENT( 20 + INV_MON )

                MDAYS = MON_DAYS( INV_MON )    ! day of months

                LYEAR = INT( 1 / YR2DAY ( INY ) )                ! convert year to days
                IF( LYEAR > 365 .AND. INV_MON == 2 ) MDAYS = 29  ! leap year (feb = 29days)

                AVEINV = STR2REAL( READDATA(1,NDY) ) / MDAYS
                WRITE( READDATA( 1,NDY ), '( E15.10 )' ) AVEINV

                IF( AVEINV < 0.0 ) THEN
                    MESG = 'ERROR: Can not process negative value'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

             END IF

        END IF

C.........  Reset annual total inventory to zero for daily/hourly FF10 processing
        IF( FF10INVFLAG ) THEN
            READPOL ( 1     ) = SEGMENT( 9 )
            READDATA( 1,NEM ) = ''
            READDATA( 1,NDY ) = ''
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

        END SUBROUTINE RDDATAFF10AR
