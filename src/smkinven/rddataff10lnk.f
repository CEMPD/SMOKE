
        SUBROUTINE RDDATAFF10LNK( LINE, READDATA, READPOL, IYEAR,
     &                            DPID, DPLAT, DPLON, ARID, ARLAT, ARLON,
     &                            HT, DM, TK, FL, VL, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an FF10 format link-level inventory
C      file and returns the inventory data values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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
        USE MODINFO, ONLY: NEM, NDY, NEF, NCE, NRE, NC1, NC2, INV_MON

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         FINDC, STR2INT
        REAL            YR2DAY, STR2REAL
        LOGICAL         CHKREAL
        
        
        EXTERNAL   CRLF, FINDC, STR2INT, STR2REAL, YR2DAY, CHKREAL

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT  (IN) :: LINE                  ! input line
        CHARACTER(*),       INTENT (OUT) :: READDATA( 1,NPTPPOL3 )! array of data values
        CHARACTER(IOVLEN3), INTENT (OUT) :: READPOL( 1 )          ! array of pollutant names
        INTEGER,            INTENT (OUT) :: IYEAR                 ! inventory year
        CHARACTER(LNKLEN3), INTENT (OUT) :: DPID                  ! link depart loc id  
        CHARACTER(LNKLEN3), INTENT (OUT) :: DPLAT                 ! link depart latitude 
        CHARACTER(LNKLEN3), INTENT (OUT) :: DPLON                 ! link depart longitude 
        CHARACTER(LNKLEN3), INTENT (OUT) :: ARID                  ! link arrival loc id  
        CHARACTER(LNKLEN3), INTENT (OUT) :: ARLAT                 ! link arrival latitude 
        CHARACTER(LNKLEN3), INTENT (OUT) :: ARLON                 ! link arrival longitude 
        CHARACTER(4),       INTENT (OUT) :: HT                    ! stack height
        CHARACTER(6),       INTENT (OUT) :: DM                    ! stack diameter
        CHARACTER(4),       INTENT (OUT) :: TK                    ! exit temperature
        CHARACTER(10),      INTENT (OUT) :: FL                    ! flow rate
        CHARACTER(9),       INTENT (OUT) :: VL                    ! exit velocity
        LOGICAL,            INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,            INTENT (OUT) :: EFLAG                 ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 60  ! arbitrary maximum pollutants in file
        INTEGER, PARAMETER :: NSEG = 80      ! number of segments in line

C...........   Other local variables
        INTEGER         I, L, L1, LL       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file
        INTEGER      :: LYEAR   !  Leap year (366 days per year)
        INTEGER      :: MDAYS   !  days of modeling inventory month


        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine is called
        LOGICAL      :: BLKFLAG  = .TRUE.  ! true when it is blank
 
        CHARACTER(40)      TMPSEG          ! tmp segments of line
        CHARACTER(40)      SEGMENT( NSEG ) ! segments of line
        CHARACTER(CASLEN3) TCAS            ! tmp cas number
        CHARACTER(300)     MESG            ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDDATAFF10LNK' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAFF10LNK

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

C......... Return if the first line is a header line
        IF( SEGMENT( 18 ) == '' .OR. .NOT. CHKREAL( SEGMENT( 18 ) ) ) THEN
            HDRFLAG = .TRUE.
            RETURN
        END IF 

C.........  Use the file format definition to parse the line into
C           the various data fields
        DPID = ADJUSTL( SEGMENT( 6 ) )   ! departure location code
        ARID = ADJUSTL( SEGMENT( 7 ) )   ! arrival location code
        HT   = SEGMENT( 8 )              ! stack height
        DM   = SEGMENT( 9 )              ! stack diameter
        TK   = SEGMENT( 10 )             ! exit temperature
        FL   = SEGMENT( 11 )             ! flow rate
        VL   = SEGMENT( 12 )             ! exit velocity
        DPLAT= SEGMENT( 13 )             ! link departure latitude        
        DPLON= SEGMENT( 14 )             ! link departure longitude
        ARLAT= SEGMENT( 15 )             ! link arrival latitude        
        ARLON= SEGMENT( 16 )             ! link arrival longitude
        
        READPOL ( 1     ) = SEGMENT( 17 )
        READDATA( 1,NEM ) = SEGMENT( 18 )  ! annual emissions
        READDATA( 1,NDY ) = ''             ! average-day emissions

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

        END SUBROUTINE RDDATAFF10LNK
