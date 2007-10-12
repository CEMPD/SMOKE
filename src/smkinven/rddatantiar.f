
        SUBROUTINE RDDATAORLAR( LINE, READDATA, READPOL, IYEAR, 
     &                          SRCTYP, EXTORL, HDRFLAG, EFLAG )

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
C      Created by C. Seppanen (01/03) based on rdntiar.f
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
        USE MODINFO, ONLY: NEM, NDY, NEF, NCE, NRE, NRP
        
        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         FINDC

        EXTERNAL    CRLF, FINDC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT  (IN) :: LINE                  ! input line
        CHARACTER(*),       INTENT (OUT) :: READDATA( 1,NARPPOL3 )! array of data values
        CHARACTER(IOVLEN3), INTENT (OUT) :: READPOL( 1 )          ! pollutant name
        INTEGER,            INTENT (OUT) :: IYEAR                 ! inventory year
        CHARACTER(STPLEN3), INTENT (OUT) :: SRCTYP                ! source type code
        CHARACTER(EXTLEN3), INTENT (OUT) :: EXTORL                ! additional ext vars
        LOGICAL,            INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,            INTENT (OUT) :: EFLAG                 ! error flag
        
C...........   Local parameters
        INTEGER, PARAMETER :: MXDATFIL = 60  ! arbitrary max no. data variables
        INTEGER, PARAMETER :: NSEG = 50      ! number of segments in line

C...........   Other local variables
        INTEGER         I        ! counters and indices

        INTEGER, SAVE:: ICC      !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY      !  inventory year
        INTEGER         IOS      !  i/o status
        INTEGER, SAVE:: NPOA     !  number of pollutants in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine is called
        LOGICAL      :: BLKFLAG  = .TRUE.  ! true when it is blank
 
        CHARACTER(25)      SEGMENT( NSEG ) ! segments of line
        CHARACTER(25)      TMPSEG          ! tmp segments of line
        CHARACTER(CASLEN3) TCAS            ! tmp cas number
        CHARACTER(300)     MESG            !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDDATAORLAR' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAORLAR

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

C.........  Use the file format definition to parse the line into
C           the various data fields
        READPOL ( 1     ) = SEGMENT( 3 )
        READDATA( 1,NEM ) = SEGMENT( 4 )
        READDATA( 1,NDY ) = SEGMENT( 5 )
        READDATA( 1,NEF ) = ' '
        READDATA( 1,NCE ) = SEGMENT( 6 )
        READDATA( 1,NRE ) = SEGMENT( 7 )
        READDATA( 1,NRP ) = SEGMENT( 8 )

        SRCTYP = ADJUSTL( SEGMENT( 9 ) )   ! source type code

C.........  Read extended orl variables and store it as string
        EXTORL = ' '
        DO I = 10, 30
            IF( SEGMENT( I ) == ' ' ) THEN
                TMPSEG = ','
            ELSE
                TMPSEG = ',' // TRIM( SEGMENT( I ) )
                BLKFLAG = .FALSE.
            ENDIF
            
            EXTORL = TRIM( EXTORL ) // TRIM( TMPSEG )
        END DO
        
        IF( BLKFLAG ) EXTORL = ' '

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

        END SUBROUTINE RDDATAORLAR
