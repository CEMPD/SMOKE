
        SUBROUTINE RDDATAEMSMB( LINE, READDATA, READPOL, 
     &                          NVARPERLN, IYEAR, X1, Y1, X2, Y2,
     &                          ZONE, LNKFLAG, CFIP, CROAD, CLNK,
     &                          HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an EMS format mobile-source inventory
C      file and returns the inventory data values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by C. Seppanen (01/03) based on rdidamb.f
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
        USE MODINFO, ONLY: TMPNAM, NPPOL, NEM
        
        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        LOGICAL                ENVYN
        
        EXTERNAL   CRLF, ENVYN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*),       INTENT (IN) :: LINE      ! input line
        CHARACTER(LEN=*),       INTENT (OUT) :: 
     &                                READDATA( NVARPERLN,NPPOL )     ! array of data values
        CHARACTER(LEN=IOVLEN3), INTENT (OUT) :: READPOL( NVARPERLN )  ! array of pollutant names
        INTEGER,                INTENT(INOUT):: NVARPERLN ! no. variables per line
        INTEGER,                INTENT(OUT) :: IYEAR     ! inventory year
        CHARACTER(LEN=25),      INTENT(OUT) :: X1        ! x-dir link coord 1
        CHARACTER(LEN=25),      INTENT(OUT) :: Y1        ! y-dir link coord 1
        CHARACTER(LEN=25),      INTENT(OUT) :: X2        ! x-dir link coord 2
        CHARACTER(LEN=25),      INTENT(OUT) :: Y2        ! y-dir link coord 2
        CHARACTER(LEN=2),       INTENT(OUT) :: ZONE      ! time zone
        LOGICAL,                INTENT(OUT) :: LNKFLAG   ! true: line contains link information
        CHARACTER(LEN=FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
        CHARACTER(LEN=RWTLEN3), INTENT(OUT) :: CROAD     ! roadway type
        CHARACTER(LEN=LNKLEN3), INTENT(OUT) :: CLNK      ! link ID
        LOGICAL,                INTENT(OUT) :: HDRFLAG   ! true: line is a header line
        LOGICAL,                INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXDATFIL = 60 ! arbitrary max data variables in file

C...........   Local allocatable arrays
        CHARACTER(LEN=25), ALLOCATABLE :: SEGMENT( : )  ! list-formatted strings
        
C...........   Other local variables
        INTEGER         I       ! counters and indices

        INTEGER, SAVE:: FMTCASE !  code for format case
        INTEGER, SAVE:: ICC = 0 !  position of CNTRY in CTRYNAM
        INTEGER         INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPRECOL !  no. src char columns for list-directed
        INTEGER, SAVE:: NVAR    !  number of variables in file
        INTEGER         NSEG    ! number of input segments

        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine is called
        LOGICAL, SAVE:: FIXED    = .FALSE. ! true: input file is fixed-format
        LOGICAL, SAVE:: LFLAG    = .FALSE. ! true: link file
 
        CHARACTER(LEN=300)     MESG            ! message buffer

        CHARACTER*16 :: PROGNAME = 'RDDATAEMSMB' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAEMSMB

        IF( FIRSTIME ) THEN
            MESG = 'Indicator for fixed-column EMS-95 format'
            FIXED = ENVYN( 'SMK_EMS95_FIXFMT', MESG, FIXED, IOS )
            
            FMTCASE = 0
            
            FIRSTIME = .FALSE.
        END IF

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXDATFIL, .FALSE., .FALSE., .TRUE., 
     &               LINE, ICC, INY, NVAR, IOS )

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

C.........  If a header line was encountered, check if file is non-link or link
        IF( IOS >= 0 ) THEN
            
            IF( LINE( 2:LEN_TRIM( LINE ) ) == 'LINK' ) THEN
                LFLAG = .TRUE.
                
                IF( FIXED ) THEN
                    FMTCASE = 3      ! link fixed
                ELSE
                    FMTCASE = 4      ! link list
                    NPRECOL = 8
                END IF
            ELSE IF( LINE( 2:LEN_TRIM( LINE ) ) == 'NONLINK' ) THEN
                LFLAG = .FALSE.
                
                IF( FIXED ) THEN
                    FMTCASE = 1       ! non-link fixed
                ELSE
                    FMTCASE = 2       ! non-link list
                    NPRECOL = 2
                END IF
            END IF
            
            HDRFLAG = .TRUE.
            NVARPERLN = NVAR
            IYEAR = INY
            RETURN
        ELSE
            HDRFLAG = .FALSE.
        END IF

C.........  Make sure file format has been set
        IF( FMTCASE == 0 ) THEN
            MESG = 'INTERNAL ERROR: Mobile EMS-95 file format ' //
     &             'is not recognized'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set pollutants for this line
        READPOL = TMPNAM
        
C.........  If not fixed format, allocate memory for number of segments
        IF( .NOT. FIXED .AND. .NOT. ALLOCATED( SEGMENT ) ) THEN
            NSEG = NPRECOL + NVAR
            ALLOCATE( SEGMENT( NSEG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )
            SEGMENT = ' '   ! array
        END IF
        
C.........  Use the file format definition to parse the LINE into
C           the various data fields        
        X1 = ' '
        Y1 = ' '
        X2 = ' '
        Y2 = ' '
        ZONE = ' '
        
        IF( FIXED ) THEN

C.............  Store source information to match with VMTMIX file
            WRITE( CFIP( 1:1 ), '(I1)' ) ICC  ! country code of FIPS
            CFIP ( 2:3 ) = ADJUSTR( LINE( 1:2 ) ) ! state
            CFIP ( 4:6 ) = ADJUSTR( LINE( 3:5 ) ) ! county
            CROAD( 1:1 ) = LINE( 6:6 )  ! area type
            CROAD( 2:3 ) = ADJUSTL( LINE( 7:10 ) ) ! facility type  <-- PROBLEM!!!!
            CLNK         = ' '          ! link

            DO I = 0,NVAR-1
                READDATA( I,NEM ) = LINE( 11+(I*8):18+(I*8) )
            END DO
        ELSE
            CALL PARSLINE( LINE, NSEG, SEGMENT )
            
            CFIP  = ADJUSTR( SEGMENT( 1 )( 1:FIPLEN3 ) )  ! fips code
            CROAD = ADJUSTL( SEGMENT( 2 ) )               ! roadway type
            
            IF( LFLAG ) THEN
                CLNK = ADJUSTL( SEGMENT( 3 ) )   ! link ID
                X1   = SEGMENT( 4 )
                Y1   = SEGMENT( 5 )
                X2   = SEGMENT( 6 )
                Y2   = SEGMENT( 7 )
                ZONE = SEGMENT( 8 )
            ELSE
                CLNK = ' '
            END IF
            
            DO I = 1, NVAR
                READDATA( I,NEM ) = SEGMENT( NPRECOL + I )
            END DO
        END IF

C.........  Replace blanks with zeros        
        DO I = 1,FIPLEN3
            IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
        END DO

C.........  Since there are 2 values per pollutant but only 1 per activity,
C           need to fill in second spot to avoid errors later on
        READDATA( :,2 ) = ' '

        LNKFLAG = LFLAG

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

        END SUBROUTINE RDDATAEMSMB
