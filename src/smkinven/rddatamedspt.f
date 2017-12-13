
        SUBROUTINE RDDATAMEDSPT( LINE, READDATA, READPOL, NPOLPERLN, 
     &                          IYEAR, CORS, BLID, DESC, HT, DM, TK,
     &                          FL, VL, SIC, LAT, LON, HDRFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an MEDS format point-source inventory
C      file and returns the inventory data values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created on 10/2013 by B.H. Baek
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
        USE MODINFO, ONLY: TMPNAM

        USE MODSOURC, ONLY: NMEDGRD, CMEDGRD

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)           CRLF
        INTEGER                INDEX1, STR2INT
        
        EXTERNAL   CRLF, INDEX1, STR2INT

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT  (IN) :: LINE                  ! input line
        CHARACTER(*),       INTENT (OUT) :: 
     &                            READDATA( NPOLPERLN,NPTPPOL3 )  ! array of data values
        CHARACTER(IOVLEN3), INTENT (OUT) :: READPOL( NPOLPERLN )  ! array of pollutant names
        INTEGER,            INTENT(INOUT):: NPOLPERLN             ! no. pollutants per line
        INTEGER,            INTENT (OUT) :: IYEAR                 ! inventory year
        CHARACTER(ORSLEN3), INTENT (OUT) :: CORS                  ! DOE plant ID
        CHARACTER(BLRLEN3), INTENT (OUT) :: BLID                  ! boiler ID
        CHARACTER(40),      INTENT (OUT) :: DESC                  ! plant description
        CHARACTER(16),      INTENT (OUT) :: HT                    ! stack height
        CHARACTER(16),      INTENT (OUT) :: DM                    ! stack diameter
        CHARACTER(16),      INTENT (OUT) :: TK                    ! exit temperature
        CHARACTER(16),      INTENT (OUT) :: FL                    ! flow rate
        CHARACTER(9),       INTENT (OUT) :: VL                    ! exit velocity
        CHARACTER(SICLEN3), INTENT (OUT) :: SIC                   ! SIC
        CHARACTER(16),      INTENT (OUT) :: LAT                   ! stack latitude
        CHARACTER(16),      INTENT (OUT) :: LON                   ! stack longitude
        LOGICAL,            INTENT (OUT) :: HDRFLAG               ! true: line is a header line

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 53  ! maximum pollutants in file

C...........   Local parameter arrays...
        
C...........   Other local variables
        INTEGER         I,J,N     ! counters and indices

        INTEGER         ROW, COL ! tmp grid row and col index
        INTEGER         INY      !  inventory year
        INTEGER         IOS      !  i/o status

        CHARACTER(300)      MESG                 !  message buffer
        CHARACTER(CHRLEN3)  ROWCOL           !  Row/Col & GAI lookup variables

        CHARACTER(16) :: PROGNAME = 'RDDATAMEDSPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAMEDSPT

        HDRFLAG = .FALSE.
        NPOLPERLN = 6

        READPOL( 1 ) = 'CO'
        READPOL( 2 ) = 'NOX'
        READPOL( 3 ) = 'SOX'
        READPOL( 4 ) = 'TOG'
        READPOL( 5 ) = 'PM'
        READPOL( 6 ) = 'NH3'

C.........  set year
        INY = STR2INT( LINE( 59:60 ) )
        IYEAR = 2000 + INY

C.........  Read source data
        CORS = ''  ! DOE plant ID
        BLID = ''  ! boiler ID
        DESC = ''  ! plant description

        HT   = ADJUSTL( LINE( 74:78 ) ) ! stack height
        DM   = '0.1'         ! stack diameter
        TK   = '273.15'      ! exit temperature
        FL   = '0.1'         ! flow rate
        VL   = '0.1'         ! exit velocity
        SIC  = ADJUSTL( LINE( 23:36 ) ) ! SIC
        
        COL  = STR2INT( LINE( 37:39 ) )
        ROW  = STR2INT( LINE( 40:42 ) )
        WRITE( ROWCOL,'( 2I3.3 )' ) COL, ROW
        N = INDEX1( ROWCOL, NMEDGRD, CMEDGRD( :,1 ) )

        LON  = ADJUSTL( CMEDGRD( N,2 ) )  ! grid cell longitude
        LAT  = ADJUSTL( CMEDGRD( N,3 ) )  ! grid cell latitude

C.........  Set all MEDS annual/avg to zero since all are daily/hourly inv
        READDATA = '0.0'

C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        END SUBROUTINE RDDATAMEDSPT
