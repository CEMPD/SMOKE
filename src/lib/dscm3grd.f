
        LOGICAL FUNCTION DSCM3GRD( GNAME, CNAME,
     &              CTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &              XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function returns true it successfully reads the properties of the
C      models-3 grid, and false otherwise.  It also returns those properties
C      using the subroutine arguments.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Copied and modified from DSCGRID, version 12/16/97
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

        IMPLICIT NONE
        
        INCLUDE 'PARMS3.EXT'               
        
C...........   ARGUMENTS and their descriptions.  All are output variables:
        
        CHARACTER*(*) GNAME	!  grid  sys name
        CHARACTER*(*) CNAME	!  coord sys name
        INTEGER       CTYPE	!  coord sys type
        REAL*8        P_ALP	!  first, second, third map
        REAL*8        P_BET	!  projection descriptive
        REAL*8        P_GAM	!  parameters
        REAL*8        XCENT	!  lon for coord-system X=0
        REAL*8        YCENT	!  lat for coord-system Y=0
        REAL*8        XORIG	!  X-coordinate origin of grid (map units)
        REAL*8        YORIG	!  Y-coordinate origin of grid
        REAL*8        XCELL	!  X-coordinate cell dimension
        REAL*8        YCELL	!  Y-coordinate cell dimension
        INTEGER       NCOLS	!  number of grid columns
        INTEGER       NROWS	!  number of grid rows
        INTEGER       NTHIK	!  BOUNDARY:  perimeter thickness (cells)
                    
C...........   EXTERNAL FUNCTIONS:
        
        INTEGER      GETEFILE
        INTEGER      TRIMLEN

        EXTERNAL     GETEFILE, TRIMLEN
        
C...........   File units and logical/physical names:
        INTEGER         IDEV    !  unit number of grid information file
        CHARACTER*16    LNAME   !  logical name for grid information file
        CHARACTER*400   PNAME   !  physical name for grid information file

C...........   Scratch local variables and their descriptions:
            
        INTEGER         I, J     !  indexes for names in lists
        INTEGER         IOS      !  I/O status return

        CHARACTER*16    CDUM     !  character dummy buffer
        CHARACTER*300   GDESC    !  grid description
        CHARACTER*300   MESG     !  message buffer

        LOGICAL      :: EFLAG = .FALSE.

        CHARACTER*16 :: PROGNAME = 'DSCM3GRD' ! program name

C***********************************************************************
C   begin body of function DSCM3GRD

C.........  Get value of Models-3 environment variable for grid and check
C           exit status

        LNAME = 'G_GRIDPATH'

C.........  Open grid file
        IDEV = GETEFILE( LNAME, .TRUE., .TRUE., PROGNAME )

C.........  Read grid information file
        READ( IDEV, 93010, END=111, ERR=222 ) 
     &        GNAME, GDESC, CNAME, CTYPE, P_ALP, P_BET, P_GAM, 
     &        XCENT, YCENT, XORIG, YORIG, CDUM, XCELL, YCELL,
     &        NCOLS, NROWS, NTHIK

C.........  Translate variables that need to be translated into I/O API

        IF ( LEN( GNAME ) .GT. 16 ) THEN
            WRITE( MESG,94010 )
     &          'Grid "', GNAME, '" Max name length 16; actual:', 
     &          LEN( GNAME )
            CALL M3WARN( 'DSCM3GRD', 0, 0, MESG )
            DSCM3GRD = .FALSE.
            RETURN
        END IF          !  if len( gname ) > 16, or if len( vname ) > 16
                
        DSCM3GRD = .TRUE.

        CLOSE( IDEV )

        RETURN

111     MESG = 'ERROR: End of file "' // LNAME( 1:TRIMLEN( LNAME ) ) //
     &         '" reached too early.'
        CALL M3MSG2( MESG )

222     DSCM3GRD = .FALSE.

        CLOSE( IDEV )

        RETURN
        
C******************  FORMAT  STATEMENTS   ******************************
                
C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( A, A, A, I8, 7E10.2, A, 2E10.2, 3I8 )
                
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )


        END

