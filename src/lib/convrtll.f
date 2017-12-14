
        SUBROUTINE CONVRTLL( NSRC, CTYPE, GDNAM, P_ALP, P_BET, 
     &                       P_GAM, XCENT, YCENT, XVALS, YVALS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine converts coordinates from the grid
C      defined by the subroutine arguments to lat-lon.
C
C  PRECONDITIONS REQUIRED:
C      
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created based on CONVRTXY 07/13 by C. Seppanen
!
!       Version 10/2016 by C. Coats:  USE M3UTILIO
C
C****************************************************************************/
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
        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C.........  SUBROUTINE ARGUMENTS
        INTEGER,            INTENT (IN) :: NSRC          !  actual source count
        INTEGER,            INTENT (IN) :: CTYPE         !  coord sys type
        CHARACTER(NAMLEN3), INTENT (IN) :: GDNAM         !  grid name
        REAL(8),            INTENT (IN) :: P_ALP         !  first, second, third map
        REAL(8),            INTENT (IN) :: P_BET         !  projection descriptive
        REAL(8),            INTENT (IN) :: P_GAM         !  parameters
        REAL(8),            INTENT (IN) :: XCENT         !  lon for coord-system X=0
        REAL(8),            INTENT (IN) :: YCENT         !  lat for coord-system Y=0
        REAL(8),            INTENT(IN OUT) :: XVALS( NSRC ) !  x location (input grid coord)
        REAL(8),            INTENT(IN OUT) :: YVALS( NSRC ) !  y location (input grid coord)

C...........   Other local variables
        INTEGER     UZONE               ! UTM zone
        INTEGER     S
        REAL        ALP, BET, GAM, XC, YC
        REAL        XLOC, XX, YLOC, YY  ! tmp x and y coordinates

        LOGICAL     EFLAG               ! true: error detected

        CHARACTER(NAMLEN3) TMPGDNAM     ! temporary grid name
        CHARACTER(256)     MESG         ! message buffer

        CHARACTER(16) :: PROGNAME = 'CONVRTLL' ! program name

C***********************************************************************
C   begin body of subroutine CONVRTLL

C.........  Copy input grid name to temporary variable since I/O API
C           routines may change it to the coordinate system name
        
        TMPGDNAM = GDNAM
        EFLAG    =.FALSE.
        ALP      = P_ALP
        BET      = P_BET
        GAM      = P_GAM
        XC       = XCENT
        YC       = YCENT

        IF ( CTYPE .EQ. LATGRD3 ) THEN
            RETURN

        ELSE IF ( CTYPE .EQ. UTMGRD3 ) THEN

            UZONE = NINT( P_ALP )

            DO S = 1, NSRC

                XLOC = XVALS( S )
                YLOC = YVALS( S )

                IF( XLOC .LT. AMISS3 .OR. YLOC .LT. AMISS3 ) CYCLE

                CALL UTM2LL( XLOC, YLOC, UZONE, XX, YY )
                XVALS( S ) = XX
                YVALS( S ) = YY

            ENDDO

        ELSE IF ( CTYPE .EQ. LAMGRD3 ) THEN

            IF( .NOT.LAMBERT( TMPGDNAM, ALP, BET, GAM, XC, YC ) ) THEN
                MESG = 'ERROR: Could not initialize Lambert grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            DO S = 1, NSRC

                XLOC = XVALS( S )
                YLOC = YVALS( S )

                IF( XLOC .LT. AMISS3 .OR. YLOC .LT. AMISS3 ) CYCLE

                IF ( .NOT. LAM2LL( XLOC, YLOC, XX, YY ) ) THEN
                    EFLAG = .TRUE.

                ELSE

                    XVALS( S ) = XX
                    YVALS( S ) = YY

                END IF

            ENDDO

            IF( EFLAG ) THEN
                MESG = 'ERROR: Problem converting coordinates ' //
     &                 'from Lambert to Lat-lon'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF

        ELSE IF ( CTYPE .EQ. POLGRD3 ) THEN

            IF( .NOT. POLSTE( TMPGDNAM, ALP, BET, GAM, XC, YC ) ) THEN
                MESG = 'ERROR: Could not initialize Polar ' //
     &                 'Stereographic grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            DO S = 1, NSRC

                XLOC = XVALS( S )
                YLOC = YVALS( S )

                IF( XLOC .LT. AMISS3 .OR. YLOC .LT. AMISS3 ) CYCLE

                IF ( .NOT. POL2LL( XLOC, YLOC, XX, YY ) ) THEN
                    EFLAG = .TRUE.

                ELSE

                    XVALS( S ) = XX
                    YVALS( S ) = YY

                END IF

            END DO

            IF( EFLAG ) THEN
                MESG = 'ERROR: Problem converting coordinates ' //
     &                 'from Polar Stereographic to Lat-lon'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2)
            
            END IF

        ELSE            !  error
            WRITE( MESG,94010 ) 'ERROR: Do not know how to convert ' //
     &             'from coordinate type number', CTYPE, 'to Lat-lon'

            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF          !  if coord type UTM, Lambert, or Polar Stereographic

        RETURN


C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93030   FORMAT( I8.8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE CONVRTLL
