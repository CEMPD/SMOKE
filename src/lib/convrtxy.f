
        SUBROUTINE CONVRTXY( NSRC, CTYPE, GDNAM, P_ALP, P_BET, 
     &                       P_GAM, XCENT, YCENT, XVALS, YVALS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine converts coordinates from lat-lon to the grid
C      defined by the subroutine arguments.
C
C  PRECONDITIONS REQUIRED:
C      
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Added POLGRD3 as a supported coord sys type - E. Giroux CNRC 03/2004
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

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C.........  SUBROUTINE ARGUMENTS
        INTEGER,          INTENT    (IN) :: NSRC          !  actual source count
        INTEGER,          INTENT    (IN) :: CTYPE         !  coord sys type
        CHARACTER(LEN=*), INTENT    (IN) :: GDNAM         !  grid name
        REAL*8 ,          INTENT    (IN) :: P_ALP         !  first, second, third map
        REAL*8 ,          INTENT    (IN) :: P_BET         !  projection descriptive
        REAL*8 ,          INTENT    (IN) :: P_GAM         !  parameters
        REAL*8 ,          INTENT    (IN) :: XCENT         !  lon for coord-system X=0
        REAL*8 ,          INTENT    (IN) :: YCENT         !  lat for coord-system Y=0
        REAL   ,          INTENT(IN OUT) :: XVALS( NSRC ) !  x location (input lon)
        REAL   ,          INTENT(IN OUT) :: YVALS( NSRC ) !  y location (input lat)

C...........   EXTERNAL FUNCTIONS
        LOGICAL       LAMBERT
        LOGICAL       LL2LAM
        LOGICAL       LL2POL
        LOGICAL       POLSTE

        EXTERNAL      LAMBERT, LL2LAM, LL2POL, POLSTE

C...........   Other local variables
        INTEGER       UZONE               ! UTM zone
        INTEGER       S

        REAL          XLOC, XX, YLOC, YY  ! tmp x and y coordinates
        REAL*4        XLOC4, YLOC4, XX4, YY4

        LOGICAL    :: EFLAG =.FALSE.   ! true: error detected

        CHARACTER*256 MESG                !  message buffer

        CHARACTER*16 :: PROGNAME = 'CONVRTXY' ! program name

C***********************************************************************
C   begin body of subroutine CONVRTXY

        IF ( CTYPE .EQ. LATGRD3 ) THEN
            RETURN

        ELSE IF ( CTYPE .EQ. UTMGRD3 ) THEN

            UZONE = NINT( P_ALP )

            DO S = 1, NSRC

                XLOC = XVALS( S )
                YLOC = YVALS( S )

                IF( XLOC .LT. AMISS3 .OR. YLOC .LT. AMISS3 ) CYCLE

                CALL LL2UTM( XLOC, YLOC, UZONE, XX, YY )
                XVALS( S ) = XX
                YVALS( S ) = YY

            ENDDO

        ELSE IF ( CTYPE .EQ. LAMGRD3 ) THEN

            IF( .NOT. LAMBERT( GDNAM, P_ALP, P_BET, P_GAM, 
     &                         XCENT, YCENT )) THEN
                MESG = 'ERROR: Could not initialize Lambert grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            DO S = 1, NSRC

                XLOC4 = XVALS( S )
                YLOC4 = YVALS( S )

                IF( XLOC4 .LT. AMISS3 .OR. YLOC4 .LT. AMISS3 ) CYCLE

                IF ( .NOT. LL2LAM( XLOC4, YLOC4, XX4, YY4 ) ) THEN
                    EFLAG = .TRUE.

                ELSE

                    XVALS( S ) = XX4
                    YVALS( S ) = YY4

                END IF

            ENDDO

            IF( EFLAG ) THEN
                MESG = 'ERROR: Problem converting coordinates ' //
     &                 'from Lat-lon to Lambert'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF

        ELSE IF ( CTYPE .EQ. POLGRD3 ) THEN

            IF( .NOT. POLSTE( GDNAM, P_ALP, P_BET, P_GAM, 
     &                        XCENT, YCENT ) ) THEN
                MESG = 'ERROR: Could not initialize Polar ' //
     &                 'Stereographic grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            DO S = 1, NSRC

                XLOC4 = XVALS( S )
                YLOC4 = YVALS( S )

                IF( XLOC4 .LT. AMISS3 .OR. YLOC4 .LT. AMISS3 ) CYCLE

                IF ( .NOT. LL2POL( XLOC4, YLOC4, XX4, YY4 ) ) THEN
                    EFLAG = .TRUE.

                ELSE

                    XVALS( S ) = XX4
                    YVALS( S ) = YY4

                END IF

            END DO

            IF( EFLAG ) THEN
                MESG = 'ERROR: Problem converting coordinates ' //
     &                 'from Lat-lon to Polar Stereographic'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2)
            
            END IF

        ELSE		!  error
            WRITE( MESG,94010 ) 'ERROR: Do not know how to convert ' //
     &             'from Lat-lon to coordinate type number', CTYPE

            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF		!  if coord type UTM, Lambert, or Polar Stereographic

        RETURN


C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93030   FORMAT( I8.8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
