
        SUBROUTINE CONVRTXY( NSRC, XVALS, YVALS, CTYPE, P_ALP, P_BET, 
     &                       P_GAM, XCENT, YCENT )

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
C
C****************************************************************************/
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

C...........   INCLUDES

        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C.........  SUBROUTINE ARGUMENTS
        INTEGER       NSRC          !  actual source count
        REAL          XVALS( NSRC ) !  x location (input lon)
        REAL          YVALS( NSRC ) !  y location (input lat)
        INTEGER       CTYPE	    !  coord sys type
        REAL*8        P_ALP	    !  first, second, third map
        REAL*8        P_BET	    !  projection descriptive
        REAL*8        P_GAM	    !  parameters
        REAL*8        XCENT	    !  lon for coord-system X=0
        REAL*8        YCENT	    !  lat for coord-system Y=0

C...........   EXTERNAL FUNCTIONS
        LOGICAL       LAMBERT
        LOGICAL       LL2LAM 

        EXTERNAL      LAMBERT, LL2LAM

C...........   Other local variables
        INTEGER       UZONE               ! UTM zone
        INTEGER       S

        REAL          XLOC, XX, YLOC, YY  ! tmp x and y coordinates

        LOGICAL       EFLAG    ! error flag

        CHARACTER*256 MESG                !  message buffer

        CHARACTER*16 :: PROGNAME = 'CONVRTXY' ! program name

C***********************************************************************
C   begin body of subroutine CONVRTXY

        IF ( CTYPE .EQ. UTMGRD3 ) THEN

            UZONE = NINT( P_ALP )

            DO S = 1, NSRC

                XLOC = XVALS( S )
                YLOC = YVALS( S )
                CALL LL2UTM( XLOC, YLOC, UZONE, XX, YY )
                XVALS( S ) = XX
                YVALS( S ) = YY

            ENDDO

        ELSE IF ( CTYPE .EQ. LAMGRD3 ) THEN

            IF( .NOT. LAMBERT( P_ALP, P_BET, P_GAM, XCENT, YCENT )) THEN
                MESG = 'ERROR: Could not initialize lambert grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            DO S = 1, NSRC

                XLOC = XVALS( S )
                YLOC = YVALS( S )

                IF ( .NOT. LL2LAM( XLOC, YLOC, XX, YY ) ) THEN
                    EFLAG = .TRUE.

                ELSE

                    XVALS( S ) = XX
                    YVALS( S ) = YY

                END IF

            ENDDO

            IF( EFLAG ) THEN
                MESG = 'ERROR: Problem converting stack coordinates ' //
     &                 'from Lat-lon to Lambert'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF

        ELSE		!  error
            WRITE( MESG,94010 ) 'ERROR: Do not know how to convert ' //
     &             'from Lat-lon to coordinate type number', CTYPE

            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF		!  if coord type UTM or Lambert 

        RETURN


C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93030   FORMAT( I8.8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
