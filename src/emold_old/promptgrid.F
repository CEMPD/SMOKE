
C Version "@(#)$Id$ $Source$ $Date$ 

        LOGICAL FUNCTION PROMPTGRID()

C***********************************************************************
C  function body starts at line
C
C  RETURNS:  TRUE iff success
C
C  FUNCTION:
C	Prompts user repeatedly for grid/coordinate system name, then
C	uses DSCGRID or DSCOORD to put grid/coordinate system description
C	into FHEAD3.EXT data structures.
C
C  PRECONDITIONS REQUIRED:
C	Valid GRIDDEFS file
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C	DSCOORD, DSCGRID, GETSTR, GETYN, M3WARN
C
C  REVISION  HISTORY:
C	prototype 11/95 by CJC
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL		DSCOORD, DSCGRID, GETYN
        INTEGER		TRIMLEN

        EXTERNAL        DSCOORD, DSCGRID, GETYN, TRIMLEN

C...........   Local Variables:

        CHARACTER*16  	ANAME, CNAME

C***********************************************************************
C   begin body of function  PROMPTGRID

11	CONTINUE
            
            CALL GETSTR( 'Enter name for output coordinate system',
     &                   'UTM_17', GDNAM3D )

            IF ( DSCOORD( GDNAM3D, GDTYP3D, 
     &                          P_ALP3D, P_BET3D, P_GAM3D, 
     &                          XCENT3D, YCENT3D ) ) THEN

                XORIG3D = DBLE( BADVAL3 )
                YORIG3D = DBLE( BADVAL3 )
                XCELL3D = DBLE( BADVAL3 )
                YCELL3D = DBLE( BADVAL3 )
                NCOLS3D = IMISS3
                NROWS3D = IMISS3
                NTHIK3D = IMISS3

            ELSE	!  dscoord() failed; retry with dscgrid()

                IF ( .NOT. DSCGRID( GDNAM3D, ANAME, GDTYP3D, 	!  retry with dscgrid()
     &                              P_ALP3D, P_BET3D, P_GAM3D, 
     &                              XCENT3D, YCENT3D, XORIG3D, YORIG3D, 
     &                              XCELL3D, YCELL3D, 
     &                              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

                    WRITE( *,92000 ) 
     &              'Could not get description for coordinate system "'
     &              // CNAME( 1 : TRIMLEN( CNAME ) ) // '"'

                    IF ( GETYN( 'Try again?', .TRUE. ) ) THEN
                        GO TO  11
                    ELSE
                        CALL M3EXIT( 'RAWAREA', 0, 0, 
     &                               'Bad coordinate system', 2 )
                    END IF	!  if retry getstr() or not

                END IF		!  if dscgrid succeeded, or failed

            END IF		!  if dscoord() failed; end of loop

        GDNAM3D = CNAME

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000	FORMAT( 5X, A )

        END

