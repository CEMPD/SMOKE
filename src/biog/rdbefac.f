
        SUBROUTINE RDBEFAC( FDEV, NLINES, VGID, EFACS, LINDX )

C***********************************************************************
C  subroutine body starts at line 69 
C
C  DESCRIPTION:
C	Reads in the biogenic emissions factors 
C       from the BFAC file.   
C
C  PRECONDITIONS REQUIRED:
C
C  REVISION  HISTORY:
C	Written  11/99 by J. Vukovich
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'     ! emissions constant parameters
        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'BIODIMS3.EXT'    ! biogenic parameters

C...........   ARGUMENTS and their descriptions: actually-occurring ASC table

        INTEGER, INTENT (IN)  :: FDEV    !  unit number for elev srcs file 
        INTEGER, INTENT (IN)  :: NLINES  !  no. veg types
 
        CHARACTER(LEN=BVGLEN3), INTENT (OUT)    :: VGID( NLINES )
        INTEGER, INTENT (OUT)    :: LINDX( NLINES )
        REAL, INTENT (OUT)    :: EFACS( NLINES, NSEF ) 
 
        LOGICAL      :: EFLAG = .FALSE.  !  error flag
        INTEGER       I, J               !  counters
        INTEGER       ISTAT              !  iostat error

        CHARACTER*300   MESG             !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDBEFAC' ! program name

C***********************************************************************
C   begin body of subroutine RDBEFAC

C.......... Read in emissions factors for each veg id

        DO I = 1, NLINES

          READ( FDEV, 93010, IOSTAT=ISTAT )
     &          VGID( I ),
     &        ( EFACS( I, J ) , J = 1, NSEF ),
     &          LINDX( I )
          IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'Error', ISTAT,
     &              'reading EMISSION FACTOR file at line', I
                CALL M3MESG( MESG )
          END IF

        ENDDO

        IF( EFLAG ) THEN
            MESG = 'Problem reading biogenic emissions factors file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( 1X, A4, 3F9.0, F6.0, I2 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

        END SUBROUTINE RDBEFAC

