
        SUBROUTINE RDBFIP( CDEV, NLINES, FCODES )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C	Reads in the counties from BFIP file.  These counties
C       will have corresponding landuse data to be read in later
C       from BCUSE.
C
C  PRECONDITIONS REQUIRED:
C
C  REVISION  HISTORY:
C	Written  9/99 by J. Vukovich
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

C...........   ARGUMENTS and their descriptions: actually-occurring ASC table

        INTEGER, INTENT (IN)  :: CDEV    !  unit number for elev srcs file 
        INTEGER, INTENT (IN)  :: NLINES  !  no. veg types
 
        REAL, INTENT (OUT)    :: FCODES( NLINES )

        INTEGER      COID
        INTEGER      STID
        INTEGER      CYID
 
        LOGICAL      :: EFLAG = .FALSE.  !  error flag

        CHARACTER*300   MESG             !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDBFIP' ! program name

C***********************************************************************
C   begin body of subroutine RDBFIP

        DO I = 1, NLINES

          READ( FDEV, 93010, IOSTAT=ISTAT ) LINE

          IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'Error', ISTAT,
     &              'reading BFIP file at line', I
                CALL M3MESG( MESG )
          ELSE
             
           IF ( LINE ( 1 ) .EQ. ' ' ) THEN
             COID  = 0 
           ELSE
             COID  = STR2INT( LINE( 1 ) )
           ENDIF

           IF ( LINE ( 2 : 3 ) .EQ. '  ' ) THEN
             STID = 0
           ELSE
             STID = STR2INT( LINE( 2 : 3 ) )
           ENDIF

           IF ( LINE ( 4 : 6 ) .EQ. '   ' ) THEN
             CYID = 0
           ELSE
             CYID = STR2INT( LINE( 4 : 6 ) )
           ENDIF

           FCODES( I ) = COID * 100000 + STID * 1000 + CYID
  
          END IF
       
        ENDDO

        IF( EFLAG ) THEN
            MESG = 'Problem reading biogenic emissions county file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

        END SUBROUTINE RDPELV

