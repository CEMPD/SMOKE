
        SUBROUTINE RDB3FAC( FDEV, NLINES, VGID, LINDX, LFAC, WNTF,
     &                      LWGT, FACS  ) 

C***********************************************************************
C  subroutine body starts at line XX 
C
C  DESCRIPTION:
C	Reads in the BEIS3 emissions factors from the BFAC file. 
C
C  PRECONDITIONS REQUIRED:
C
C  REVISION  HISTORY:
C	03/01 protoype by J. Vukovich
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C MCNC-Environmental Programs Group
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     !
        INCLUDE 'B3DIMS3.EXT'    ! biogenic parameters

C...........   ARGUMENTS and their descriptions: actually-occurring ASC table

        INTEGER, INTENT (IN)  :: FDEV    !  unit number for elev srcs file 
        INTEGER, INTENT (IN)  :: NLINES  !  no. veg types

        CHARACTER*16, INTENT (OUT)   :: VGID( NLINES )       ! veg ids
        INTEGER, INTENT (OUT)        :: LINDX( NLINES )      ! leaf area index
        REAL, INTENT (OUT)           :: LFAC( NLINES )       ! leaf biomass
        REAL, INTENT (OUT)           :: WNTF( NLINES )       ! winter factor
        REAL, INTENT (OUT)           :: LWGT( NLINES )       ! specific leaf wgt
        REAL, INTENT (OUT)           :: FACS( NLINES, NSEF ) ! emis facs
 
        LOGICAL      :: EFLAG = .FALSE.  !  error flag
        INTEGER       I, J               !  counters
        INTEGER       ISTAT              !  iostat error

        CHARACTER*300   MESG             !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDB3FAC' ! program name

C***********************************************************************
C   begin body of subroutine RDB3FAC

C.......... Read in emissions factors for each veg id

        DO I = 1, NLINES

          READ( FDEV, 93010, IOSTAT=ISTAT )
     &          VGID( I ), LINDX( I ), LFAC( I ), WNTF( I ),
     &          LWGT( I ) , ( FACS( I, J ) , J = 1, NSEF )   

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

93010   FORMAT( 8X, A16, 8x, I1, F8.0, F8.1, 5F8.0 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

        END SUBROUTINE RDB3FAC
