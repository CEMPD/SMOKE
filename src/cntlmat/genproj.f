
        SUBROUTINE GENPROJ( PYEAR, RDEV, ENAME, USEPOL )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine processes the projection data and writes out
C      the projection matrix.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     
C
C************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         PROMPTFFILE

        EXTERNAL   PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: PYEAR  ! projection year
        INTEGER     , INTENT (IN) :: RDEV   ! unit numnber of report file
        CHARACTER(*), INTENT (IN) :: ENAME  ! emission inventory file name
        LOGICAL     , INTENT (IN) :: USEPOL( NIPPA )  ! true: pol used in pkt

C...........   Locally allocated arrays 
        INTEGER, ALLOCATABLE :: ISPRJ ( : ) ! projection control data table index
        REAL   , ALLOCATABLE :: PRJFAC( : ) ! projection factor

C...........  Local static arrays
        LOGICAL          LF   ( MXCHRS )      !  true: column should be output
        CHARACTER*20     CHARS( MXCHRS )      !  source fields for output

C...........   Logical names

        CHARACTER*16     PNAME      ! logical name for projection matrix


C...........   Other local variables

        INTEGER          J, K, L, S    ! counters and indices
        INTEGER          IDUM          ! dummy integer
        INTEGER          IOS           ! i/o error status
        INTEGER          NC            ! local number src chars

        LOGICAL       :: EFLAG    = .FALSE.   ! true: error has occurred
        LOGICAL, SAVE :: APPLFLAG = .FALSE.  ! true: something has been applied

        CHARACTER*16     VARNAM 
        CHARACTER*300    MESG                 ! message buffer

        CHARACTER*16  :: PROGNAME = 'GENPROJ' ! program name

C***********************************************************************
C   begin body of subroutine GENPROJ

C.........  Allocate memory for the projection matrix.  Use data
C           structures for point sources, but this routine can be used for area
C           sources or mobile sources as well. 

        ALLOCATE( ISPRJ( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPRJ', PROGNAME )
        ALLOCATE( PRJFAC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRJFAC', PROGNAME )

C.........  Determine which projection packet goes to each source.
C           Since the projection packet is not pollutant specific,
C           simply send the first pollutant in the pollutant list
C           to ASGNCNTL in order to avoid looping through all pollutants.

         CALL ASGNCNTL( NSRC, 1, 'PROJECTION', USEPOL, EANAM(1), 
     &                  IDUM, ISPRJ )

C..........  Write header for report.
         WRITE( RDEV, 93000 ) 'Processed as '// CATDESC// ' sources'
         WRITE( RDEV, 93000 ) 
     &          'Projection factors applied with /PROJECTION/ packet'

         WRITE( RDEV, 93390 ) '      from base year    ', BYEAR
         WRITE( RDEV, 93390 ) '      to projected year ', PYEAR

         WRITE( RDEV, 93000 ) '      to all pollutants uniformly'
         WRITE( RDEV, 93000 ) ' '

C.........  Loop through all sources and store projection information for
C           those that have it.  Otherwise, set projection factor=1.

c        ISPRJ = 1  ! array

        DO S = 1, NSRC

            K = ISPRJ( S )       ! index to projection data tables

            IF( K .GT. 0 ) THEN

C.................  Store projection factor
                PRJFAC( S ) = PRJFC( K )

C.................  Format source characteristic information
                CALL PARSCSRC( CSOURC( S ), MXCHRS, SC_BEGP, SC_ENDP, 
     &                         LF, NC, CHARS )
                NC = MIN( NC, NCHARS )

C.................  Write out projection information for all sources
C                   that are getting projected
                WRITE( MESG, 94015 ) 
     &               ( CHARS( J )( 1:SC_ENDP(J)-SC_BEGP(J)+1 ), J=1,NC )
                L = LEN_TRIM( MESG )
                WRITE( RDEV, 94020 ) MESG( 1:L ), PRJFAC( S )
                APPLFLAG = .TRUE.

            ELSE
C.................  If source does not have projection info., set to 1.

                PRJFAC( S ) = 1.0

            END IF

        END DO

        IF( .NOT. APPLFLAG ) THEN

            MESG = 'WARNING: No PROJECTION packet entries match ' //
     &             'inventory.'
            CALL M3MSG2( MESG )

            MESG = 'WARNING: Projection matrix will not be created!'
            CALL M3MSG2( MESG )

C.............  Write not into report file
            WRITE( RDEV, 93000 ) 
     &             'No projection packet entries matched the inventory.'

            RETURN

        END IF

C.........  Set up and open output projection matrices

        CALL OPENPMAT( ENAME, BYEAR, PYEAR, PNAME )

C.........  Write the projection matrix

C.........  Initialize message to use in case there is an error

        MESG = 'Problem writing to output file "' //
     &         PNAME( 1:LEN_TRIM( PNAME ) ) // '"'

        L = LEN_TRIM( MESG )

C.........  Write the I/O API variables for the non-speciation data

        IF( .NOT. WRITE3( PNAME, 'pfac', 0, 0, PRJFAC ) ) THEN
            MESG = MESG( 1:L ) // ' for variable "PRJFAC"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        DEALLOCATE( PRJFAC )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93390   FORMAT( A, I4.4 )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94015   FORMAT( 10( A, :, 1X ) )

94020   FORMAT( A, 1X, E13.5 )

C******************  INTERNAL SUBPROGRAMS  *****************************

        END SUBROUTINE GENPROJ
