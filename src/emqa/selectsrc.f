
        SUBROUTINE SELECTSRC

c The SELECTSRC routine is responsible for selecting the sources from the
c inventory records based on the settings for a given report.

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains Smkreport-specific settings
        USE MODREPRT

C.........  This module contains report arrays for each output bin
        USE MODREPBN

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  Local variables
        INTEGER         S          ! counters and indices

        INTEGER         IOS        ! i/o status


        CHARACTER*16 :: PROGNAME = 'SELECTSRC' ! program name

C***********************************************************************
C   begin body of subroutine SELECTSRC

c The inventory data are already required to have been read in prior to
c this routine being called

c The routine should only be called if the report contains SELECT statements 
c other than SUBGRID

c Processing steps:

c       - Loop through sources
c		- If SELECT REGION, search for source Region number in 
c                 group and if found, EXCLUDE it
c		- If ELEVATED used, EXCLUDE if a low-level source
c		- If NOELEVATED used, EXCLUDE if an elevated source
c		- If other SELECT (future), exclude as appropriate
c		- Keep a count of the number of selected sources
c       - End loop
c
c	- Allocate memory for a list of source IDs
c	- Create a list of source IDs of the selected sources.

        IF( ALLOCATED( OUTSRC ) ) DEALLOCATE( OUTSRC, OUTBIN, INDEXA )

        NOUTREC = NSRC
        NSRCDROP = 0
        ALLOCATE( OUTSRC( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTSRC', PROGNAME )
        ALLOCATE( OUTBIN( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTBIN', PROGNAME )

        ALLOCATE( INDEXA( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )

C.........  Select all sources
        DO S = 1, NSRC

            OUTSRC( S ) = S
            INDEXA( S ) = 1

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE SELECTSRC
