
        SUBROUTINE SELECTSRC( RCNT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C    The SELECTSRC routine is responsible for selecting the sources from the
C    inventory records based on the settings for a given report.  If no 
C    groups are used in the current report, then all sources are selected. 
C    Selected sources will have INDEXA( S ) = 1, while unselected sources will
C    have INDEXA( S ) = 0.
C
C  PRECONDITIONS REQUIRED:
C    The inventory data are already required to have been read in prior to
C    this routine being called
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 4/2001 by M Houyoux
C
C***********************************************************************
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
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CIFIP, SPPNLO, SPPNHI, SPPROF, SPFRAC, INDEXA

C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: RPT_, LREGION, NREGNGRP, REGNNAM,
     &                      PINGOUT3, ELEVOUT3, NOELOUT3, PSFLAG,
     &                      ALLRPT, SLFLAG, SSFLAG, NREGREC, EXCLDRGN

C.........  This module contains report arrays for each output bin
        USE MODREPBN, ONLY: NOUTREC, NSRCDROP, OUTSRC, OUTBIN, OUTSPRO, OUTSFAC

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: LMAJOR, LPING

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC

        USE MODSPRO, ONLY: MXSPEC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2) CRLF
        INTEGER      INDEX1
        INTEGER      FINDC

        EXTERNAL    CRLF, INDEX1, FINDC

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: RCNT    ! current report number

C...........  Local variables
        INTEGER         J, L, N, S      ! counters and indices
        INTEGER         IOS             ! i/o status
        INTEGER         REGNIDX         ! index to list of region groups for this rpt 

        LOGICAL      :: EFLAG = .FALSE.  ! True: error has been detected

        CHARACTER(FIPLEN3) CFIP                ! tmp country/state/county code
        CHARACTER(256)     MESG                ! message buffer

        CHARACTER(16) :: PROGNAME = 'SELECTSRC' ! program name

C***********************************************************************
C   begin body of subroutine SELECTSRC

C.........  Set report-specific local settings
        RPT_    = ALLRPT( RCNT )
        LREGION = .FALSE.
 
C.........  Determine current report has any groups
        REGNIDX = 0
        IF( RPT_%REGNNAM .NE. ' ' ) THEN
            REGNIDX = INDEX1( RPT_%REGNNAM, NREGNGRP, REGNNAM( : ) )
            LREGION = .TRUE.

            IF( REGNIDX .LE. 0 ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( RPT_%REGNNAM )
                MESG = 'INTERNAL ERROR: Group "'// RPT_%REGNNAM( 1:L )//
     &             '"' // CRLF() // BLANK10 // 'not found in list '//
     &             'of groups.'
                CALL M3MSG2( MESG )
            END IF

        END IF

C.........  NOTE - if adding group types, do it here and in errors below

C.........  Report internal errors if group name is not found
        IF ( EFLAG ) THEN
            MESG = 'Problems using groups'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF 

C.........  Deallocate arrays if they've been allocated before
        ALLOCATE( INDEXA( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )

C......... Initialize list of selected sources
        DO S = 1, NSRC
            INDEXA( S ) = 1
        END DO

C.........  If need to select sources, loop through sources to identify 
C.........  the ones not to output
        IF( LREGION .OR. RPT_%ELEVSTAT .GT. 0 ) THEN

            DO S = 1, NSRC

C.................  If using a region group, search for FIPS code in list
                IF( LREGION ) THEN
                    CFIP = CIFIP( S )
                    J = FINDC( CFIP, NREGREC( REGNIDX ), 
     &                         EXCLDRGN ( 1,REGNIDX )   )

                    IF ( J .GT. 0 ) INDEXA( S ) = 0

                END IF

C.................  If selecting elevated srcs...
                SELECT CASE( RPT_%ELEVSTAT )
                CASE( PINGOUT3 )
                    IF ( .NOT. LPING( S ) ) INDEXA( S ) = 0

                CASE( ELEVOUT3 )
                    IF ( .NOT. LMAJOR( S ) ) INDEXA( S ) = 0

                CASE( NOELOUT3 )
                    IF ( LMAJOR( S ) ) INDEXA( S ) = 0

                END SELECT

C.................  NOTE - this should always be at end of select statements
C.................  Keep a count of the number of selected sources
                IF( INDEXA( S ) .EQ. 0 ) NSRCDROP = NSRCDROP + 1

            END DO   ! End loop through sources

        END IF

C.........  Now create compressed list of sources, but leave INDEXA as is 
C           for use by REPMRGGRD

        IF( PSFLAG ) THEN
            NOUTREC = 0
            DO S = 1, NSRC

                IF ( INDEXA( S ) .EQ. 1 ) THEN
                    NOUTREC = NOUTREC + SPPNHI( S ) - SPPNLO( S ) + 1
               END IF

            END DO
        
        ELSE
        
            NOUTREC = NSRC - NSRCDROP
        
        END IF
        
        IF ( ALLOCATED( OUTSRC  ) )  DEALLOCATE( OUTSRC )
        IF ( ALLOCATED( OUTBIN  ) )  DEALLOCATE( OUTBIN )
        IF ( ALLOCATED( OUTSPRO ) )  DEALLOCATE( OUTSPRO )
        IF ( ALLOCATED( OUTSFAC ) )  DEALLOCATE( OUTSFAC )

        ALLOCATE( OUTSRC( NOUTREC ), 
     &            OUTBIN( NOUTREC ), 
     &           OUTSPRO( NOUTREC ), 
     &           OUTSFAC( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTSRC,OUTSFAC', PROGNAME )

C.........  Now create compressed list of sources, but leave INDEXA as is 
C           for use by REPMRGGRD

        IF( PSFLAG ) THEN

            J = 0            
            DO S = 1, NSRC

                IF( INDEXA( S ) .EQ. 0 ) CYCLE

                DO N = SPPNLO( S ), SPPNHI( S )
                    J = J + 1
                    OUTSRC(  J ) = S
                    OUTSPRO( J ) = SPPROF( N )
                    OUTSFAC( J ) = SPFRAC( N )
                END DO

            END DO
        
        ELSE

            J = 0            
            DO S = 1, NSRC

                IF( INDEXA( S ) .EQ. 0 ) CYCLE

                J = J + 1
                OUTSRC(  J ) = S
                OUTSPRO( J ) = 'NA'
                OUTSFAC( J ) = 1.0

            END DO
        
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE SELECTSRC
