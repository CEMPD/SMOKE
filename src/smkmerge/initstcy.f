
        SUBROUTINE INITSTCY

C***********************************************************************
C  subroutine INITSTCY body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to initialize the necessary fields
C      for performing state and county totals.  The first call sets up the
C      indices from each source to each county.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 8/99 by M. Houyoux
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

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         FIND1  

        EXTERNAL   CRLF, ENVYN, FIND1

C...........   Other local variables

        INTEGER          IOS      ! i/o status
        INTEGER          J        ! counter

        LOGICAL, SAVE :: FIRSTIME = .TRUE. ! true: first time routine called

        CHARACTER*300    MESG     ! message buffer

        CHARACTER*16  :: PROGNAME = 'INITSTCY' ! program name

C***********************************************************************
C   begin body of subroutine INITSTCY

C.........  Read surrogates (if needed) and state/county names
        IF( FIRSTIME ) THEN
        
C.............  Allocate memory for indices from Co/st/cy codes to counties
            ALLOCATE( AICNY( NASRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AICNY', PROGNAME )
            ALLOCATE( MICNY( NMSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MICNY', PROGNAME )
            ALLOCATE( PICNY( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PICNY', PROGNAME )

C.............  Create indices to counties from Co/st/cy codes for each source 
C               category
            IF( AFLAG ) THEN                
                CALL SET_COUNTY_INDEX( NASRC, AIFIP, AICNY )
            END IF

            IF( MFLAG ) THEN                
                CALL SET_COUNTY_INDEX( NMSRC, MIFIP, MICNY )
            END IF

            IF( PFLAG ) THEN                
                CALL SET_COUNTY_INDEX( NPSRC, PIFIP, PICNY )
            END IF

            FIRSTIME = .FALSE.

        END IF

C.........  Initialize totals to zero...
C.........  State totals...
        IF( LREPSTA ) THEN

            IF( AFLAG ) THEN
                AEBSTA = 0.   ! state total inven or speciated emissions
        	AEUSTA = 0.   ! state total multiplicative-controlled emissions
        	AEASTA = 0.   ! state total additive- controlled emissions
        	AERSTA = 0.   ! state total reactivity-controlled emissions
        	AECSTA = 0.   ! state total all-controlled emissions
            END IF

            IF( BFLAG ) THEN
                BEBSTA = 0.
            END IF

            IF( MFLAG ) THEN
        	MEBSTA = 0.
        	MEUSTA = 0.
        	MEASTA = 0.
        	MERSTA = 0.
        	MECSTA = 0.
            END IF

            IF( PFLAG ) THEN
        	PEBSTA = 0.
        	PEUSTA = 0.
        	PEASTA = 0.
        	PERSTA = 0.
        	PECSTA = 0.
            END IF

            IF( XFLAG ) THEN
        	TEBSTA = 0.
        	TEUSTA = 0.
        	TEASTA = 0.
        	TERSTA = 0.
        	TECSTA = 0.
            END IF

        END IF

C.........  County totals...
        IF( AFLAG ) THEN
            AEBCNY = 0.   ! county total inven or speciated emissions
            AEUCNY = 0.   ! county total multiplicative-controlled emissions
            AEACNY = 0.   ! county total additive- controlled emissions
            AERCNY = 0.   ! county total reactivity-controlled emissions
            AECCNY = 0.   ! county total all-controlled emissions
        END IF

        IF( BFLAG ) THEN
            BEBCNY = 0.
        END IF

        IF( MFLAG ) THEN
            MEBCNY = 0.
            MEUCNY = 0.
            MEACNY = 0.
            MERCNY = 0.
            MECCNY = 0.
        END IF

        IF( PFLAG ) THEN
            PEBCNY = 0.
            PEUCNY = 0.
            PEACNY = 0.
            PERCNY = 0.
            PECCNY = 0.
        END IF

        IF( XFLAG ) THEN
            TEBCNY = 0.
            TEUCNY = 0.
            TEACNY = 0.
            TERCNY = 0.
            TECCNY = 0.
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS   ******************************

        CONTAINS

C.............  This subroutine creates the county to src FIPs index
            SUBROUTINE SET_COUNTY_INDEX( NSRC, IFIP, ICNY )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: NSRC
            INTEGER, INTENT (IN) :: IFIP( NSRC )
            INTEGER, INTENT(OUT) :: ICNY( NSRC )

C.............  Local variables
            INTEGER     J, S     ! counters and indices
            INTEGER     FIP      ! tmp cy/st/co code
            INTEGER     PFIP     ! previous cy/st/co code

C----------------------------------------------------------------------------

            PFIP = -9
            DO S = 1, NSRC

                FIP = IFIP( S )

                IF( FIP .NE. PFIP ) THEN

                    J = MAX( FIND1( FIP, NCOUNTY, CNTYCOD ), 0 )
                    PFIP = FIP

                END IF

                ICNY( S ) = J

            END DO

            END SUBROUTINE SET_COUNTY_INDEX

        END SUBROUTINE INITSTCY
