
        SUBROUTINE EFCOMBO( NNAME, DNAME, NCNT, NDATE, DDATE,
     &                      CNTRBPSI, FACTORS                 )

C***********************************************************************
C  program body starts at line 108
C
C  DESCRIPTION:
C      This subroutine creates emission factors using multiple previously 
C      calculated emission factors and percentages for each
C
C  PRECONDITIONS REQUIRED:
C      Must already have pure emission factors from which can calculate
C      combination factors.  Must have arrays that specify how many and
C      which pure EFs and percentages for each.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Copied from efcombo.F 7/99 by M. Houyoux
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

C...........   MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

C.........  This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

        IMPLICIT NONE

        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures.

C.........  External function

        INTEGER       SEC2TIME

        EXTERNAL      SEC2TIME

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN):: NNAME     ! non-diurnal EF file
        CHARACTER(*), INTENT (IN):: DNAME     ! diurnal EF file
        INTEGER     , INTENT (IN):: NCNT      ! no. contributing PSIs in combo
        INTEGER,      INTENT (IN):: NDATE     ! output date for non-di
        INTEGER,      INTENT (IN):: DDATE     ! output date for diurnl
        INTEGER     , INTENT (IN):: CNTRBPSI( NCNT ) ! contributing PSIs
        REAL        , INTENT (IN):: FACTORS ( NCNT ) ! combo factors

C.........  Input emission factor arrays
        REAL    NDEF( NTMPR  , NVTYPE, NNDI )  ! Non-diurnal EF input array
        REAL    DEF ( NVLDTMM, NVTYPE, NDIU )  ! Diurnal     EF input array

C.........  Local variables

        INTEGER         I            ! Counters and pointers

        INTEGER         CPSI         ! tmp contributing PSI
        INTEGER         IPTIM        ! PSI converted to seconds

        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER*300   MESG

        CHARACTER*16 :: PROGNAME = 'EFCOMBO'   ! program name

C***********************************************************************
C   begin body of program EFCOMBO

C.........  Initialize emission factor output variables
        EFACNDI = 0.   ! array
        EFACDIU = 0.   ! array

C.........  Loop over component PSI for combo PSI
        DO I = 1, NCNT

C.............  Convert PSI code to seconds to read from file
            CPSI = CNTRBPSI( I )
            IPTIM = SEC2TIME( CPSI )

C.............  Read in each required existing Non-diurnal EF set
            IF( .NOT. READ3
     &          ( NNAME, 'ALL', 1, NDATE, IPTIM, NDEF ) ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &                 'Could not read non-diurnal EFs for PSI',
     &                 CPSI, 'from file "' // NNAME // '"'
                CALL M3MESG( MESG )

            ENDIF 

C.............  Read in each required existing Diurnal EF set
            IF( .NOT. READ3
     &          ( DNAME, 'ALL', 1, DDATE, IPTIM, DEF ) ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &                 'Could not read diurnal EFs for PSI',
     &                 CPSI, 'from file "' // DNAME // '"'
                CALL M3MESG( MESG )

            ENDIF 

            IF( EFLAG ) THEN
                MESG = 'Problem creating combination emission factor'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            EFACNDI = EFACNDI + NDEF * FACTORS( I )
            EFACDIU = EFACDIU + DEF  * FACTORS( I )
 
        END DO

        RETURN 

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 20 ( A, :, I5, :, 2X ) )

        END SUBROUTINE EFCOMBO
