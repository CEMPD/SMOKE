
        SUBROUTINE PROCSPRO( NMSPC, SPCNAM )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Abridge full pollutants table to a data structure that can be used
C      in ASGNSPRO
C
C  PRECONDITIONS REQUIRED:
C      Expects cross-reference tables to be set to IMISS3 if not defined
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 2/99 by M. Houyoux
C
C****************************************************************************/
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
C***************************************************************************

C...........   MODULES for public variables   
C...........   This module contains the speciation profile tables
        USE MODSPRO, ONLY: NSPFUL, NSPROF, NSPECIES, SPROFN,
     &                     IDXSPRO, IDXSSPEC, INPRF, SPECID

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         FINDC
        EXTERNAL        FINDC

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NMSPC           ! No. of species for 1 pol
        CHARACTER(*), INTENT (IN) :: SPCNAM( NMSPC ) ! Species names for 1 pol
      
C.........  Other local variables
        INTEGER          I, J, K, L1      !  counters and indices
        INTEGER          IOS              !  i/o status

        CHARACTER(300)     MESG       ! message buffer

        CHARACTER(IOVLEN3) SBUF       ! species name buffer
        CHARACTER(SPNLEN3) PCODE      ! current speciation profile code
        CHARACTER(SPNLEN3) PREVCODE   ! previous speciation profile code

        CHARACTER(16) :: PROGNAME = 'PROCSPRO' ! program name

C***********************************************************************
C   begin body of subroutine PROCSPRO

C.........  Loop through sorted, unprocessed speciation profiles table to count
C           the number of unique profiles

        PREVCODE = EMCMISS3
        J = 0
        DO I = 1, NSPFUL

            PCODE = INPRF( I )

            IF( PCODE .NE. PREVCODE ) THEN
                J = J + 1
                PREVCODE = PCODE
            ENDIF

        END DO        !  end loop on speciation profile table

        NSPROF = J

C.........  Allocate memory for unique profiles arrays
        ALLOCATE( SPROFN( NSPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPROFN', PROGNAME )
        ALLOCATE( IDXSPRO( NSPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSPRO', PROGNAME )
        ALLOCATE( NSPECIES( NSPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NSPECIES', PROGNAME )
        ALLOCATE( IDXSSPEC( NSPROF,NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSSPEC', PROGNAME )

C.........  Initialize processed speciation profiles counter array
        NSPECIES = 0   ! array

C.........  Loop through sorted, unprocessed speciation profiles table for 
C           current pollutant to process it.  This algorithm works because we
C           have stored (in INPRF) only those profiles for the current
C           pollutant.

        PREVCODE = EMCMISS3
        J = 0
        DO I = 1, NSPFUL

            PCODE = INPRF ( I )
            SBUF  = SPECID( I )

            IF( PCODE .NE. PREVCODE ) THEN

                J = J + 1

                SPROFN  ( J ) = PCODE
                IDXSPRO ( J ) = I

                PREVCODE = PCODE

            END IF

            K = FINDC( SBUF, NMSPC, SPCNAM )

            IF( K .LE. 0 ) THEN
                L1 = LEN_TRIM( SBUF )
                MESG = 'INTERNAL ERROR: model species "' //
     &                 SBUF( 1:L1 ) // '" not found in sorted ' //
     &                 'names list developed for this pollutant!'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            ELSE

                NSPECIES( J ) = NSPECIES( J ) + 1
                IDXSSPEC( J,NSPECIES( J ) ) = K

            END IF

        END DO        !  end loop on speciation profile table

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE PROCSPRO
