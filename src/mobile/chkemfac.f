
        LOGICAL FUNCTION CHKEMFAC( EMFACTYP, FNAME, JDATE, JTIME,
     &                             REUSE, READFLAG )
   
C***********************************************************************
C  function CHKEMFAC body starts at line < >
C
C  DESCRIPTION:
C      
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
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
C...........   This module contains emission factor tables and related
        USE MODEMFAC

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'NETCDF.EXT'    !  NetCDF parameters

C...........   EXTERNAL FUNCTIONS
        INTEGER      TIME2SEC
        EXTERNAL     TIME2SEC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN):: EMFACTYP        ! DIURNAL or NON-DIURNAL
        CHARACTER(*), INTENT (IN):: FNAME           ! logical file name
        INTEGER,      INTENT (IN):: JDATE           ! julian date
        INTEGER,      INTENT (IN):: JTIME           ! time for read (PSI)
        LOGICAL,      INTENT (IN):: REUSE           ! true: factors to be reused
        LOGICAL,      INTENT (IN):: READFLAG        ! true: read EMFACS

C...........   Local variables
        INTEGER         I, L    ! counters and indices
        INTEGER         PSI     ! tmp parameter scheme index

        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER*300   MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'CHKEMFAC'

C***********************************************************************
C   begin body of function CHKEMFAC

C.........  Turn off verbose NetCDF to elimiate error messages when using
C           CHECK3 to check for available time steps
        CALL NCPOPT( NCNOERR )

C.........  Initialize exit status as 
        CHKEMFAC = .TRUE.

        IF( .NOT. REUSE ) THEN

            CALL SET_BADVAL3
            CHKEMFAC = .FALSE.

C.........  Check if the current PSI is in the file
        ELSE IF( .NOT. CHECK3( FNAME, ALLVAR3, JDATE, JTIME ) ) THEN

            CALL SET_BADVAL3
            CHKEMFAC = .FALSE.

C.............  If diurnal factors already exist, initialize by reading
C               existing values
        ELSE IF ( READFLAG ) THEN

            SELECT CASE( EMFACTYP )
            CASE( 'NON-DIURNAL' )

C.................  This is a placeholder which should be used when non-diurnal
C                   emission factors are computed only for the temperatures
C                   of interest and not for all of the temperatures as is 
C                   done now.

c               DO I = 1, NNDI

c               END DO

C.............  Read diurnal emission factors
            CASE( 'DIURNAL' )

                DO I = 1, NDIU

                    IF( .NOT. READ3( FNAME, DIUNAM( I ), ALLAYS3,
     &                               JDATE, JTIME, EFACDIU( 1,1,I ) ) )
     &                  THEN

                        EFLAG = .TRUE.
                        PSI   = TIME2SEC( JTIME )
                        L     = LEN_TRIM( DIUNAM( I ) )
                        WRITE( MESG,94010 ) 
     &                         'ERROR: Could not read PSI', PSI, 
     &                         'for "' // DIUNAM( I )( 1:L ) // 
     &                         '" in diurnal EFs file'
                        CALL M3MSG2( MESG )

                    END IF

                END DO

            CASE DEFAULT

                CALL BAD_CALL_EXIT

            END SELECT

        END IF  ! End if diurnal EFs already computed or not.

C.........  Set exit status depending on the error flag
        IF( EFLAG ) THEN
            MESG = 'Problem checking or reading emission factors file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG,  2 )
        END IF

C.........  Turn back on verbose NetCDF 
        CALL NCPOPT( NCVERBOS )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

        CONTAINS

            SUBROUTINE SET_BADVAL3

            SELECT CASE( EMFACTYP )
            CASE( 'NON-DIURNAL' )
                EFACNDI = BADVAL3

            CASE( 'DIURNAL' )
                EFACDIU = BADVAL3  ! array

            CASE DEFAULT

                CALL BAD_CALL_EXIT

            END SELECT

            END SUBROUTINE SET_BADVAL3

C------------------------------------------------------------------------

            SUBROUTINE BAD_CALL_EXIT

            MESG = 'INTERNAL ERROR: Emission factor type "' // 
     &             EMFACTYP // '" not recognized by routine' // 
     &             PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END SUBROUTINE BAD_CALL_EXIT

        END FUNCTION CHKEMFAC
