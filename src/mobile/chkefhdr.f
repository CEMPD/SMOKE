
        SUBROUTINE CHKEFHDR( NNAME, DNAME )
   
C***********************************************************************
C  subroutine CHKEFHDR body starts at line < >
C
C  DESCRIPTION:
C      Check the headers of the diurnal and non-diurnal emission factor
C      files
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

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'M5CNST3.EXT'   !  Mobile5a/b constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        REAL                   GETRFDSC
        CHARACTER(LEN=IODLEN3) GETCFDSC

        EXTERNAL CRLF, GETRDSC, GETCFDSC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: NNAME  ! Name of non-diurnal file
        CHARACTER(*), INTENT (IN) :: DNAME  ! Name of diurnal file

C...........   Local variables
        INTEGER         I, L, L2       ! counters and indices

        INTEGER         JYEAR   ! 4-digit year

        REAL            NDTMIN    ! non-diurnal minimum temperature
        REAL            NDTMAX    ! non-diurnal maximum temperature
        REAL            NDTINT    ! non-diurnal temperature interval
        REAL            MINT_MIN  ! diurnal minimum temperature
        REAL            MAXT_MAX  ! diurnal maximum temperature
        REAL            TMMINVL   ! diurnal temperature interval

        LOGICAL      :: EFLAG = .FALSE.  ! true: processing error found

        CHARACTER(LEN=IODLEN3) BUF
        CHARACTER(LEN=IODLEN3) TUNIT   ! temperature units from non-diurnal file
        CHARACTER*300          MESG    ! message buffer

        CHARACTER*16 :: PROGNAME = 'CHKEFHDR' ! program name

C***********************************************************************
C   begin body of subroutine CHKEFHDR

C.........  Get information from non-diurnal file header
        IF( .NOT. DESC3( NNAME ) ) THEN
            MESG = 'Could not get description of ' // NNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Store important information for non-diurnal file header
        JYEAR  = SDATE3D / 1000
        NDTMIN = GETRFDSC( FDESC3D, '/MINT_MIN/'  , .TRUE. )
        NDTMAX = GETRFDSC( FDESC3D, '/MAXT_MAX/'  , .TRUE. )
        NDTINT = GETRFDSC( FDESC3D, '/T_INTERVAL/', .TRUE. )
        TUNIT  = GETCFDSC( FDESC3D, '/T_UNITS/'   , .TRUE. )

C.........  Get information from diurnal file header
        IF( .NOT. DESC3( DNAME ) ) THEN
            MESG = 'Could not get description of ' // DNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Check modeling year to ensure same as ENAME
        I = SDATE3D / 1000
        IF( I .NE. JYEAR ) THEN

            WRITE( MESG,94010 ) 
     &            'Year of non-diurnal emission factors in file "' //
     &             NNAME( 1:LEN_TRIM( NNAME ) ) // '": ', JYEAR,
     &             CRLF() // BLANK5 //
     &            'Year of diurnal     emission factors in file "' //
     &             DNAME( 1:LEN_TRIM( DNAME ) ) // '": ', I
            CALL M3MSG2( MESG )

            MESG = 'Incompatible input files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF 

C.........  Check that units are the same
        BUF = GETCFDSC( FDESC3D, '/T_UNITS/', .TRUE. )
        IF( TUNIT .NE. BUF ) THEN

            EFLAG = .TRUE.
            L  = LEN_TRIM( TUNIT )
            L2 = LEN_TRIM( BUF )
            MESG = 'Non-diurnal temperature units: ' // TUNIT( 1:L ) //
     &             CRLF() // BLANK5 //
     &             '    Diurnal temperature units: ' // BUF( 1:L2 )
            CALL M3MSG2( MESG )

        END IF

C.........  Check that file min/max temperatures and interval are consistent
C           in both files
        MINT_MIN = GETRFDSC( FDESC3D, '/MINT_MIN/', .TRUE. )
        MAXT_MAX = GETRFDSC( FDESC3D, '/MAXT_MAX/', .TRUE. )
        TMMINVL  = GETRFDSC( FDESC3D, '/T_INTERVAL/', .TRUE. )

        IF( NDTMIN .NE. MINT_MIN .OR.
     &      NDTMAX .NE. MAXT_MAX .OR.
     &      NDTINT .NE. TMMINVL       ) THEN

            EFLAG = .TRUE.
            L = LEN_TRIM( TUNIT )
            WRITE( MESG,94100 ) 
     &           'Non-diurnal temperatures [' // TUNIT(1:L) // ']:' //
     &           'Min=', NDTMIN, 'Max=', NDTMAX, 'Interval=', NDTINT,
     &            CRLF() // BLANK5 //
     &           '    Diurnal temperatures [' // TUNIT(1:L) // ']:' //
     &           'MinT_min=', MINT_MIN, 'MaxT_max=', MAXT_MAX,
     &           'Interval=', TMMINVL

            CALL M3MSG2( MESG )

        END IF

C.........  If there were errors, abort
        IF( EFLAG ) THEN

            WRITE( MESG,94010 ) 
     &             'Inconsistent headers in non-diurnal and diurnal ' //
     &             'emission factors files.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

94100   FORMAT( 10( A, :, F9.2, :, 1X ) )

        END SUBROUTINE CHKEFHDR
