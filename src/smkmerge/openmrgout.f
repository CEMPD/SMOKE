
        SUBROUTINE OPENMRGOUT( NGRP )

C***********************************************************************
C  subroutine OPENMRGOUT body starts at line 86
C
C  DESCRIPTION:
C      The purpose of this subroutine is to open all of the necessary
C      output files for the merge routine (both I/O API and ASCII files)
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
C
C***********************************************************************
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        INTEGER         INDEX1
        CHARACTER*16    MULTUNIT
        INTEGER         PROMPTFFILE  
        CHARACTER*16    PROMPTMFILE  
        CHARACTER*16    VERCHAR

        EXTERNAL  CRLF, INDEX1, MULTUNIT, PROMPTFFILE, PROMPTMFILE

C...........  SUBROUTINE ARGUMENTS
       INTEGER, INTENT (IN) :: NGRP     ! Actual number of groups

C...........  Local parameters
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C.........  Base and future year per 

C.........  Other local variables

        INTEGER         I, J, K, L, L1, L2, LD, LJ, N, V

        CHARACTER*50    RTYPNAM      ! name for type of state/county report file
        CHARACTER*50    UNIT         ! tmp units buffer
        CHARACTER*300   BUFFER       ! tmp buffer
        CHARACTER*300   MESG         ! message buffer
        CHARACTER(LEN=IODLEN3) DESCBUF ! variable description buffer

        CHARACTER*16 :: PROGNAME = 'OPENMRGOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENMRGOUT

C.........  Set default output file names
        CALL MRGONAMS     

C.........  Set up header for I/O API output files
        FTYPE3D = GRDDED3
        P_ALP3D = DBLE( P_ALP )
        P_BET3D = DBLE( P_BET )
        P_GAM3D = DBLE( P_GAM )
        XCENT3D = DBLE( XCENT )
        YCENT3D = DBLE( YCENT )
        XORIG3D = DBLE( XORIG )
        YORIG3D = DBLE( YORIG )
        XCELL3D = DBLE( XCELL )
        YCELL3D = DBLE( YCELL )
        SDATE3D = SDATE
        STIME3D = STIME
        TSTEP3D = TSTEP
        NCOLS3D = NCOLS
        NROWS3D = NROWS
        NTHIK3D = 1
        GDTYP3D = GDTYP
        VGTYP3D = VGTYP
        VGTOP3D = VGTOP
        GDNAM3D = GRDNM

        FDESC3D = ' '   ! array
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        WRITE( FDESC3D( 5 ),94010 ) '/BASE YEAR/ ', BYEAR 
        IF( PYEAR .NE. BYEAR ) 
     &      WRITE( FDESC3D( 6 ),94010 ) '/PROJECTED YEAR/ ', PYEAR
        
C.........  Set ozone-season description buffer
        DESCBUF = ' '
        IF( LO3SEAS ) DESCBUF = 'Ozone-season value, '
        LD = MAX( LEN_TRIM( DESCBUF ), 1 )

C.........  Set up and open I/O API output file
        IF( LGRDOUT ) THEN
          
C.............  Prompt for and gridded open file(s)
            IF( AFLAG ) THEN
                CALL SETUP_VARIABLES( ANIPOL, ANMSPC, AEINAM, AEMNAM )
                NLAYS3D = 1
                FDESC3D( 1 ) = 'Area source emissions data'
                AONAME = PROMPTMFILE(  
     &            'Enter name for AREA-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, AONAME, PROGNAME )
            END IF 

            IF( BFLAG ) THEN
                CALL SETUP_VARIABLES( BNIPOL, BNMSPC, BEINAM, BEMNAM )
                NLAYS3D = 1
                FDESC3D( 1 ) = 'Biogenic source emissions data'
                BONAME = PROMPTMFILE(  
     &            'Enter name for BIOGENIC-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, BONAME, PROGNAME )
            END IF 

            IF( MFLAG ) THEN
                CALL SETUP_VARIABLES( MNIPPA, MNMSPC, MEANAM, MEMNAM )
                NLAYS3D = 1
                FDESC3D( 1 ) = 'Mobile source emissions data'
                MONAME = PROMPTMFILE(  
     &            'Enter name for MOBILE-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, MONAME, PROGNAME )
            END IF 

            IF( PFLAG ) THEN
                CALL SETUP_VARIABLES( PNIPOL, PNMSPC, PEINAM, PEMNAM )
                NLAYS3D = EMLAYS
                IF( ALLOCATED( VGLVS ) ) THEN
                    J = LBOUND( VGLVS3D,1 )
                    DO V = 0, EMLAYS
                        VGLVS3D( J ) = VGLVS( V )
                        J = J + 1
                    END DO
                ELSE
                    VGLVS3D = 0  ! array
                ENDIF
                FDESC3D( 1 ) = 'Point source emissions data'
                PONAME = PROMPTMFILE(  
     &            'Enter name for POINT-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, PONAME, PROGNAME )
            END IF 

            IF( XFLAG ) THEN
                CALL SETUP_VARIABLES( NIPPA, NMSPC, EANAM, EMNAM )
                NLAYS3D = EMLAYS
                IF( ALLOCATED( VGLVS ) ) THEN
                    J = LBOUND( VGLVS3D,1 )
                    DO V = 0, EMLAYS
                        VGLVS3D( J ) = VGLVS( V )
                        J = J + 1
                    END DO
                ELSE
                   VGLVS3D = 0  ! array
                ENDIF
                FDESC3D( 1 ) = 'Multiple category emissions data'
                TONAME = PROMPTMFILE(  
     &            'Enter name for MULTI-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, TONAME, PROGNAME )
            END IF 

        END IF  ! End of gridded output

C.........  Open plume-in-grid output
        IF( PINGFLAG ) THEN

C.............  Override gridded file settings
            NCOLS3D = 1
            NROWS3D = NGROUP
            NLAYS3D = 1
            GDTYP3D = GDTYP
            VGTYP3D = IMISS3
            VGTOP3D = BADVAL3

            FDESC3D = ' '   ! array
            
            PINGNAME = PROMPTMFILE( 
     &                     'Enter name for PING EMISSIONS OUTPUT file',
     &                     FSUNKN3, PINGNAME, PROGNAME )
        END IF

C.........  Open plume-in-grid output
        IF( ELEVFLAG .AND. NGRP .GT. 1 ) THEN

            WRITE( MESG, 94010 ) 'The number of processing groups is ',
     &             NGRP, ', but it cannot be greater than' // CRLF() //
     &             BLANK10 // '1 for ASCII elevated output.  Not ' //
     &             'enough memory.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE IF( ELEVFLAG ) THEN
            
            EVDEV = PROMPTFFILE(
     &              'Enter name for ASCII ELEVATED SOURCES file', 
     &              .FALSE., .TRUE., PELVNAME, PROGNAME )

        END IF

C.........  Open report file(s)

        IF( LREPSTA .OR. LREPCNY ) THEN

            IF( LREPSTA .AND. LREPCNY ) THEN
                RTYPNAM = 'STATE AND COUNTY'
            ELSE IF ( LREPSTA ) THEN
                RTYPNAM = 'STATE'
            ELSE IF ( LREPCNY ) THEN
                RTYPNAM = 'COUNTY'
            END IF

            L = LEN_TRIM( RTYPNAM )

            IF( AFLAG ) THEN
                ARDEV  = PROMPTFFILE(
     &                  'Enter name for AREA-SOURCE ' // 
     &                  RTYPNAM( 1:L ) // ' REPORT', 
     &                  .FALSE., .TRUE., AREPNAME, PROGNAME )
            END IF 

            IF( BFLAG ) THEN
                BRDEV  = PROMPTFFILE(
     &                  'Enter name for BIOGENIC ' // 
     &                  RTYPNAM( 1:L ) // ' REPORT', 
     &                  .FALSE., .TRUE., BREPNAME, PROGNAME )
            END IF 

            IF( MFLAG ) THEN
                MRDEV  = PROMPTFFILE(
     &                  'Enter name for MOBILE-SOURCE ' // 
     &                  RTYPNAM( 1:L ) // ' REPORT', 
     &                  .FALSE., .TRUE., MREPNAME, PROGNAME )
            END IF 

            IF( PFLAG ) THEN
                PRDEV  = PROMPTFFILE(
     &                  'Enter name for POINT-SOURCE ' // 
     &                  RTYPNAM( 1:L ) // ' REPORT', 
     &                  .FALSE., .TRUE., PREPNAME, PROGNAME )
            END IF 

            IF( XFLAG ) THEN
                TRDEV  = PROMPTFFILE(
     &                  'Enter name for TOTAL ' // 
     &                  RTYPNAM( 1:L ) // ' REPORT', 
     &                  .FALSE., .TRUE., TREPNAME, PROGNAME )
            END IF 

        END IF  ! End of state and/or county output

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

C*****************  INTERNAL SUBPROGRAMS  ******************************

        CONTAINS

C.............  This internal subprogram uses WRITE3 and exits gracefully
C               if a write error occurred
            SUBROUTINE SETUP_VARIABLES( NIPPA_L, NMSPC_L, 
     &                                  EANAM_L, EMNAM_L  )

C.............  MODULES for public variables
C.............  This module contains the major data structure and control flags
            USE MODMERGE

C.............  Internal subprogram arguments
            INTEGER     , INTENT (IN) :: NIPPA_L
            INTEGER     , INTENT (IN) :: NMSPC_L
            CHARACTER(*), INTENT (IN) :: EANAM_L( NIPPA_L )
            CHARACTER(*), INTENT (IN) :: EMNAM_L( NMSPC_L )

C.............  Local subprogram varibles
            INTEGER     I, J, K, LJ, M, V

            CHARACTER(LEN=IOVLEN3) CBUF

C------------------------------------------------------------------------

C.............  Set constants number and values for variables
C.............  Do this regardless of whether we have outputs or not
C.............  For speciation...
            IF( SFLAG ) THEN
        	NVARS3D = NMSPC_L

        	K = 0
        	LJ = -1
C.................  Loop through global species index
        	DO N = 1, NGRP
                    DO V = 1, VGRPCNT( N )

C.........................  Access global indices
                	I = SIINDEX( V,N )
                	J = SPINDEX( V,N )
                	IF( J .EQ. LJ ) CYCLE    ! Do not repeat species

C.........................  Make sure current species is in local array
                        CBUF = EMNAM( J )
                        M = INDEX1( CBUF, NMSPC_L, EMNAM_L )
                        IF( M .LE. 0 ) CYCLE

                	DESCBUF= DESCBUF(1:LD)//' Model species '// CBUF

                	K = K + 1
                	VNAME3D( K ) = CBUF
                	UNITS3D( K ) = GRDUNIT( I )
                	VDESC3D( K ) = ADJUSTL( DESCBUF )
                	VTYPE3D( K ) = M3REAL

                	LJ = J

                    END DO
        	END DO

C.............  For no speciation...
            ELSE

        	NVARS3D = NIPPA_L

        	K = 0
C.................  Loop through global list
        	DO V = 1, NIPPA

C.....................  Determine if variable is an emission type, pollutant,
C                       or activity, and set variable description accordingly
                    I = INDEX1( EANAM( V ), NIPOL, EINAM )
                    J = INDEX ( EANAM( V ), ETJOIN )
                    IF( J .GT. 0 ) THEN
                	DESCBUF = DESCBUF(1:LD) // 
     &                           ' Emission type ' // EANAM( V )
                    ELSE IF( I .GT. 0 ) THEN
                	DESCBUF = DESCBUF(1:LD) // 
     &                            ' Pollutant ' // EANAM( V )
                    ELSE
                	DESCBUF = DESCBUF(1:LD) // 
     &                            ' Activity ' // EANAM( V )
                    END IF

C.....................  Make sure current data variable is in local array
                    CBUF = EANAM( V )
                    M = INDEX1( CBUF, NIPPA_L, EANAM_L )
                    IF( M .LE. 0 ) CYCLE

                    K = K + 1

C.....................  Define variable information
                    VNAME3D( K ) = EANAM  ( V )
                    UNITS3D( K ) = GRDUNIT( V )
                    VDESC3D( K ) = ADJUSTL( DESCBUF )
                    VTYPE3D( K ) = M3REAL

        	END DO

            END IF

            RETURN

            END SUBROUTINE SETUP_VARIABLES

        END SUBROUTINE OPENMRGOUT
