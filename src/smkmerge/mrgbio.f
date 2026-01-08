
        SUBROUTINE MRGBIO( VNAME, FNAME, JDATE, JTIME, NGRID, 
     &                     UNITFAC, BIOARR, ALLARR )

C***********************************************************************
C  subroutine MRGBIO body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the 
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 11/99 by M. Houyoux
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C****************************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C...........  EXTERNAL FUNCTIONS
C       CHARACTER(2)    CRLF   
C        EXTERNAL        CRLF

C...........  SUBROUTINE ARGUMENTS

        CHARACTER(*), INTENT  (IN) :: VNAME           ! variable name to read
        CHARACTER(*), INTENT  (IN) :: FNAME           ! file name to read
        INTEGER     , INTENT  (IN) :: JDATE           ! Julian date
        INTEGER     , INTENT  (IN) :: JTIME           ! time
        INTEGER     , INTENT  (IN) :: NGRID           ! no. grid cells
        REAL        , INTENT  (IN) :: UNITFAC         ! units conv factor
        REAL        , INTENT (OUT) :: BIOARR( NGRID ) ! biogenic emissions
        REAL        , INTENT (OUT) :: ALLARR( NGRID ) ! merged emissions

C...........   Other local variables
        INTEGER         C, L    !  counters and indices

        CHARACTER(300)         :: MESG ! message buffer

        CHARACTER(16) :: PROGNAME = 'MRGBIO' ! program name

C***********************************************************************
C   begin body of subroutine MRGBIO

C.........  Read biogenic emissions
        IF( .NOT. READ3( FNAME, VNAME, ALLAYS3, 
     &                   JDATE, JTIME, BIOARR      ) ) THEN
            L = LEN_TRIM( VNAME )
            MESG = 'Could not read "'// VNAME( 1:L ) // 
     &             '" from ' // FNAME
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF   ! if read3() failed

C.........  Update all emissions data array with biogenics
        DO C = 1, NGRID
            BIOARR( C ) = BIOARR( C ) * UNITFAC
            ALLARR( C ) = ALLARR( C ) + BIOARR( C )
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE MRGBIO
