
        SUBROUTINE WREMFACS( NNAME, DNAME, NDATE, DDATE, JTIME, 
     &                       WFLAG_NDI, WFLAG_DIU )
   
C***********************************************************************
C  subroutine WREMFACS body starts at line < >
C
C  DESCRIPTION:
C      Ths subroutine write the emission factors to the non-diurnal and 
C      diurnal emission factor files.
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
C.........  This module contains emission factor tables and related
        USE MODEMFAC

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS
        INTEGER     TIME2SEC
        EXTERNAL    TIME2SEC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN):: NNAME     ! name of non-diur file
        CHARACTER(*), INTENT (IN):: DNAME     ! name of diur file
        INTEGER,      INTENT (IN):: NDATE     ! output date for non-di
        INTEGER,      INTENT (IN):: DDATE     ! output date for diurnl
        INTEGER,      INTENT (IN):: JTIME     ! output time (PSI)
        LOGICAL,      INTENT (IN):: WFLAG_NDI ! true: reusing nondiur
        LOGICAL,      INTENT (IN):: WFLAG_DIU ! true: reusing diur

C...........   Local variables
        INTEGER         L            ! counters and indices

        INTEGER         PSI          ! tmp parameter scheme index

        CHARACTER*300   MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'WREMFACS' ! program name

C***********************************************************************
C   begin body of subroutine WREMFACS

C.........  Convert second for file to parameter scheme index value
        PSI = TIME2SEC( JTIME )

C.........  Write all non-diurnal EFs to non-diurnal file
        L = LEN_TRIM( NNAME )
        IF( WFLAG_NDI ) THEN

            IF( .NOT. WRITE3( NNAME,ALLVAR3,NDATE,JTIME,EFACNDI ) ) THEN

                WRITE( MESG, 94010 ) 
     &                 'Could not write non-diurnal EFs for PSI', 
     &                 PSI, 'to file "' // NNAME( 1:L ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

        ELSE

            WRITE( MESG, 94010 ) 'NOTE: Non-diurnal EFs for PSI', PSI,
     &             'already exist in "' // NNAME( 1:L ) // '"'
            CALL M3MESG( MESG )

        END IF

C.........  Write all diurnal EFs to diurnal file
        L = LEN_TRIM( DNAME )
        IF( WFLAG_DIU ) THEN

            IF( .NOT. WRITE3( DNAME,ALLVAR3,DDATE,JTIME,EFACDIU ) ) THEN

                WRITE( MESG, 94010 ) 
     &                 'Could not write diurnal EFs for PSI',
     &                 PSI, 'to file "' // DNAME( 1:L ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

        ELSE

            WRITE( MESG, 94010 ) 'NOTE: Diurnal     EFs for PSI', PSI,
     &             'already exist in "' // DNAME( 1:L ) // '"'
            CALL M3MESG( MESG )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

        END SUBROUTINE WREMFACS
