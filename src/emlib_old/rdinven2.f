
        SUBROUTINE RDINVEN2( CATEGORY, FNAME, MXSRC, NSRCS, NIVARS, 
     &                       VNAMES, IFIP, VAR2, VAR3, VAR4, VAR5, 
     &                       TZONES, TPFLAG, INVYR, REAL1, REAL2, REAL3,
     &                       REAL4, REAL5, REAL6, INVVAL )

C***********************************************************************
C  subroutine body starts at line 98
C
C  DESCRIPTION:
C       Reads in point, area, or mobile source SMOKE inventory files
C
C  PRECONDITIONS REQUIRED:
C       File opened and exists
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C
C  REVISION  HISTORY:
C       Initial version 8/98 by M Houyoux
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
 
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C.........  SUBROUTINE ARGUMENTS

        CHARACTER*6  CATEGORY         ! emissions category (point,area,mobile)
        CHARACTER*16 FNAME            ! inventory file name
        INTEGER      MXSRC            ! maximum number of sources (incoming)
        INTEGER      NSRCS            ! actual number of sources (incoming)
        INTEGER      NIVARS           ! number of inventory pols or 1 for VMT
        CHARACTER*3  VNAMES( NIVARS ) ! names of inventory pols or "VMT"
        INTEGER      IFIP( NSRCS )    ! FIPs state/county codes
        INTEGER      VAR2( NSRCS )    ! PT: SCC     , AR: ASC7, MB: Road class
        INTEGER      VAR3( NSRCS )    ! PT: SIC     , AR: ASC3, MB: Veh type
        INTEGER      VAR4( NSRCS )    ! PT: plant ID,         , MB: Emis type
        INTEGER      VAR5( NSRCS )    ! PT: stack ID
        INTEGER      TZONES( NSRCS )  ! Time zone
        INTEGER      TPFLAG( NSRCS )  ! Temporal basis flag
        INTEGER      INVYR ( NSRCS )  ! Inventory year
        REAL         REAL1 ( NSRCS )  ! PT: stack X location, , MB: Link beg X
        REAL         REAL2 ( NSRCS )  ! PT: stack Y location, , MB: Link beg Y
        REAL         REAL3 ( NSRCS )  ! PT: stack height      , MB: Link end X
        REAL         REAL4 ( NSRCS )  ! PT: stack diameter    , MB: Link end Y
        REAL         REAL5 ( NSRCS )  ! PT: stack exit temperature
        REAL         REAL6 ( NSRCS )  ! PT: stack exit velocity
        REAL         INVVAL( MXSRC,NIVARS )! PT & AR: inv pols, MB: VMT 

C.........  EXTERNAL FUNCTIONS
        INTEGER      TRIMLEN

        EXTERNAL     TRIMLEN

C.........  LOCAL VARIABLES
        INTEGER       S, V
        INTEGER       NVARI

        CHARACTER*16  FIPNAM
        CHARACTER*16  VR2NAM
        CHARACTER*16  VR3NAM
        CHARACTER*16  VR4NAM
        CHARACTER*16  VR5NAM
        CHARACTER*16  ZONNAM
        CHARACTER*256 MESG  

C***********************************************************************
C   begin body of subroutine RDINVEN2

        IF( CATEGORY .EQ. 'AREA' ) THEN
            NVARI  = 3
            FIPNAM = 'FIP'
            VR2NAM = 'ASC7'
            VR3NAM = 'ASC3'
            ZONNAM = 'ZONES'

        ELSEIF( CATEGORY .EQ. 'MOBILE' ) THEN
            NVARI  = 3
            FIPNAM = 'IFIP'
            VR2NAM = 'IRCLAS'
            VR3NAM = 'ILINK'
            ZONNAM = 'TZONES'

        ELSEIF( CATEGORY .EQ. 'POINT' ) THEN
            NVARI  = 5
            FIPNAM = 'IFIP'
            VR2NAM = 'ISCC'
            VR3NAM = 'ISIC'
            VR4NAM = 'IPLANT'
            VR5NAM = 'ISTACK'
            ZONNAM = 'TZONES'

        ELSE  

            MESG = 'Emissions category "' // CATEGORY // '" not known.'
            CALL M3EXIT( 'RDINVEN2', 0, 0, MESG, 2 )

        ENDIF

        IF( NVARI .GE. 1 ) THEN
            IF( .NOT. READ3( FNAME, FIPNAM, ALLAYS3, 0, 0, IFIP ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading FIP from file ' //FNAME, 2 )
            END IF
        END IF

        IF( NVARI .GE. 2 ) THEN
            IF( .NOT. READ3( FNAME, VR2NAM, ALLAYS3, 0, 0, VAR2 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading ' //VR2NAM// ' from file '// 
     &                       FNAME, 2 )
            END IF
        ELSE
            DO 55 S = 1, NSRCS
                VAR2( S ) = 0
55          CONTINUE
        END IF

        IF( NVARI .GE. 3 ) THEN
            IF( .NOT. READ3( FNAME, VR3NAM, ALLAYS3, 0, 0, VAR3 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading ' //VR3NAM// ' from file '// 
     &                       FNAME, 2 )
            END IF
        ELSE
            DO 66 S = 1, NSRCS
                VAR3( S ) = 0
66          CONTINUE
        END IF

        IF( NVARI .GE. 4 ) THEN
            IF( .NOT. READ3( FNAME, VR4NAM, ALLAYS3, 0, 0, VAR4 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading ' //VR4NAM// ' from file '// 
     &                       FNAME, 2 )
            END IF
        ELSE
            DO 77 S = 1, NSRCS
                VAR4( S ) = 0
77          CONTINUE
        END IF

        IF( NVARI .GE. 5 ) THEN
            IF( .NOT. READ3( FNAME, VR5NAM, ALLAYS3, 0, 0, VAR5 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading ' //VR5NAM// ' from file '//
     &                       FNAME, 2 )
            END IF
        ELSE
            DO 88 S = 1, NSRCS
                VAR5( S ) = 0
88          CONTINUE
        END IF

C.........  General reads (common to all source categories)

        IF( .NOT. READ3( FNAME, ZONNAM, ALLAYS3, 0, 0, TZONES ) ) THEN
            CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                   'Error reading TZONES from file '// FNAME, 2 )
        END IF

        IF( .NOT. READ3( FNAME, 'TPFLAG', ALLAYS3, 0, 0, TPFLAG ) ) THEN
            CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                   'Error reading TPFLAG from file '// FNAME, 2 )
        END IF

        IF( .NOT. READ3( FNAME, 'INVYR', ALLAYS3, 0, 0, INVYR ) ) THEN
            CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                   'Error reading INVYR from file '// FNAME, 2 )
        END IF

        DO 101 V = 1, NIVARS
            IF( .NOT. READ3( FNAME, VNAMES( V ), ALLAYS3,
     &                           0, 0, INVVAL( 1,V )         ) ) THEN
                MESG = 'Error reading "' //
     &                 VNAMES( V )( 1:TRIMLEN( VNAMES( V ) ) ) //
     &                 '" from file "' //
     &                 FNAME( 1:TRIMLEN( FNAME ) ) // '".'
                CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )
            ENDIF
101     CONTINUE

C.........  Category-specific reads and adjustments

        IF( CATEGORY .EQ. 'MOBILE' ) THEN

            IF( .NOT. READ3( FNAME,'XLOC1',ALLAYS3,0,0,REAL1 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading XLOC1 from file '// 
     &                       FNAME, 2 )
            END IF

            IF( .NOT. READ3( FNAME,'YLOC1',ALLAYS3,0,0,REAL2 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading YLOC1 from file '// 
     &                       FNAME, 2 )
            END IF

            IF( .NOT. READ3( FNAME,'XLOC2',ALLAYS3,0,0,REAL3 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading XLOC2 from file '// 
     &                       FNAME, 2 )
            END IF

            IF( .NOT. READ3( FNAME,'YLOC2',ALLAYS3,0,0,REAL4 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading YLOC2 from file '// 
     &                       FNAME, 2 )
            END IF        

        ELSEIF( CATEGORY .EQ. 'POINT' ) THEN

            IF( .NOT. READ3( FNAME,'XLOCA',ALLAYS3,0,0,REAL1 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading XLOCA from file '// 
     &                       FNAME, 2 )
            END IF

            IF( .NOT. READ3( FNAME,'YLOCA',ALLAYS3,0,0,REAL2 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading YLOCA from file '// 
     &                       FNAME, 2 )
            END IF

            IF( .NOT. READ3( FNAME,'STKHT',ALLAYS3,0,0,REAL3 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading STKHT from file '// 
     &                       FNAME, 2 )
            END IF

            IF( .NOT. READ3( FNAME,'STKDM',ALLAYS3,0,0,REAL4 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading STKDM from file '// 
     &                       FNAME, 2 )
            END IF

            IF( .NOT. READ3( FNAME,'STKTK',ALLAYS3,0,0,REAL5 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading STKTK from file '// 
     &                       FNAME, 2 )
            END IF

            IF( .NOT. READ3( FNAME,'STKVE',ALLAYS3,0,0,REAL6 ) ) THEN
                CALL M3EXIT( 'RDINVEN2', 0, 0,
     &                       'Error reading STKVE from file '// 
     &                       FNAME, 2 )
            END IF

        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )
 
94020   FORMAT ( 10 ( A, :, E10.3 :, 2X ) )
 
        END
 
