 
        INTEGER FUNCTION CSVNAME( IDEV  , LINVEN, LMODEL, LPROCESS,
     &                            VNAMES )

C***********************************************************************
C  function body starts at line 116
C
C  DESCRIPTION:
C       Creates variable names from combinations of inventory species,
C       model species, and MV process names.  It reads the
C       file SPCLIST, checks for consistency, checks input parameters
C       to function, and outputs variable names based on the TYPE
C       specified.  
C       ******* Returned names are in compiled order! *******
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Prototype  5/96 by MRH.
C       5/31/96 MRH changed to allow names only in order allowed by EXT file
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
 
        INCLUDE 'CHDGEN3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'MBDIMS3.EXT'   !  motor vehicles parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER       INDEX1
        INTEGER       TRIMLEN

        EXTERNAL      INDEX1, TRIMLEN

C...........   FUNCTION PARAMETERS and their descriptions:
        INTEGER       IDEV      ! SPCLIST unit number
        LOGICAL       LINVEN    ! Indicates if the name has inventory species
        LOGICAL       LMODEL    ! Indicates if the name has model species 
        LOGICAL       LPROCESS  ! Indicates if the name has MV process
        CHARACTER*(*) VNAMES( * )
                                ! Array of the names of the variables returned

C...........   LOCAL VARIABLES and their descriptions:
        INTEGER       INVIDX( MXIPOL )   ! Pointer to indiv. inventory species 
        INTEGER       NUMINVA            ! Number of actual inventory species
        INTEGER       NUMMODA            ! Number of actual model species
        INTEGER       NUMAMVP            ! Number of actual MV processes
        INTEGER       NP                 ! Temporary number of MV processes
        INTEGER       NUMPC  ( MXMPOL )  ! No. MV processes per model species
        INTEGER       T_NUMPC( MXMPOL )  ! Unsorted no. MV proc per model spcs
    
        CHARACTER*3   ACTINV ( MXIPOL )  ! Actual inventory species list
        CHARACTER*3   ACTMVP ( MXEPRC )  ! Actual MV process list
        CHARACTER*3   EM     ( MXEPRC )  ! Temporary MV process names
        CHARACTER*3   INV    ( MXMPOL )  ! Inv species list per model species
        CHARACTER*3   ISPEC              ! Temporary inventory species names
        CHARACTER*3   PROC   ( MXMPOL, MXEPRC ) ! Actual MV procs per model spcs
        CHARACTER*3   T_ACTMVP( MXEPRC ) ! Actual MV process list
        CHARACTER*3   T_INV  ( MXMPOL )  ! Unsorted inventory species array
        CHARACTER*3   T_PROC ( MXEPRC, MXMPOL ) ! Unsort MV procs per model spcs

        CHARACTER*8   MSPEC              ! Temporary model species names
        CHARACTER*8   MODEL  ( MXMPOL )  ! Actual model species list
        CHARACTER*8   T_MODEL( MXMPOL )  ! Unsorted model species array

        SAVE          ACTINV, ACTMVP,
     &                INV, INVIDX, MODEL, NUMINVA, NUMMODA, 
     &                NUMAMVP, NUMPC, PROC

        LOGICAL       FIRSTIME           ! Indicates if first function call
        DATA          FIRSTIME / .TRUE. /
        SAVE          FIRSTIME

C...........   Other local variables

        INTEGER       COUNT, C1, C2
        INTEGER       I, ICNT, II, IOS, J, K, L , LEN
        LOGICAL       EFLAG

        CHARACTER*256 MESG      ! Error mesg returned if exit status ne 0
        CHARACTER*50  NAME

C***********************************************************************
C   begin body of function CSVNAM

C.........  Read SPCLIST file from unit number IDEV if firstime function called

        IF( FIRSTIME ) THEN

            FIRSTIME = .FALSE.
            CSVNAME = 0

C.............  Read past header
            READ( IDEV, * ) 

C.............  Read table of name combinations
            EFLAG = .FALSE.
            ICNT  = 1 ! Already past header
            COUNT = 0
   33       CONTINUE  ! Head of table read

               ICNT = ICNT + 1
               READ( IDEV, *, END=44, IOSTAT= IOS ) 
     &               ISPEC, MSPEC, NP, ( EM( I ), I=1,NP )

               IF( IOS .NE. 0 ) THEN
                   WRITE( MESG,94010 )
     &                     'Error ', IOS,
     &                     'reading SPECIES LIST file at line ', ICNT
                   CALL M3MESG( MESG )
                   CSVNAME = IMISS3
                   RETURN

               ELSEIF( COUNT .GE. MXMPOL ) THEN
                   WRITE( MESG,94010 )
     &                  'Max number of model species', MXMPOL,
     &                  'exceeded by SPECIES LIST file at line', ICNT
                   CALL M3MESG( MESG )
                   CSVNAME = IMISS3
                   RETURN

               ELSEIF( NP .GT. MXEPRC ) THEN
                   WRITE( MESG,94010 )
     &                  'Max number of MV processes', MXEPRC,
     &                  'exceeded by SPECIES LIST file at line', ICNT
                   CALL M3MESG( MESG )
                   CSVNAME = IMISS3
                   RETURN

               ELSE

C....................  Look for species names in valid lists
                   K = INDEX1( ISPEC, MXIPOL, INVSPEC )
                   L = INDEX1( MSPEC, MXMPOL, MODSPEC )
 
C....................  Reports error if names not in lists
                   IF( K .EQ. 0 ) THEN
                       EFLAG = .TRUE.
                       WRITE( MESG,94010 )
     &                        'Invalid inventory species "' //
     &                        ISPEC( 1:TRIMLEN( ISPEC ) ) //
     &                        '" in SPECIES LIST file at line', ICNT
                       CALL M3MESG( MESG )
 
                   ELSEIF( L .EQ. 0 ) THEN
                       EFLAG = .TRUE.
                       WRITE( MESG,94010 )
     &                        'Invalid model species "' //
     &                        MSPEC( 1:TRIMLEN( MSPEC ) ) //
     &                        '" in SPECIES LIST file at line', ICNT
                       CALL M3MESG( MESG )
 
                   ENDIF

C....................  Store table of names in arrays
                   COUNT = COUNT + 1 
                   T_INV  ( COUNT ) = ISPEC
                   T_MODEL( COUNT ) = MSPEC
                   T_NUMPC( COUNT ) = NP

C....................  Loop through all MV process for storage and checks
                   DO 37 J = 1,NP
 
                       K = INDEX1( EM( J ), MXEPRC, MVPROCS )
 
C........................  Reports error if MV process not in list
                       IF( K .EQ. 0 ) THEN
                           EFLAG = .TRUE.
                           WRITE( MESG,94010 )
     &                            'Invalid MV process "' // EM( J ) //
     &                            EM(J)( 1:TRIMLEN( EM(J) ) ) //
     &                            '" in SPECIES LIST file at line', ICNT
                           CALL M3MESG( MESG )

C........................  Store MV process names
                       ELSE
                           T_PROC( J, COUNT ) = EM( J )
 
                       ENDIF

   37              CONTINUE

               ENDIF

            GO TO 33  ! End of table read

   44       CONTINUE
            NUMMODA = COUNT

C.............  Abort if found error in list
            IF( EFLAG ) THEN
                CSVNAME = IMISS3
                RETURN
            ENDIF

C.............  Sort the input list based on compiled order of model species
C.............  Do this by searching for compiled names in temp list

            C1 = 0
            DO 49 I = 1, MXMPOL

                K = INDEX1( MODSPEC( I ), NUMMODA, T_MODEL )
                IF( K .GT. 0 ) THEN
                    C1 = C1 + 1
                    MODEL( C1 ) = MODSPEC( I )
                    INV  ( C1 ) = T_INV  ( K )
                    NUMPC( C1 ) = T_NUMPC( K )

                    C2 = 0
                    DO 47 J = 1, MXEPRC

                       L = INDEX1( MVPROCS( J ), NUMPC( C1 ), 
     &                             T_PROC( 1,K )             )
                       IF( L .GT. 0 ) THEN
                           C2 = C2 + 1
                           PROC( C1,C2 ) = T_PROC( L,K )
                       ENDIF

   47               CONTINUE

                ENDIF

   49       CONTINUE

C.............  Create index of actual MV processes
            COUNT = 1
            T_ACTMVP( 1 ) = PROC( 1,1 )
            DO 55 I = 1, NUMMODA
                DO 51 J = 1, NUMPC( I )
 
                    K = INDEX1( PROC( I,J ), COUNT, T_ACTMVP )
                    IF( K .EQ. 0 ) THEN
                        COUNT = COUNT + 1
                        T_ACTMVP( COUNT ) = PROC( I,J )
                    ENDIF
 
   51           CONTINUE
   55       CONTINUE
            NUMAMVP = COUNT

C.............  Sort actual MV procs using compiled list
            COUNT = 0
            DO 56 I = 1, MXEPRC
 
                K = INDEX1( MVPROCS( I ), NUMAMVP, T_ACTMVP )
                IF( K .GT. 0 ) THEN
                    COUNT = COUNT + 1
                    ACTMVP( COUNT ) = MVPROCS( I )
                ENDIF
 
   56       CONTINUE

C.............  Create array of actual inventory species and index to table
            COUNT = 1
            ACTINV( 1 ) = INV( 1 )
            INVIDX( 1 ) = 1
            DO 66 I = 1, NUMMODA

                K = INDEX1( INV( I ), COUNT, ACTINV )
                IF( K .EQ. 0 ) THEN
                    COUNT = COUNT + 1
                    ACTINV( COUNT ) = INV( I )
                    INVIDX( COUNT ) = I
                ENDIF

   66       CONTINUE
            NUMINVA = COUNT

        ENDIF ! firstime

C.........  Initialize function return status
        CSVNAME = 0

C.........  Names have inventory species...
        IF ( LMODEL ) THEN

            COUNT = 0
            DO 111 I = 1, NUMMODA

C.................   Names have MV process...
                IF( LPROCESS ) THEN

                    DO 106 J = 1, NUMPC( I )

                        COUNT = COUNT + 1 

C.........................   Names have MV process, inv., and model species
                        IF( LINVEN ) THEN
                            VNAMES( COUNT )= 
     &                       PROC (I,J)(1:TRIMLEN( PROC(I,J))) // '_' // 
     &                       INV  ( I )(1:TRIMLEN( INV ( I ))) // '_' // 
     &                       MODEL( I )

C.........................   Names have MV process and model species
                        ELSE
                            VNAMES( COUNT )= 
     &                       PROC (I,J)(1:TRIMLEN( PROC(I,J))) // '_' //
     &                       MODEL( I )
                        ENDIF

  106               CONTINUE

C.................   Names have inventory & model species, but no MV process
                ELSEIF( LINVEN ) THEN
                    COUNT = COUNT + 1 
                    VNAMES( COUNT ) = 
     &                       INV  ( I )(1:TRIMLEN( INV ( I ))) // '_' //
     &                       MODEL( I )

C.................   Names are just model species
                ELSE                   
                    COUNT = COUNT + 1 
                    VNAMES( COUNT ) =  MODEL( I )

                ENDIF

  111       CONTINUE

C.........  Names do not have inventory species...
        ELSE

C.............  Names have MV process and inventory species
            IF( LPROCESS .AND. LINVEN ) THEN

                COUNT = 0
                DO 122 I = 1, NUMINVA
                    II = INVIDX( I )

                    DO 116 J = 1, NUMPC( II )

                        COUNT = COUNT + 1 
                        VNAMES( COUNT )=
     &                     PROC (II,J)(1:TRIMLEN(PROC(II,J)))// '_' //
     &                     ACTINV( I )

C.........................  Check to make sure valid
                        K = INDEX1( VNAMES( COUNT ), MXEMIS, EMTLST )

                        IF( K .EQ. 0 ) THEN
                            LEN = TRIMLEN( VNAMES( COUNT ) )
                            NAME = VNAMES( COUNT )( 1:LEN )
                            WRITE( MESG, 94010 ) 
     &                          'Invalid emission type name: ' // NAME
                            CALL M3MESG( MESG )
                        ENDIF

  116               CONTINUE

  122           CONTINUE

C.............  Names have just MV process
            ELSEIF( LPROCESS ) THEN

                DO 133 I = 1, NUMAMVP
                    VNAMES( I ) = ACTMVP( I )
  133           CONTINUE
                COUNT = NUMAMVP

C.............  Names have just inventory species
            ELSEIF( LINVEN ) THEN

                DO 144 I = 1, NUMINVA
                    VNAMES( I ) = ACTINV( I )
  144           CONTINUE
                COUNT = NUMINVA

            ENDIF

        ENDIF  ! Use model species name or not

        CSVNAME = COUNT

        RETURN

C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Informational (LOG) message formats... 92xxx

C...........   Formatted file I/O formats............ 93xxx

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

        END
