
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  DAYREC( NDSRC, INDEXA, 
     &                      IFIPA, ISCCA, IPLANA, ISTACA, 
     &                      ICODEA, EMISDA, JDATE, JTIME,
     &                      ENAME, DNAME )

C***********************************************************************
C  subroutine body starts at line
C
C  FUNCTION:
C       Take records for JDATE,JTIME, sort them, construct subscript
C       table for them, and write them out to DAY-POINT file DNAME
C
C  PRECONDITIONS REQUIRED:
C       DNAME, ENAME already opened.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C
C  REVISION  HISTORY:
C       Prototype 10/95 by CJC for PDPOINT
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'    ! I/O API constants
        INCLUDE 'FDESC3.EXT'    ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem info (both inventory and model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  time related parameters


C...........   ARGUMENTS and their descriptions:

        INTEGER       NDSRC            !  number of sources
        INTEGER       INDEXA( NDSRC )  !  subscript table for SORTI4()
        INTEGER       IFIPA ( NDSRC )  !  source FIPS (county) ID
        INTEGER       ISCCA ( NDSRC )  !  source SCC
        INTEGER       IPLANA( NDSRC )  !  Plant ID
        INTEGER       ISTACA( NDSRC )  !  Stack ID
        INTEGER       ICODEA( NDSRC )  !  spcecies-subscript
        REAL          EMISDA( NDSRC )  !  emissions values (tons/yr)
        INTEGER       JDATE, JTIME     !  current date and time
        CHARACTER*16  ENAME     !  logical name for PNTS file from RAWPOINT
        CHARACTER*16  DNAME     !  " DAY-POINTS output file

C...........   PARAMETERS and their descriptions:


C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*24    DT2STR
        LOGICAL         ENVYN
        INTEGER         FIND1, FIND4   !  return -1 for failure
        INTEGER         TRIMLEN
        
        EXTERNAL DT2STR, ENVYN, FIND1, FIND4, TRIMLEN


C...........   SAVED LOCAL VARIABLES and their descriptions:

        LOGICAL         FIRSTIME
        DATA            FIRSTIME / .TRUE. /

        LOGICAL         VFLAG   !  input verification:  TRUE iff ERROR
        LOGICAL         WFLAG   !  input verification:  TRUE iff ERROR
        LOGICAL         SFLAG   !  input verification:  report missing species

C...........   RAWPOINT point Sources Table from ENAME

        INTEGER       IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER       ISCC  ( NPSRC )  !  source SCC
        INTEGER       IPLANT( NPSRC )  !  Plant ID
        INTEGER       ISTACK( NPSRC )  !  Stack ID

        SAVE          FIRSTIME, SFLAG, IFIP, ISCC, IPLANT, ISTACK

C...........   Day-specific point Sources Table input unsorted copy (*A); 
C...........   sorted final version (image EMISREC of data record)

        INTEGER       NPDPT           !  current source-count
        INTEGER       INDXD( NPDAY )   !  subscript into netCDF inventory file
        REAL          EMISD( NPDAY, NIPOL )  !  emissions values (tons/yr)

        COMMON / EMISREC / NPDPT, INDXD, EMISD

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         I, J, K, L, S, V    !  loop counters
        INTEGER         IOS 
        INTEGER         FIP
        INTEGER         SCC
        INTEGER         PLANT
        INTEGER         STACK
        REAL            EMIS
        CHARACTER*256   MESG    !  text for M3EXIT()

C***********************************************************************
C   begin body of subroutine  DAYREC

C.......   At setup:  get environment flag for missing record-checks;
C.......   read in FIP,SCC,PLANT,STACK tables from ENAME:

        IF ( FIRSTIME ) THEN

            FIRSTIME = .FALSE.

            SFLAG = ENVYN( 'RAW_SRC_CHECK', 
     &                     'PDPOINT check for missing species-records',
     &                     .FALSE.,
     &                     IOS )

            IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3, 0, 0, 
     &                        IFIP ) ) THEN
                CALL M3EXIT( 'PDPOINT:DAYREC', 0, 0, 
     &                       'Could not read IFIP from ' // ENAME, 2 )
            END IF

            IF ( .NOT. READ3( ENAME, 'ISCC', ALLAYS3, 0, 0, 
     &                        ISCC ) ) THEN
                CALL M3EXIT( 'PDPOINT:DAYREC', 0, 0, 
     &                       'Could not read ISCC from ' // ENAME, 2 )
            END IF

            IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3, 0, 0, 
     &                        IPLANT ) ) THEN
                CALL M3EXIT( 'PDPOINT:DAYREC', 0, 0, 
     &                       'Could not read IPLANT from ' // ENAME, 2 )
            END IF

            IF ( .NOT. READ3( ENAME, 'ISTACK', ALLAYS3, 0, 0, 
     &                        ISTACK ) ) THEN
                CALL M3EXIT( 'PDPOINT:DAYREC', 0, 0, 
     &                       'Could not read ISTACK from ' // ENAME, 2 )
            END IF

        END IF          !  if firstime



C.......   Use SORTI4() to perform an indirect sort by FIPS,SCC,PLANT,STACK
C.......   then permute the records according to the result:

        CALL SORTI4( NDSRC, INDEXA, IFIPA, ISCCA, IPLANA, ISTACA )

        DO  102  I = 1, NIPOL
        DO  101  S = 1, NPDAY
            EMISD( S,I ) = BADVAL3
101     CONTINUE
102     CONTINUE

        L  = -1
        J  =  0
        VFLAG = .FALSE.
        WFLAG = .FALSE.

        DO  133  S = 1, NDSRC

            I = INDEXA( S )

            FIP   = IFIPA ( I )
            SCC   = ISCCA ( I )
            PLANT = IPLANA( I )
            STACK = ISTACA( I )
            
            K = FIND4( FIP, SCC, PLANT, STACK, 
     &                 NPSRC, IFIP, ISCC, IPLANT, ISTACK )

            V    = ICODEA( I )

            IF ( K .LT. 0 ) THEN

                WFLAG = .TRUE.
                WRITE( MESG,94020 )
     &              'Invalid source -- FIP:', FIP,
     &              'SCC:',     SCC, 
     &              'Plant:',   PLANT,
     &              'Stack:',   STACK,
     &              'Species:', PCODE( V )
                GO TO  133

            ELSE IF ( L .NE. K ) THEN

                L = K
                J = J + 1
                INDXD( J ) = K

            END IF      !  if first encounter with this source

            EMIS = EMISDA( I )
            IF ( EMISD( J,V ) .LT. AMISS3 ) THEN
                EMISD ( J,V ) = EMIS
            ELSE
                VFLAG = .TRUE.
                WRITE( MESG,94020 )
     &              'Duplicate emissions record:  FIP:', FIP,
     &              'SCC:',     SCC, 
     &              'Plant:',   PLANT,
     &              'Stack:',   STACK,
     &              'Species:', PCODE( V )
                CALL M3MESG( MESG )
            END IF

133     CONTINUE

        NPDPT = J

        IF ( WFLAG ) THEN

            CALL M3EXIT( 'PDPOINT:DAYREC', 0, 0, 
     &                   'Invalid sources found in input file', 2 )

        ELSE IF ( VFLAG ) THEN

            CALL M3EXIT( 'PDPOINT:DAYREC', 0, 0, 
     &                   'Duplicate records found in input file', 2 )

        END IF

        IF ( SFLAG ) THEN       !  if missing records are fatal:

            VFLAG = .FALSE.

            DO  146  V = 1, NIPOL
            DO  145  S = 1, NPDPT

               IF ( EMISD( S,V ) .LT. AMISS3 ) THEN

                    VFLAG = .TRUE.
                    WRITE( MESG,94020 )
     &                  'Missing emissions record:  FIP:', IFIP( S ),
     &                  'SCC:',     ISCC( S ), 
     &                  'Plant:',   IPLANT( S ),
     &                  'Stack:',   ISTACK( S ),
     &                  'Species:', PCODE( V )
                    CALL M3MESG( MESG )

                END IF

145         CONTINUE
146         CONTINUE

            IF ( VFLAG ) THEN
                CALL M3EXIT( 'PDPOINT:DAYREC', 0, 0,
     &                       'Missing records in input file', 2 )
            END IF

        ELSE            !  "missing" ~~~> 0

            DO  148  V = 1, NIPOL
            DO  147  S = 1, NPDAY
                IF ( EMISD( S,V ) .LT. AMISS3 ) EMISD( S,V ) = 0.0
147         CONTINUE
148         CONTINUE

        END IF          !  sflag or not:  check for non-set EMIS values


        IF ( .NOT. WRITE3( DNAME, ALLVAR3, JDATE, JTIME, 
     &                     NPDPT ) ) THEN
            CALL M3EXIT( 'PDPOINT/DAYREC', JDATE, JTIME,
     &           'Error writing output file "' //
     &           DNAME( 1:TRIMLEN( DNAME ) ) // '"',  2 )
        END IF
     
      RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx

91000   FORMAT ( //5X , '*** ERROR ABORT in subroutine DAYREC ***',
     &            /5X , A ,
     &           // )        !  generic error message format


C...........   Informational (LOG) message formats... 92xxx


C...........   Formatted file I/O formats............ 93xxx


C...........   Internal buffering formats............ 94xxx

94020   FORMAT( A, X, I5, X, A, I7, X, A, X, I5, X, A, X, A, I3, A, I5 )


C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

