C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File Version: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C * All Rights Reserved
C *
C * See file COPYRIGHT for conditions of use.
C *
C * Environmental Programs Group
C * MCNC--North Carolina Supercomputing Center
C * P.O. Box 12889
C * Research Triangle Park, NC  27709-2889
C *
C * env_progs@mcnc.org
C *
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C
****************************************************************************/

        PROGRAM SFCINTERP
 
C***********************************************************************
C  program body starts at line 148 
C
C  DESCRIPTION:
C       This program interpolates MM5 output (MCIP file) to desired heights, 
C       and writes the result to a netCDF file. Heights and Output Variables 
C       can be defined by user using the script file
C       - up to 99 different interpolation heights 
C       - var's can be chosen from WIND, TEMP, Mixing Ratio, Pressure, RH,
C         UWIND, and VWIND. 
C       - FDESC3D needs to be used to specify scenario and scheme for the
C         subsequent SMOKE programs to read 
C
C         NOTE: 
C         1.  Wind's U, V components (UWIND and VWIND) are only approximations 
C             using the direction of the first level wind. If the interpolation
C             height is beyond (10m,30m) above sfc, the result will be dubious. 
C         2.  The output variable's name has height information embedded and 
C             is rounded to a integer (PAVE doesn't take "point" in a variable 
C             name . Therefore, interpolations to different heights but with 
C             the same integer before the floating point should be run 
C             separately to prevent the construction of the same output 
C             variable name.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       PBLPWR and DIRECT     
C       by Daewon Byun and Kiran Alapaty
C
C  REVISION  HISTORY:
C       Prototype  7/96 by C. SUNNY LO
C       Modified  12/96 by C. SUNNY LO
C
C***********************************************************************
 
       IMPLICIT NONE
 
       INCLUDE 'GRDIMS3.EXT'   !  grid parameters
       INCLUDE 'PARMS3.EXT'    !  I/O API parameters
       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
       INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
       INCLUDE 'PARAMS.EXT'    !  Parameters.
 
C........  Enternal functions
 
       LOGICAL         DSCGRID
       INTEGER         PROMPTFFILE
       CHARACTER*16    PROMPTMFILE
       REAL            GETREAL
       LOGICAL         GETYN
 
       EXTERNAL        DSCGRID, PROMPTFFILE,
     &                 PROMPTMFILE, GETREAL, GETYN
 
C........  Local variables
       INTEGER         LDEV      !  UNIT NUMBER FOR LOG FILE 
       CHARACTER*16    MC3NAME   !  logical name for MET_CRO_3D
       CHARACTER*16    MD3NAME   !  logical name for MET_DOT_3D
       CHARACTER*16    GC3NAME   !  logical name for GRID_CRO_3D
       CHARACTER*16    MC2NAME   !  logical name for MET_CRO_2D
       CHARACTER*16    GC2NAME   !  logical name for GRID_CRO_2D
       CHARACTER*16    ONAME     !  logical name for MET_OUT
       CHARACTER*256   MESG      !  message buffer 
       CHARACTER*16    GNAME     !  scratch copy of grid name
       CHARACTER*16    CNAME     !  coordinate system name from DSCGRID()
       CHARACTER*16    STRINGA   !  Intermediate Character String
       CHARACTER*256   STRINGB   !  Intermediate Character String
       CHARACTER*16    TEMPNAME  !  Intermediate VNAME3D
       CHARACTER*256   TEMPDESC  !  Intermediate VDESC3D
       CHARACTER*16    TEMPUNIT  !  Intermediate UNITS3D
       INTEGER         INT_HGT   !  Integer height in I5.5 format
       CHARACTER*5     CHAR_HGT  !  Integer height in character format
       INTEGER         ITEMOUT   !  Total items to be written out
       LOGICAL         YN          

       INTEGER         NTP       !  Iterator for interpolation height
       INTEGER         INDEX     !  Iterator
       INTEGER         TOTALHTS  !  Total heights to be interpolated
       INTEGER         I, J, K, T, C, JDATE, JTIME, TSTEP
       REAL            XPBL, XZ0, XTG, XPS  !  Intermediate Variable
       REAL            BUF3D( NCOLS+1, NROWS+1, NL, 2 ) ! Input Data Buffer
C                                1. UWIND 
C                                2. VWIND
       REAL            BUF3C( NCOLS, NROWS, NL, 5 ) ! Input Data Buffer
C                                1. TA
C                                2. QV
C                                3. PRES
C                                4. ZH
C                                5. ZF
       REAL            BUF2C( NCOLS, NROWS, 5 )     ! Input Data Buffer
C                                1. PBL
C                                2. MOL
C                                3. ZZERO
C                                4. TGD
C                                5. PSTAR0

       LOGICAL         SELECT(OUTVAR, NTPHTMAX)
       CHARACTER*5     VARNAME ( OUTVAR )
       DATA            VARNAME / 'WIND', 'TEMP', 'QV', 'PRES', 'RH',
     &                           'UWIND', 'VWIND',
     &                           'U_CRO', 'V_CRO' /
       REAL            PI180
       PARAMETER       (PI180 = 3.141592654/180.0)
C........   Intermediate variables

       REAL            PSFC( NCOLS, NROWS ) ! Surface Pressure calculated from
C                                             MCIP output 


C........   Arguments for PBLPWR

       REAL            UL(NL),TL(NL),QL(NL),
     &                 PL(NL),ZL(NL),ZF(0:NL)
       REAL            UUL,VVL,WINDDIR
       REAL            AMOL 
       REAL            UHGT, THGT, QHGT, PHGT, RHGT
       REAL            HGT(NTPHTMAX)

C........   Output array for IO/API

       REAL            VARHGT  ( NCOLS, NROWS, OUTVAR , NTPHTMAX)
C                                1. Wind Speed
C                                2. Temp 
C                                3. Mixing Ratio
C                                4. Pressure
C                                5. Relative Humidity
C                                6. U Wind component
C                                7. V Wind component
C                                8. U Wind component on CROSS point
C                                9. V Wind component on CROSS point

C***********************************************************************
C   begin body of program SFCINTERP
     
       LDEV =INIT3()

C........  Initialize Selection List

       DO C = 1, OUTVAR
          DO I = 1, NTPHTMAX
             SELECT(C, I) = .FALSE.
             HGT( I ) = -9999.0
          END DO
       END DO


C........  Process User Input

       ITEMOUT = 0        ! Total item number that will be written to output

       I = 0
10     CONTINUE
       I =  I + 1
       HGT(I) = GETREAL( 0.0, 16000.00, 1.5, 
     &          'Enter a Height (M) between 0.0 and 16000.0' )

C........  Select Variables needed at this height

       DO 30 C = 1, OUTVAR
          WRITE(MESG, 5000) VARNAME(C)
5000      FORMAT('WANT ',A5 ' AT THIS HEIGHT ?')
          SELECT(C, I) = GETYN ( MESG, .TRUE. )
          IF ( SELECT(C, I) ) ITEMOUT = ITEMOUT + 1
30     END DO

       YN = GETYN('Want Interpolation at Another Height?', 
     &            .TRUE.)

       IF ( YN ) THEN
          GO TO 10
       ELSE
          GO TO 50
       END IF


C........  Process Input File (MCIP output) 

50     CONTINUE 
       TOTALHTS = I

       MC3NAME = PROMPTMFILE(
     &        'Enter logical name for MET_CRO_3D file',
     &         FSREAD3, 'MC3NAME', 'SFCINTERP' )
 
       MD3NAME = PROMPTMFILE(
     &        'Enter logical name for MET_DOT_3D file',
     &         FSREAD3, 'MD3NAME', 'SFCINTERP' )
 
       GC3NAME = PROMPTMFILE(
     &        'Enter logical name for MET input file',
     &         FSREAD3, 'GC3NAME', 'SFCINTERP' )
 
       MC2NAME = PROMPTMFILE(
     &        'Enter logical name for MET input file',
     &         FSREAD3, 'MC2NAME', 'SFCINTERP' )
 
       GC2NAME = PROMPTMFILE(
     &        'Enter logical name for MET input file',
     &         FSREAD3, 'GC2NAME', 'SFCINTERP' )

       IF ( .NOT. DESC3( MC3NAME ) ) THEN
           MESG = 'Could not get description of file  ' //
     &            MC3NAME
           CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
       END IF
 
       IF ( .NOT. DSCGRID( GRDNM, CNAME, GDTYP3D,
     &            P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &            XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &            NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
 
 
          GNAME = GRDNM
          MESG  = 'COULD NOT GET DESCRIPTION OF GRID'
          CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
 
       END IF          !  IF DSCGRID() FAILED

       JDATE = SDATE3D
       JTIME = STIME3D
       TSTEP = 10000

       GDNAM3D = GRDNM
       FTYPE3D = GRDDED3
c      SDATE3D = 1995188       !  n/a
c      STIME3D = 130000        !  n/a
c      TSTEP3D = 10000         !  time independent
       NVARS3D = ITEMOUT
       NLAYS3D = 1
       NTHIK3D = 1
       VGTYP3D = IMISS3
       VGTOP3D = AMISS3
       FDESC3D( 1 ) = 'MM5 Variables Interpolated to a Specific Height'

       DO K = 2, MXDESC3
          FDESC3D( K ) = ' '
       END DO 

       FDESC3D( 2 ) = 'Scenario: SMRAQ95'
       FDESC3D( 3 ) = 'Scheme:   KUO' 
       
       INDEX = 0
       INT_HGT = 0
       DO J = 1, TOTALHTS

          INT_HGT = INT( HGT(J) )
          WRITE(CHAR_HGT, 9990)  INT_HGT 

          DO I = 1, OUTVAR

             IF ( SELECT(I,J) ) THEN
                IF ( I.EQ.1 ) THEN
                   WRITE(TEMPNAME,10000) 'WIND_', CHAR_HGT
                   WRITE(TEMPDESC,10010) 'WIND at ', HGT(J), ' m'
                   TEMPUNIT = 'M/S'
                ELSE IF ( I.EQ.2 ) THEN
                   WRITE(TEMPNAME,10000) 'TEMP_', CHAR_HGT 
                   WRITE(TEMPDESC,10010) 'TEMP at ', HGT(J), ' m'
                   TEMPUNIT = 'K'
                ELSE IF ( I.EQ.3 ) THEN
                   WRITE(TEMPNAME,10000) 'QV_', CHAR_HGT
                   WRITE(TEMPDESC,10010) 'QV at ', HGT(J), ' m'
                   TEMPUNIT = 'G/G'
                ELSE IF ( I.EQ.4 ) THEN
                   WRITE(TEMPNAME,10000) 'PRES_', CHAR_HGT
                   WRITE(TEMPDESC,10010) 'PRES at ', HGT(J), ' m'
                   TEMPUNIT = 'MB'
                ELSE IF ( I.EQ.5 ) THEN
                   WRITE(TEMPNAME,10000) 'RH_', CHAR_HGT
                   WRITE(TEMPDESC,10010) 'RH at ', HGT(J), ' m'
                   TEMPUNIT = '-'
                ELSE IF ( I.EQ.6 ) THEN
                   WRITE(TEMPNAME,10000) 'UWIND_', CHAR_HGT
                   WRITE(TEMPDESC,10010) 'UWIND at ', HGT(J), ' m'
                   TEMPUNIT = 'M/S'
                ELSE IF ( I.EQ.7 ) THEN
                   WRITE(TEMPNAME,10000) 'VWIND_', CHAR_HGT
                   WRITE(TEMPDESC,10010) 'VWIND at ', HGT(J), ' m'
                   TEMPUNIT = 'M/S'
                ELSE IF ( I.EQ.8 ) THEN
                   WRITE(TEMPNAME,10000) 'U_CRO_', CHAR_HGT
                   WRITE(TEMPDESC,10010) 'U_CRO at ', HGT(J), ' m'
                   TEMPUNIT = 'M/S'
                ELSE IF ( I.EQ.9 ) THEN
                   WRITE(TEMPNAME,10000) 'V_CRO_', CHAR_HGT
                   WRITE(TEMPDESC,10010) 'V_CRO at ', HGT(J), ' m'
                   TEMPUNIT = 'M/S'
                ELSE
                   WRITE(TEMPNAME,10000) 'N/A_', CHAR_HGT
                   TEMPDESC = 'N/A'
                   TEMPUNIT = 'N/A'
                END IF

                INDEX = INDEX + 1 
                VNAME3D( INDEX ) = TEMPNAME
                VDESC3D( INDEX ) = TEMPDESC
                UNITS3D( INDEX ) = TEMPUNIT
                VTYPE3D( INDEX ) = M3REAL
             END IF     ! SELECT = .TRUE.

          END DO
       END DO   

C 9980   FORMAT( I5 )
9990   FORMAT( I5.5 )
10000  FORMAT( A, A5 ) 
10010  FORMAT( A, F8.2, A )

       ONAME = PROMPTMFILE(
     &        'Enter logical name for MET output file',
     &         FSUNKN3, 'MET_OUT', 'SFCINTERP' )

C........    Read Variables one timestep at a time

       DO 300 T = 1, MXREC3D

             IF ( .NOT. READ3( MC3NAME, 'TA', ALLAYS3, JDATE, 
     &            JTIME, BUF3C(1,1,1,1) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF

             IF ( .NOT. READ3( MC3NAME, 'QV', ALLAYS3, JDATE, 
     &            JTIME, BUF3C(1,1,1,2) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF

             IF ( .NOT. READ3( MC3NAME, 'PRES', ALLAYS3, JDATE, 
     &            JTIME, BUF3C(1,1,1,3) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF

             IF ( .NOT. READ3( MD3NAME, 'UWIND', ALLAYS3, JDATE, 
     &            JTIME, BUF3D(1,1,1,1) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF

             IF ( .NOT. READ3( MD3NAME, 'VWIND', ALLAYS3, JDATE, 
     &            JTIME, BUF3D(1,1,1,2) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF

             IF ( .NOT. READ3( GC3NAME, 'ZH', ALLAYS3, JDATE, 
     &            JTIME, BUF3C(1,1,1,4) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF

             IF ( .NOT. READ3( GC3NAME, 'ZF', ALLAYS3, JDATE, 
     &            JTIME, BUF3C(1,1,1,5) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF

             IF ( .NOT. READ3( GC2NAME, 'ZZERO', ALLAYS3, JDATE, 
     &            JTIME, BUF2C(1,1,3) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF

             IF ( .NOT. READ3( MC2NAME, 'PBL', ALLAYS3, JDATE, 
     &            JTIME, BUF2C(1,1,1) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF

             IF ( .NOT. READ3( MC2NAME, 'MOL', ALLAYS3, JDATE, 
     &            JTIME, BUF2C(1,1,2) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF

             IF ( .NOT. READ3( MC2NAME, 'TGD', ALLAYS3, JDATE, 
     &            JTIME, BUF2C(1,1,4) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF
             
             IF ( .NOT. READ3( GC2NAME, 'PSTAR0', ALLAYS3, JDATE, 
     &            JTIME, BUF2C(1,1,5) ) ) THEN
                  MESG = 'Could not read from file'
                 CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
             END IF


C........    Prepare Variables to pass into PBLPWR, and Call PBLPWR
C........    to do the interpolation for each grid - one grid at a time

             DO 90 I = 1, NCOLS
                DO 80 J = 1, NROWS

                   PSFC(I,J) = 0.01 * ( 0.0025 * BUF2C( I,J,5 )
     &                           + BUF3C( I, J, 1, 3 ) ) ! Calc. SFC Pres.

                   DO 70 K = 1, NL
                      UUL 
     &                = 0.25*( BUF3D(I,J,K,1)+BUF3D(I+1,J,K,1)
     &                       + BUF3D(I,J+1,K,1)+BUF3D(I+1,J+1,K,1) ) 
                      VVL 
     &                = 0.25*( BUF3D(I,J,K,2)+BUF3D(I+1,J,K,2)
     &                       + BUF3D(I,J+1,K,2)+BUF3D(I+1,J+1,K,2) )  

                      UL(K)            ! Calculate Wind at Cross Position
     &                = SQRT( UUL**2 + VVL**2 )

                      IF ( K.EQ.1 ) THEN ! Calculate Wind direction at 1st layer
                         VARHGT(I,J,8,2) = UUL  ! hardwired for test
                         VARHGT(I,J,9,2) = VVL
                         CALL DIRECT( UUL,VVL,WINDDIR )
                      END IF 

                      TL ( K ) = BUF3C ( I, J, K, 1 )   
                      QL ( K ) = BUF3C ( I, J, K, 2 )   
                      PL ( K ) = 0.01 * BUF3C ( I, J, K, 3 )   
                      ZL ( K ) = BUF3C ( I, J, K, 4 )   
                      ZF ( K ) = BUF3C ( I, J, K, 5 )   

70                 END DO

                   ZF(0) = 0.0
                   AMOL = 1.0 / BUF2C(I,J,2)
                   XPBL = BUF2C(I,J,1)
                   XZ0  = BUF2C(I,J,3)
                   XTG  = BUF2C(I,J,4)
                   XPS  = PSFC(I,J)

                   DO 75 NTP = 1, TOTALHTS

                   CALL PBLPWR(UL,TL,QL,PL,ZL,ZF,NL, HGT(NTP),
     &                      XPBL, AMOL, XZ0, XTG, XPS,
     &                      UHGT,THGT,QHGT,PHGT,RHGT)

                   VARHGT(I,J,1,NTP) = UHGT
                   VARHGT(I,J,2,NTP) = THGT
                   VARHGT(I,J,3,NTP) = QHGT
                   VARHGT(I,J,4,NTP) = PHGT
                   VARHGT(I,J,5,NTP) = RHGT
                   VARHGT(I,J,6,NTP) = -1.0*UHGT*SIN(WINDDIR*PI180)
                   VARHGT(I,J,7,NTP) = -1.0*UHGT*COS(WINDDIR*PI180)

75                 END DO
80              END DO
90           END DO


C........    Write the interpolated field into netCDF file

             C = 0
             DO 110 NTP = 1, TOTALHTS
             DO 100 I = 1, OUTVAR

                IF ( SELECT(I,NTP) ) THEN
                   C = C + 1
                   IF ( .NOT. WRITE3( ONAME, VNAME3D(C), 
     &                  JDATE, JTIME, VARHGT(1,1,I,NTP) ) ) THEN
                      MESG = 'Could not write ' // VNAME3D(C)
                      CALL M3EXIT( 'SFCINTERP', 0, 0, MESG, 2 )
                   END IF
                END IF

100             END DO
110          END DO

C........    Finished processing one time step, call nextime

             CALL NEXTIME( JDATE, JTIME, TSTEP )

300    END DO    ! T = 1, MXREC3D

       CALL M3EXIT( 'SFCINTERP', 0, 0,
     &              'Normal completion of PROGRAM SFCINTERP', 0 )

       STOP
       END
