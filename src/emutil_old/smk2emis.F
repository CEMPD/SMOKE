
C Version @(#)$Id$ $Source$ $Date$  
C
      PROGRAM SMK2EMIS

C***********************************************************************
C  program body starts at line 90
C
C  DESCRIPTION:
C       Use NetCDF input and EPS subroutines to produce UAM input files
C       Originally developed for EMISSIONS input, but potential use for
C       other UAM input files.
C
C  PRECONDITIONS REQUIRED:
C       M3IO single layer gridded data
C       EPS USERIN file relevant to model run for output data
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       GETYN, PROMPTFFILE, PROMPTMFILE, NEXTIME
C       EPS: RDCNTL, RDREGN, TRANS
C       CHKVAR, CHKTMP, CHKGRD
C
C  REVISION  HISTORY:
C       Prototype  9/95 by MRH
C
C***********************************************************************

        IMPLICIT NONE

C.........  INCLUDES:
        INCLUDE 'GRDIMS3.EXT'
        INCLUDE 'PARMS3.EXT'
        INCLUDE 'IODECL3.EXT'
        INCLUDE 'FDESC3.EXT'

C.........  UAM includes (in /200p_scratch/uamguides/eps2.0b/src/include)
        INCLUDE 'param.inc'
        INCLUDE 'uamcntl.inc'
        INCLUDE 'uamregn.inc'
        INCLUDE 'mrguam.inc' 

C.........  EXTERNAL FUNCTIONS and their descriptions:
        INTEGER      PROMPTFFILE
        CHARACTER*16 PROMPTMFILE
        LOGICAL      GETYN
 
        EXTERNAL     GETYN, PROMPTFFILE, PROMPTMFILE

C.........  LOCAL FUNCTIONS and their descriptions:
        LOGICAL  CHKVAR  ! Checks for the output vars in the input var list

C.........  Local variables
        INTEGER*4 EDDATE        ! Ending date of output file (YYYYDDD)
        INTEGER*4 EDTIME        ! Ending time of output file (HHMMSS)
        INTEGER*4 IBD           ! Start date for each output time period (YYDDD)
        INTEGER*4 IED           ! End date for output time period (YYDDD)
        INTEGER*4 IERR          ! Error flag, 0 = no error
        INTEGER   IFTYPE(10)    ! Output filename in integer format
        INTEGER   IFNOTE(60)    ! Output note in integer format
        INTEGER*4 ISPTR(MXSPEC) ! location of output vars in input vars list
        INTEGER*4 JDATE         ! Date variable (YYYYDDD)
        INTEGER*4 JDATEP1       ! Date variable plus one time step (YYYYDDD)
        INTEGER*4 JTIME         ! Time variable (HHMMSS)
        INTEGER*4 JTIMEP1       ! Time variable plus one time step (HHMMSS)
        INTEGER*4 NSEG          ! Always 1
        INTEGER*4 ODEV          ! UAM file desigator for output file
        INTEGER*4 TSTEP         ! Time step (HHMMSS)
        INTEGER*4 UDEV          ! Userin file desigator for applicable UAM run

        INTEGER*4 I, ILOOP, IVAR, IXY, J, LDEV

        REAL*4    BT            ! Beginning time for a time step period
        REAL*4    ET            ! Ending time for a time step period
        REAL      VAR(NCOLS, NROWS) ! Input and output variable
  
        CHARACTER*16   INAME   !  NetCDF logical input file name
        CHARACTER*10   ONAME   !  UAM file name label for output file
        CHARACTER*1    TEMPCH  !  For integer conversion
        CHARACTER*256  MESG    !  Temporary message array

        LOGICAL   STOPLOOP     !  For stopping time loop

C.........  Data statements
        DATA IXY  /0/
        DATA NSEG /1/

C***********************************************************************
C   begin body of program SMK2EMIS

        LDEV = INIT3()

        WRITE( *,92000 )
     &  ' ',
     &  'Program SMK2EMIS to take the NetCDF 1 layer gridded data file',
     &  'and a USERIN FILE to produce a UAM input file of single, ',
     &  'layer data.  The program ensures that the spatial and',
     &  'temporal properties of the input file and the information.',
     &  'in the USERIN file are compatible.',
     &  ' '
        WRITE( *,92000 )
     &  'You will need to have set the logical names for the input and',
     &  'output files (using "setenv <logicalname> <pathname>") prior',
     &  'to program launch, and will need to enter these logical names',
     &  'at the appropriate prompts.  Input files must have been ',
     &  'sorted as indicated, prior to program execution.',
     &  ' '
        WRITE( *,92000 )
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses to prompts',
     &  'are indicated in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program SMK2EMIS', ' '
            CALL EXIT( 0 )
        ENDIF

C.........  Prompt for name of NetCDF input file

        INAME = PROMPTMFILE(
     &          'Enter logical name for input (NetCDF) file',
     &          FSREAD3, 'PGTSC', 'SMK2EMIS' )

C.........  Prompt for name of USERIN file

        UDEV = PROMPTFFILE(
     &         'Enter logical name for USERIN file applying to UAM run',
     &         .TRUE., .TRUE., 'USERIN', 'SMK2EMIS' )

C.........  Prompt for descriptive name of output file     

        MESG = 'Enter UAM file label for output file'
        CALL GETSTR( MESG, 'EMISSIONS', ONAME )

C.........  Prompt for name of output file

        ODEV = PROMPTFFILE(
     &         'Enter logical name for ' // ONAME // 
     &         ' file applying to UAM run',
     &         .FALSE., .FALSE., 'UAMGTSCA', 'SMK2EMIS' )

C.........  Read NetCDF header information:
C.........  coordinate info, date, time, timestep , and variables

        IF( .NOT. DESC3( INAME ) ) THEN
            CALL M3ERR( 'SMK2EMIS', 0, 0,
     &                  'Could not get description of file ' // INAME,
     &                  .TRUE. )
        END IF

C.........  Read the EPISODE packet from the USERIN file

        CALL RDCNTL( IERR, UDEV )                    ! Note: EPS function
        IF ( IERR .NE. 0 ) 
     &      CALL M3ERR( 'SMK2EMIS', 0, 0,
     &                  'Error reading EPISODE packet' ,
     &                  .TRUE. )

C.........  Read the UAMREGN packet from the USERIN file

        CALL RDREGN( IERR, UDEV )                    ! Note: EPS function
        IF (IERR .NE. 0)
     &      CALL M3ERR( 'SMK2EMIS', 0, 0,
     &                  'Error reading UAMREGN packet' ,
     &                  .TRUE. )

C.........  Compare the NetCDF header information with the USERIN information
C..a) Check variables and set up pointer for correct order

        DO 11 IVAR = 1,NSPECS

            IF( .NOT. CHKVAR(SPCNAM(IVAR),VNAME3D,
     &                       NVARS3D,ISPTR(IVAR))) THEN
         
                IF( .NOT. GETYN( SPCNAM(IVAR) // ' not found '//
     &                           'in input.  Zero and continue?', 
     &                           .TRUE. ) ) THEN
                    WRITE( *,92000 ) ' ', 'Ending program CSGAREA', ' '
                    CALL EXIT( 0 ) 
                ENDIF

            ENDIF
   11   CONTINUE

C..b) Check remaining parameters and initialize JTIME and JDATE

        CALL CHKGRD( IERR )
        IF (IERR .NE. 0) 
     &       CALL M3ERR( 'SMK2EMIS', 0, 0,
     &                   'Fatal difference between input and output ' //
     &                   'grid file specifications', .TRUE. )

        CALL CHKTMP( IERR, JDATE, EDDATE, JTIME, EDTIME, TSTEP )
        IF (IERR .NE. 0) 
     &      CALL M3ERR( 'SMK2EMIS', 0, 0,
     &                  'Fatal difference between input and output ' //
     &                  'time file specifications', .TRUE. )

C.........  Write UAM header information:
C..a) Convert output file name to integer format (as in UAM)

        DO 22 ILOOP= 1, 10
            TEMPCH = ONAME(ILOOP:ILOOP)
            READ(TEMPCH, 91000) IFTYPE(ILOOP)
   22   CONTINUE

C..b) Convert FNOTE (from EPISODE packet) to integer format (as in UAM)

        DO 33 ILOOP= 1,60
            TEMPCH = FNOTE(ILOOP:ILOOP)
            READ(TEMPCH, 91000) IFNOTE(ILOOP)
   33   CONTINUE

C..c) Convert output species names to integer format

        CALL TRANS(IERR)  !  set INTNAM              ! note: EPS subroutine
        IF (IERR .NE. 0)
     &      CALL M3ERR( 'SMK2EMIS', 0, 0,
     &                  'Error converting species from character to ' //
     &                  'integer format', .TRUE. )

C..d) Actually do writing

        WRITE(ODEV,ERR=7010) IFTYPE, IFNOTE, NSEG, NSPECS, IBGDAT,
     &                       BEGTIM, IENDAT, ENDTIM

        WRITE(ODEV,ERR=7010) UTMX, UTMY, IZONE,
     &                       XORIG, YORIG, DXCELL,
     &                       DYCELL, NXCELL, NYCELL, NZCELL,
     &                       NZLOWR, NZUPPR, HTSUR, HTLOWR, HTUPPR 


        WRITE(ODEV,ERR=7010) IXY, IXY, NXCELL, NYCELL

C.........  Write out integer form of output variable names

        WRITE(ODEV,ERR=7010) ((INTNAM(J,I),I=1,10),J=1,NSPECS)

C.........  Initialize time and date for time step after JTIME and JDATE

        JDATEP1 = JDATE
        JTIMEP1 = JTIME
        CALL NEXTIME (JDATEP1, JTIMEP1, TSTEP)
        STOPLOOP = .FALSE.

C.........  Loop over time steps

        DOWHILE ( .NOT. STOPLOOP)

C.............  Initialize output file special time parameters

            BT  = REAL(JTIME/10000)
            ET  = BT + 1.
            IBD = JDATE   - 1900000
            IED = JDATEP1 - 1900000 

            WRITE (*,93000) JDATE, JTIME

C.............  Set UAM formatted times and dates

            WRITE(ODEV,ERR=7011) IBD, BT, IED, ET

C.............  Loop over output variables

            DO 50 IVAR = 1, NSPECS     
        
C.................  Ensure that variable is in input file

                IF ( ISPTR(IVAR) .NE. 0 ) THEN

C.....................  Read input file for time and species of interest

                    IF ( .NOT. READ3( INAME, VNAME3D(ISPTR(IVAR)), 1, 
     &                            JDATE, JTIME, VAR ) ) THEN
                        CALL M3ERR( 'SMK2EMIS', 0, 0,
     &                       'Error reading ' // VNAME3D(ISPTR(IVAR)) // 
     &                       ' from file ' // INAME , .TRUE. )
                    ENDIF 

                ELSE

C.....................  Set variable to zero

                    DO 48 J = 1, NYCELL
                        DO 44 I = 1, NXCELL
                            VAR( I,J ) = 0.e0
   44                   CONTINUE
   48               CONTINUE

                ENDIF

C.............  Write output data for all variables

                WRITE(ODEV,ERR=7011) NSEG, (INTNAM(IVAR,J),J=1,10),
     &                               ((VAR(I,J),I=1,NXCELL), J=1,NYCELL)

   50       CONTINUE       ! Output variables loop

C.............  Check for stopping loop: 
C.............  last output is (output end date/time) - (1 step)

            IF( JDATE .EQ. EDDATE ) THEN
                IF( JTIME .EQ. EDTIME ) STOPLOOP = .TRUE.
            ENDIF

C.............  Update times and dates for next pass through loop

            JDATE = JDATEP1
            JTIME = JTIMEP1
            CALL NEXTIME( JDATEP1, JTIMEP1, TSTEP) 

        ENDDO ! Time step loop ( DOWHILE )

        CALL M3EXIT( 'SMK2EMIS', 0, 0,
     &               'Normal completion of PROGRAM SMK2EMIS', 0 )


C******************    ERROR MESSAGES     ******************************

7010  WRITE (*,*) 'ERROR:  WRITING HEADER OF OUTPUT'
      STOP
7011  WRITE (*,*) 'ERROR:  WRITING RECORD OUTPUT'
      STOP

C******************  FORMAT  STATEMENTS   ******************************

91000 FORMAT (A1)
92000 FORMAT (5X, A )
93000 FORMAT (/, 5X, 'Write out for date ',i7,' and time ',i6,'.')
      END

