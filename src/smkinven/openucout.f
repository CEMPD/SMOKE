
        SUBROUTINE OPENUCOUT( UONAME, UPNAME, UENAME )

C*************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine opens the uncertainty output files.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines
C      Functions: I/O API functions
C
C  REVISION  HISTORY:
C      Created 9/2001 by A. Holland
C
C*************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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

C.........  MODULES for public variables
C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO
        
C.........  This module contains uncertainty-specific settings        
        USE MODUNCERT

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=NAMLEN3), INTENT(OUT) :: UONAME ! uncertainty i/o api
        CHARACTER(LEN=NAMLEN3), INTENT(OUT) :: UPNAME ! parametric output 
        CHARACTER(LEN=NAMLEN3), INTENT(OUT) :: UENAME ! empirical output

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag
 

C...........   Other local variables

        INTEGER       I, J, L, LV           ! counter and indices

        INTEGER       IOS                   ! i/o status

        LOGICAL, SAVE :: FIRSTIME = .TRUE.    ! true:  first time
        LOGICAL, SAVE :: SECONDTIME = .TRUE. ! true:  second time

        CHARACTER*300 MESG               ! message buffer 
        CHARACTER*10  VBUFFER            ! variable buffer

        CHARACTER(LEN=NAMLEN3)  NAMBUF   ! file name buffer
        CHARACTER(LEN=IOVLEN3)  VNAME
        CHARACTER(LEN=IOULEN3)  UNITS    ! tmp units name
        CHARACTER(LEN=IOVLEN3)  CBUF     ! tmp pollutant name
        
        CHARACTER*16 :: PROGNAME = 'OPENUCOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENUCOUT

        IF( FIRSTIME ) THEN
        
            MESG = 'Opening uncertainty output file...'
            CALL M3MSG2( MESG )

C.........  Set up for opening I/O API output file header

            CALL HDRMISS3  ! Initialize for emissions 

            NVARS3D = ( 4 * NUOVAR ) + 1
            NROWS3D = UCOUNT  !  number of rows = # of sources.

            FDESC3D( 1 ) = CATDESC // ' Uncertainty'
            FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
            FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
            WRITE( FDESC3D( 4 ),94010 ) '/NON POLLUTANT/ ', 1

            IF( NIPOL .GT. 0 ) THEN
               WRITE( FDESC3D( 5 ),94010 ) '/VARIABLES/', NUOVAR
               WRITE( FDESC3D( 6 ),94010 ) '/PER VARIABLE/ ', 4
            END IF

            IF( NIACT .GT. 0 ) THEN
                WRITE( FDESC3D( 7 ),94010 ) '/ACTIVITIES/', NIACT
                WRITE( FDESC3D( 8 ),94010 ) '/PER ACTIVITY/ ', 4
            END IF


C.........  Define source characteristic variables that are not strings

            J = 1
            VNAME3D( J ) = 'SRCNUM'
            VTYPE3D( J ) = M3INT
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Source number'
            J = J + 1

            DO I = 1, NUOVAR
            
                CBUF = UONAMES( I )
                L = LEN_TRIM( CBUF )

                VNAME3D( J ) = 'MTH_'//CBUF( 1:L )       
                VTYPE3D( J ) = M3INT
                UNITS3D( J ) = 'n/a'
                VDESC3D( J ) = 'Empirical (0) or Parametric (1)'
                J = J + 1

                VNAME3D( J ) = 'TYP_'//CBUF( 1:L )       
                VTYPE3D( J ) = M3INT
                UNITS3D( J ) = 'n/a'
                VDESC3D( J ) = 'Normal=1,Logarithmic=2,Gamma=3'//
     &          ',Weibul=4,Beta=5,Stepwise=6,Linear Interpretation=7'
                J = J + 1 

                VNAME3D( J ) = 'NEP_'//CBUF( 1:L )       
                VTYPE3D( J ) = M3INT
                UNITS3D( J ) = 'n/a'
                VDESC3D( J ) = 'Number of parameters or empirical'//
     &                         ' entries'                
                J = J + 1
        
                VNAME3D( J ) = 'UIX_'//CBUF( 1:L )
                VTYPE3D( J ) = M3INT
                UNITS3D( J ) = 'n/a'
                VDESC3D( J ) = 'Row number of empirical or parametric'//
     &                     ' entry in corresponding output file'
                J = J + 1
                
            END DO
        
        
C.........  Prompt for and open I/O API output file

            MESG = 'Enter logical name for the I/O API ' //
     &             'UNCERTAINTY output file'
     

            NAMBUF = PROMPTMFILE( MESG, FSUNKN3, CRL//'UCOUT',
     &                            PROGNAME )

            UONAME = NAMBUF
            
            FIRSTIME = .FALSE.
            
        ELSE IF( SECONDTIME ) THEN
        
            MESG = 'Opening uncertainty output parametric file...'
            CALL M3MSG2( MESG )

C.........  Set up for opening I/O API output file header    
        
            CALL HDRMISS3  ! Initialize for emissions 

            NVARS3D = 1
            NROWS3D = NPPCKT    !  number of rows = # of parametric packets.
            NCOLS3D = MXPARDAT  !  number of cols = max # of parametric entries

            FDESC3D( 1 ) = CATDESC // ' Uncertainty Parametric Data'
            FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
            FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

C.........  Define source characteristic variables that are not strings

            J = 1
            VNAME3D( J ) = 'PARMS'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Parameters'
            J = J + 1

        
C.........  Prompt for and open I/O API output file

            MESG = 'Enter logical name for the I/O API ' //
     &             'UNCERTAINTY parameters output file'
     

            NAMBUF = PROMPTMFILE( MESG, FSUNKN3, CRL//'UCPOUT',
     &                            PROGNAME )

            UPNAME = NAMBUF
            
            SECONDTIME = .FALSE.
            
        ELSE
        
            MESG = 'Opening uncertainty output empirical file...'
            CALL M3MSG2( MESG )
            
C.........  Set up for opening I/O API output file header    
        
            CALL HDRMISS3  ! Initialize for emissions 

            NVARS3D = 2
            NROWS3D = NEPCKT    !  number of rows = # of empirical packets.
            NCOLS3D = MXEMPDAT  !  number of cols = max # of empirical entries

            FDESC3D( 1 ) = CATDESC // ' Uncertainty Empirical Data'
            FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
            FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

C.........  Define source characteristic variables that are not strings

            J = 1
            VNAME3D( J ) = 'EMFVAL'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Emission factor value'
            J = J + 1

            VNAME3D( J ) = 'PROBVAL'       
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Probability'
            J = J + 1

        
C.........  Prompt for and open I/O API output file

            MESG = 'Enter logical name for the I/O API ' //
     &             'UNCERTAINTY empirical output file'
     

            NAMBUF = PROMPTMFILE( MESG, FSUNKN3, CRL//'UCEOUT',
     &                            PROGNAME )

            UENAME = NAMBUF
            
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENUCOUT

