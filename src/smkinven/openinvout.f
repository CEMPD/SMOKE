
        SUBROUTINE OPENPNTS( NPSRC, NIPOL, EINAM, ENAME, SDEV )

C***********************************************************************
C  subroutine body starts at line 99
C
C  DESCRIPTION:
C      This subroutine sets up the header and variables for the I/O API 
C      inventory file, and opens the I/O API and ASCII files for the SMOKE
C      point source inventory.
C
C  PRECONDITIONS REQUIRED:
C      Correct number of sources NPSRC is set
C      Correct number of pollutants and names EINAM are set
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, BLDENAMS
C      Functions: I/O API functions, VERCHAR
C
C  REVISION  HISTORY:
C      Created 10/98 by M. Houyoux
C
C****************************************************************************/
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        CHARACTER*2     CRLF
        INTEGER         PROMPTFFILE
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE
        CHARACTER*16    VERCHAR

        EXTERNAL CRLF, PROMPTFFILE, PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        INTEGER       NPSRC          ! Actual source count
        INTEGER       NIPOL          ! Actual no of inv pollutants
        CHARACTER(LEN=IOVLEN3) EINAM( NIPOL ) ! Name of actual pollutants
        CHARACTER(LEN=NAMLEN3) ENAME ! Emissions output inventory logical name
        INTEGER       SDEV           ! For ASCII output inventory file

C...........   LOCAL PARAMETERS
        CHARACTER*50  SCCSW          ! SCCS string with version number at end

        PARAMETER   ( SCCSW   = '@(#)$Id$'
     &              )

C...........   Names, Units, types, & descriptions for pollutant-specific 
C              output variables

        CHARACTER(LEN=IOVLEN3) EONAMES( NIPOL,NPTPPOL3 ) ! Names 
        INTEGER                EOTYPES( NIPOL,NPTPPOL3 ) ! Types (Real or Int) 
        CHARACTER(LEN=IOULEN3) EOUNITS( NIPOL,NPTPPOL3 ) ! Units  
        CHARACTER(LEN=IODLEN3) EODESCS( NIPOL,NPTPPOL3 ) ! Dscriptions  

C...........   Other local variables

        INTEGER       I, J, L1, L2, V     ! counter and indices
        INTEGER       NIOVARS   ! Number of I/O API file non-emis variables
        INTEGER       NPOLMAX   ! Max no of pollutants, based on I/O API

        CHARACTER*300 MESG      ! message buffer 

        CHARACTER*16 :: PROGNAME = 'OPENPNTS' ! program name

C***********************************************************************
C   begin body of subroutine OPENPNTS

C.........  Check number of output variables against I/O API maximum

        NIOVARS = NPTVAR3 + NPTPPOL3 * NIPOL
        NPOLMAX = INT( ( MXVARS3 - NPTVAR3 ) / NPTPPOL3 )

C.........  If there are too many output variables, reset NIPOL
        IF( NIOVARS .GT. MXVARS3 ) THEN

            WRITE( MESG,94010 ) 
     &             'WARNING: Maximum number of pollutants that can ' //
     &             'be written to the' // CRLF() // BLANK5 //
     &             '         I/O API file is', NPOLMAX, 
     &             '. This limitation is caused by' // CRLF()// BLANK5//
     &             '         the I/O API variable limit of', MXVARS3,'.'
            CALL M3MSG2( MESG )
 
            WRITE( MESG,94010 ) 
     &             'WARNING: Reseting number of output pollutants to ',
     &             NPOLMAX
            CALL M3MSG2( MESG )

            NIPOL   = NPOLMAX
            NIOVARS = NPTVAR3 + NPTPPOL3 * NIPOL

        ENDIF

C.........  Set up for opening I/O API output file header

        FTYPE3D = GRDDED3
        P_ALP3D = DBLE( AMISS3 )
        P_BET3D = DBLE( AMISS3 )
        P_GAM3D = DBLE( AMISS3 )
        XCENT3D = 0.0D0
        YCENT3D = 0.0D0
        XORIG3D = DBLE( AMISS3 )
        YORIG3D = DBLE( AMISS3 )
        SDATE3D = 0       !  n/a
        STIME3D = 0       !  n/a
        TSTEP3D = 0       !  time independent
        NVARS3D = NIOVARS
        NCOLS3D = 1
        NROWS3D = NPSRC   !  number of rows = # of point sources.
        NLAYS3D = 1
        NTHIK3D = 1
        GDTYP3D = IMISS3
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3
        GDNAM3D = ' '

        FDESC3D( 1 ) = 'Point source inventory'
        FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )
        WRITE( FDESC3D( 4 ),94010 ) '/NON POLLUTANT/ ', NPTVAR3
        WRITE( FDESC3D( 5 ),94010 ) '/PER POLLUTANT/ ', NPTPPOL3 

        DO I = 6, MXDESC3
            FDESC3D( I ) = ' '
        ENDDO

C.........  Define source characteristic variables that are not strings

        J = 1
        VNAME3D( J ) = 'IFIP'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'State and county FIPS code'
        J = J + 1

        VNAME3D( J ) = 'ISIC'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Source Industrial Code'
        J = J + 1

        VNAME3D( J ) = 'ISCC'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Source Classification Code'
        J = J + 1

        VNAME3D( J )= 'IORIS'
        VTYPE3D( J )= M3INT
        UNITS3D( J )= 'n/a'
        VDESC3D( J )= 'Office of the Regulatory Information System code'
        J = J + 1

        VNAME3D( J ) = 'TZONES'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Time zone for site'
        J = J + 1

        VNAME3D( J ) = 'TPFLAG'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'T|2? T|3?'
        VDESC3D( J ) = 'Use week(2), month(3) temporal profiles or not'
        J = J + 1

        VNAME3D( J ) = 'INVYR'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'year AD'
        VDESC3D( J ) = 'Year of inventory for this record'
        J = J + 1

        VNAME3D( J ) = 'XLOCA'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'degrees'
        VDESC3D( J ) = 'longitude'
        J = J + 1

        VNAME3D( J ) = 'YLOCA'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'degrees'
        VDESC3D( J ) = 'latitude'
        J = J + 1

        VNAME3D( J ) = 'STKHT'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'Stack height'
        J = J + 1

        VNAME3D( J ) = 'STKDM'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'Stack diameter'
        J = J + 1

        VNAME3D( J ) = 'STKTK'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'deg K'
        VDESC3D( J ) = 'Stack exhaust temperature'
        J = J + 1

        VNAME3D( J ) = 'STKVE'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm/s'
        VDESC3D( J ) = 'Stack exhaust velocity'
        J = J + 1

C.........  Get names, units, etc. of output pollutant-specific records
        CALL BLDENAMS( 'POINT', NIPOL, NPTPPOL3, EINAM, 
     &                 EONAMES, EOUNITS, EOTYPES, EODESCS )

        DO V = 1 , NIPOL
            
            DO I = 1, NPTPPOL3 ! Loop through number of variables per pollutant

                VNAME3D( J ) = EONAMES( V, I )
                VTYPE3D( J ) = EOTYPES( V, I )
                UNITS3D( J ) = EOUNITS( V, I )
                VDESC3D( J ) = EODESCS( V, I )
                J = J + 1

            ENDDO    !  end loop on number of variables per pollutant

        ENDDO        !  end loop on inventory pollutants V

C.........  Prompt for and open I/O API output file
        ENAME= PROMPTMFILE( 
     &       'Enter logical name for the I/O API INVENTORY output file',
     &       FSUNKN3, 'PNTS', PROGNAME )
        
C.........  Prompt for and open ASCII output file
        SDEV= PROMPTFFILE( 
     &      'Enter logical name for the ASCII INVENTORY output file',
     &      .FALSE., .TRUE., 'PSRC', PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END

