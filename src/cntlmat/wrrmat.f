
        SUBROUTINE WRRMAT( NSRC, NMSPC, FDEV, FILNAM, INDX, REPEM, 
     &                     PRJFAC, MKTPEN, RMTX, CSCC, 
     &                     SPROF, VNAMES )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine writes the reactivity matrices
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 3/99 by M. Houyoux
C
C****************************************************************************/
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
C***************************************************************************

        USE M3UTILIO

      IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.........  Subroutine arguments and their descriptions:
        INTEGER     , INTENT (IN) :: NSRC           ! no. of source
        INTEGER     , INTENT (IN) :: NMSPC          ! no. model species
        INTEGER     , INTENT (IN) :: FDEV           ! supplement file unit no.
        CHARACTER(*), INTENT (IN) :: FILNAM         ! I/O API file name
        INTEGER     , INTENT (IN) :: INDX  ( NSRC ) ! pollutant count
        REAL        , INTENT (IN) :: REPEM ( NSRC ) ! replacement emissions
        REAL        , INTENT (IN) :: PRJFAC( NSRC ) ! projection factors
        REAL        , INTENT (IN) :: MKTPEN( NSRC ) ! market penetration
        REAL        , INTENT (IN) :: RMTX  ( NSRC,NMSPC ) ! speciation data
        CHARACTER(*), INTENT (IN) :: CSCC  ( NSRC ) ! reactivity SCCs
        CHARACTER(*), INTENT (IN) :: SPROF ( NSRC ) ! speciation profile IDs
        CHARACTER(*), INTENT (IN) :: VNAMES( NMSPC ) ! variable names for spcs

C...........   Other local variables
        INTEGER                 I, L, L2, S

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time subroutine called

        CHARACTER(IOVLEN3) VAR     !  tmp variable names

        CHARACTER(100) OUTFMT           !  format buffer
        CHARACTER(300) MESG             !  message buffer

        CHARACTER(16) :: PROGNAME = 'WRRMAT' !  program name

C***********************************************************************
C   begin body of program WRRMAT

C.........  The first time this routine is called, write out the supplement file
        IF( FIRSTIME ) THEN

C.............  Generate format
            WRITE( OUTFMT, 94020 ) SCCLEN3, SPNLEN3
            
C.............  Write out header
            WRITE( FDEV, 93010 ) 3, OUTFMT( 1:LEN_TRIM( OUTFMT ) )
            WRITE( FDEV, 93000 ) 'SRCID   Source index'
            WRITE( FDEV, 93000 ) 'CSCC    Source category code '
            WRITE( FDEV, 93000 ) 'SPROF   Speciation profile number'

C.............  Write out records for the ASCII supplement file
            DO S = 1, NSRC
                WRITE( FDEV, OUTFMT ) INDX( S ), CSCC( S ), SPROF( S )
            END DO

            FIRSTIME = .FALSE.

        END IF

C.........  Initialize message to use in case there is an error

        MESG = 'Problem writing to output file "' //
     &         FILNAM( 1:LEN_TRIM( FILNAM ) ) // '"'

        L = LEN_TRIM( MESG )

C.........  Write the I/O API variables for the non-speciation data

        IF( .NOT. WRITESET( FILNAM,'SRCID',ALLFILES, 0,0,INDX ) ) THEN
            MESG = MESG( 1:L ) // ' for variable "SRCID"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( .NOT. WRITESET( FILNAM,'REPEMIS',ALLFILES,0,0,REPEM ) ) THEN
            MESG = MESG( 1:L ) // ' for variable "REPEMIS"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( .NOT. WRITESET( FILNAM,'PRJFAC',ALLFILES,0,0,PRJFAC ) ) THEN
            MESG = MESG( 1:L ) // ' for variable "PRJFAC"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( .NOT. WRITESET( FILNAM,'MKTPEN',ALLFILES,0,0,MKTPEN ) ) THEN
            MESG = MESG( 1:L ) // ' for variable "MKTPEN"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Write the I/O API variables for the speciation factors

        DO I = 1, NMSPC

            VAR = VNAMES( I )
            L2  = LEN_TRIM( VAR )

            IF( .NOT. WRITESET( FILNAM, VAR, ALLFILES, 
     &                          0, 0, RMTX( 1,I )     ) ) THEN
                MESG = MESG( 1:L  ) // ' for variable "' //
     &                 VAR ( 1:L2 ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( I2, A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( '(I6,1X,A', I2.2, ',1X,A', I2.2, ')' )

        END SUBROUTINE WRRMAT
