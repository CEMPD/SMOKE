OPTIONS LINESIZE=110;
/* This program reads in fuel use data for a particular year
   and a cross reference file to find matching OTAG id info for
   a device.   Then the program takes these matches and finds out
   if there is actual fuel use data for all of the summer season
   months.   It keeps only the matches that have fuel use info 
   for ALL months specified.

Prototype :  Jeff Vukovich 10/97


*/
/*
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1997, MCNC--North Carolina Supercomputing Center
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
*/

/* read in filenames */

DATA _NULL_;
CALL SYMPUT ('filnm',((SCAN("&sysparm",1,";"))));
CALL SYMPUT ('filnm2',((SCAN("&sysparm",2,";"))));
stop;
RUN;

/* read in cross reference file */

DATA xref (keep=stid cyid pltid dvid doeid blrid );

   INFILE "&filnm" DLM=',';
   INPUT  stid   2.       /* STID  code */
          cyid   3.     /* CYID  code */
          pltid  $15.     /* Round3 zone*/
          dvid   $12.
          doeid  $5.
          blrid  $5.;
RUN;

PROC SORT DATA=xref ;
by stid cyid blrid doeid pltid dvid; 
RUN;


/* read in fuel use data for a specific year */

DATA fu95 (keep= stid cyid name doeid blrid scc may jun jul aug sep);

   INFILE "&filnm2" DLM=',';
   INPUT  stid   2.       /* STID  code */
          cyid   3.     /* CYID  code */
          name   $60.     /* Round3 zone*/
          doeid  $5.
          blrid  $5.
          scc    8.
          jan   11.
          feb   11.
          mar   11.
          apr   11.
          may   11.
          jun   11.
          jul   11.
          aug   11.
          sep   11.
          oct   11.
          nov   11.
          dec   11.;
RUN;

PROC SORT DATA=fu95;
by stid cyid blrid doeid ;
RUN;


DATA fu95a(keep=stid cyid blrid doeid);
set fu95;
RUN;

/* find matches */

DATA fu95x ;
 MERGE fu95a (IN=a) xref (IN=b);
 by stid cyid blrid doeid;
 if b and not a then delete;
/* if blrid eq . then delete;
 if blrid eq ' ' then delete;
 if doeid eq . then delete;
 if doeid eq ' ' then delete;
 if may eq . then k=k+1;
 if jun eq . then k=k+1;
 if jul eq . then k=k+1;
 if aug eq . then k=k+1;
 if sep eq . then k=k+1;
*/
RUN;

PROC SORT DATA=fu95x;
by stid cyid blrid doeid ;
RUN;

DATA fu95x2;
set fu95x;
if pltid eq ' ' then delete;
if blrid eq ' ' then delete;
RUN;

/* find matches that have fuel use info for all months specified */

DATA fu95x2;
MERGE fu95x2(IN=a) fu95(IN=b);
by stid cyid blrid doeid ;
if b and not a then delete;
k=0;
kb=0;
if may eq 0 or may eq . then k=k+1;
if may eq . then kb=kb+1;
if jun eq 0 or jun eq . then k=k+1;
if jun eq . then kb=kb+1;
if jul eq 0 or jul eq . then k=k+1;
if jul eq . then kb=kb+1;
if aug eq 0 or aug eq . then k=k+1;
if aug eq . then kb=kb+1;
if sep eq 0 or sep eq . then k=k+1;
if sep eq . then kb=kb+1;
if kb ge 1 then delete;
if k ge 1 then delete; 

RUN;

PROC SORT DATA=fu95x2 NODUPREC FORCE;
BY stid cyid pltid scc dvid ;
RUN;

/* print out info */
 
DATA _null_;
 SET fu95x2 END=EOF;
 FILE print   NOTITLES;

  PUT   @1 stid 2.
        @3 cyid 3.  
        @6 pltid $15.
        @21 dvid $12.
        @33 scc   8. 
        @41 may  10.2
        @51 jun  10.2
        @61 jul  10.2
        @71 aug  10.2
        @81 sep  10.2;
RETURN;

 
