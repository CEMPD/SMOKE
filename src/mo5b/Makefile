# Version "@(#)$Id$ $Source$ $Date$ 

#
# MKMF template makefile for protected mode executables.
#
FC        = f77
LINKER    = f77
PROGRAM   = mobile5b.exe
DEST      = .
EXTHDRS   =
FFLAGS    =
HDRS      =
LDFLAGS   = -STUB RUNB -MAXDATA 0
LDMAP     = nul
LIBS      =
MAKEFILE  = Makefile
OBJS            = a8to11.o adj12.o adjust.o atpaer.o atpeff.o \
                  bageme.o bd01.o bd02.o bd03.o bd04.o bd05.o \
                  bd10.o bd11.o bd12.o bd13.o bd14.o bd15.o \
                  bd16.o bd17.o bd18.o bd19.o bd20.o bd21.o \
                  bd22.o bd23.o bd24.o bd25.o bd26.o bd27.o \
                  bd28.o bd29.o bd30.o bd31.o bd32.o bd33.o \
                  bd34.o bd35.o bd36.o bd37.o bd38.o bd39.o \
                  bdsave.o bef.o befun.o bigalh.o bigcfx.o \
                  bigin1.o bigin2.o bigmet.o bigomt.o bigrvp.o \
                  bigsal.o bigsc1.o bigsc2.o bigsc3.o bigtam.o \
                  bigtcf.o c1to7.o c8to12.o calaer.o calcum.o \
                  caltem.o caludi.o ccevrt.o cevmpd.o cevtpd.o \
                  chkatp.o chkfun.o chkpgm.o consec.o crankc.o \
                  cvtpdf.o disatp.o discal.o disimp.o diurnl.o \
                  driver.o efcalx.o effgrp.o efnftp.o emirat.o \
                  enforc.o etpcal.o etpred.o ev81.o ev85.o \
                  evhi.o evmain.o evpeme.o evpre.o evstd.o \
                  evunc.o exmptr.o fail.o findpb.o fuel.o \
                  fuel2.o getalh.o getatp.o getbym.o getcol.o \
                  getcum.o getefr.o geterq.o getetp.o getfyr.o \
                  gethor.o getim2.o getimc.o getimp.o getiou.o \
                  getlev.o getmrn.o getnat.o getppd.o getprg.o \
                  getprs.o getsc1.o gettam.o getttc.o getver.o \
                  getvmt.o gsfidx.o hccalx.o hdcgpm.o hdcpic.o \
                  hddmym.o hdrcsv.o hotsok.o ibfptr.o iccopt.o \
                  idlcal.o ierptr.o ievptr.o igeptr.o iglptr.o \
                  igrvpt.o igsfpt.o igtptr.o iguptr.o imovlp.o \
                  imptr.o initap.o irulpt.o ispptr.o itampt.o \
                  itcptr.o levef.o local.o mobile.o myratp.o \
                  myrimt.o nparam.o obyhdr.o obyimp.o obytcn.o \
                  onesec.o otcalc.o outatp.o outby1.o outby2.o \
                  outby3.o outby4.o outcsv.o outdt1.o outdt2.o \
                  outdt3.o outdt4.o outdt5.o outdtg.o outetp.o \
                  outgsf.o outhd1.o outhd2.o outhd3.o outhd4.o \
                  outhdr.o outimp.o outlap.o outmth.o outnew.o \
                  outpol.o output.o outstg.o outtam.o outzev.o \
                  outztt.o ovlp2.o p3wpox.o parsec.o pcleft.o \
                  pctlev.o piv.o quiter.o rdcity.o refuel.o \
                  regmod.o reinit.o restor.o rstlos.o ruloss.o \
                  rvpexh.o savepb.o saver.o scfsub.o setegs.o \
                  stoerq.o tamper.o tech0.o tfcalx.o togfid.o \
                  usdptr.o vocfid.o wtlrul.o yrtest.o zcalc.o
SRCS            = a8to11.F adj12.F adjust.F atpaer.F atpeff.F \
                  bageme.F bd01.F bd02.F bd03.F bd04.F bd05.F \
                  bd10.F bd11.F bd12.F bd13.F bd14.F bd15.F \
                  bd16.F bd17.F bd18.F bd19.F bd20.F bd21.F \
                  bd22.F bd23.F bd24.F bd25.F bd26.F bd27.F \
                  bd28.F bd29.F bd30.F bd31.F bd32.F bd33.F \
                  bd34.F bd35.F bd36.F bd37.F bd38.F bd39.F \
                  bdsave.F bef.F befun.F bigalh.F bigcfx.F \
                  bigin1.F bigin2.F bigmet.F bigomt.F bigrvp.F \
                  bigsal.F bigsc1.F bigsc2.F bigsc3.F bigtam.F \
                  bigtcf.F c1to7.F c8to12.F calaer.F calcum.F \
                  caltem.F caludi.F ccevrt.F cevmpd.F cevtpd.F \
                  chkatp.F chkfun.F chkpgm.F consec.F crankc.F \
                  cvtpdf.F disatp.F discal.F disimp.F diurnl.F \
                  driver.F efcalx.F effgrp.F efnftp.F emirat.F \
                  enforc.F etpcal.F etpred.F ev81.F ev85.F \
                  evhi.F evmain.F evpeme.F evpre.F evstd.F \
                  evunc.F exmptr.F fail.F findpb.F fuel.F \
                  fuel2.F getalh.F getatp.F getbym.F getcol.F \
                  getcum.F getefr.F geterq.F getetp.F getfyr.F \
                  gethor.F getim2.F getimc.F getimp.F getiou.F \
                  getlev.F getmrn.F getnat.F getppd.F getprg.F \
                  getprs.F getsc1.F gettam.F getttc.F getver.F \
                  getvmt.F gsfidx.F hccalx.F hdcgpm.F hdcpic.F \
                  hddmym.F hdrcsv.F hotsok.F ibfptr.F iccopt.F \
                  idlcal.F ierptr.F ievptr.F igeptr.F iglptr.F \
                  igrvpt.F igsfpt.F igtptr.F iguptr.F imovlp.F \
                  imptr.F initap.F irulpt.F ispptr.F itampt.F \
                  itcptr.F levef.F local.F mobile.F myratp.F \
                  myrimt.F nparam.F obyhdr.F obyimp.F obytcn.F \
                  onesec.F otcalc.F outatp.F outby1.F outby2.F \
                  outby3.F outby4.F outcsv.F outdt1.F outdt2.F \
                  outdt3.F outdt4.F outdt5.F outdtg.F outetp.F \
                  outgsf.F outhd1.F outhd2.F outhd3.F outhd4.F \
                  outhdr.F outimp.F outlap.F outmth.F outnew.F \
                  outpol.F output.F outstg.F outtam.F outzev.F \
                  outztt.F ovlp2.F p3wpox.F parsec.F pcleft.F \
                  pctlev.F piv.F quiter.F rdcity.F refuel.F \
                  regmod.F reinit.F restor.F rstlos.F ruloss.F \
                  rvpexh.F savepb.F saver.F scfsub.F setegs.F \
                  stoerq.F tamper.F tech0.F tfcalx.F togfid.F \
                  usdptr.F vocfid.F wtlrul.F yrtest.F zcalc.F

###   For make to take advantage of compiler response files, uncomment
###   the next four lines, and the additional three lines listed after
###   the '$(PROGRAM):   $(OBJS) $(LIBS)' line.
###.BEFORE:
###       DEL F77L3.RSP
###.F.o:
###       ECHO $<,$@ $(FFLAGS) >>F77L3.RSP


$(PROGRAM):     $(OBJS) $(LIBS)
###       !IF -e F77L3.RSP
###           $(FC) @F77L3.RSP
###       !ENDIF
          $(LINKER) $(OBJS) -EXE $@ -MAP $(LDMAP) -LIB $(LIBS) $(LDFLAGS)

clean:;         @del -f $(OBJS)

depend:;   @mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

install:   $(PROGRAM)
                @echo Installing $(PROGRAM) in $(DEST)
                @if not $(DEST)x==.x copy $(PROGRAM) $(DEST)
### OPUS MKMF:  Do not remove this line!  Automatic dependencies follow.

a8to11.o: LOOKUP.I MYCODE.I SIZCAL.I SIZPAR.I

adj12.o: LOOKUP.I MYCODE.I SIZCAL.I

adjust.o: EVAL.I EVAPGR.I MAXIMA.I RESUL1.I RESUL2.I RESUL3.I VMXCOM.I \
         X1.I X2.I X3.I X4.I X5.I

atpaer.o: ATPAR1.I ATPAR2.I ATPAR3.I ATPAR4.I FLAGS3.I IMPAR1.I IMPAR2.I \
         IMPAR6.I LOOKUP.I MAXIMA.I MYCODE.I

atpeff.o: EGSCAL.I LOOKUP.I MYCODE.I SIZCAL.I

bageme.o: EGSCAL.I EMECAL.I FLAGS4.I LOOKUP.I MAXIMA.I MYCODE.I TAMOUT.I

bd01.o: SPEED1.I

bd02.o: PARAME.I SPEED2.I SPEED3.I SPEED4.I SPEED5.I SPEED6.I SPEED7.I \
         SPEED8.I

bd03.o: BASEQ2.I BASEQ4.I BASEQ6.I PARAME.I

bd04.o: BASEQ1.I BASEQ5.I

bd05.o: BASE10.I BASEQ9.I TIER1S.I

bd10.o: NMETH1.I NMETH2.I PARAME.I

bd11.o: ETPPAR.I EVAHIA.I EVAHS1.I EVAPAR.I EVAPGR.I EVAPHD.I EVASTD.I \
         PARAME.I

bd12.o: SAVE01.I SAVE02.I SAVE03.I SAVE04.I SAVE05.I SAVE06.I

bd13.o: EVADU1.I EVADU2.I

bd14.o: MPDCOM.I PARAME.I TPDCOM.I

bd15.o: CRACOM.I PARAME.I

bd16.o: IOUCOM.I IVPCOM.I MAXIMA.I PARAME.I PARAMS.I QUITXQ.I YEARS4.I

bd17.o: PARAME.I TEMPC1.I TEMPC2.I TEMPC3.I TEMPC4.I TEMPC5.I TEMPC6.I \
         TEMPC7.I TEMPS.I

bd18.o: OPMOD1.I

bd19.o: OPMOD2.I OPMOD3.I PARAME.I

bd20.o: MYRSAV.I VMXCOM.I

bd21.o: ALUHAC.I ALUHLT.I

bd22.o: ALUPAC.I ALUPLT.I PARAME.I

bd23.o: GSFCOM.I PARAME.I

bd24.o: NAMES1.I STRING.I USDATA.I

bd25.o: ATPAR1.I ATPAR2.I ATPAR3.I ATPAR4.I EGSCAL.I EMECAL.I EVPDAT.I \
         LOOKUP.I MYCODE.I MYUB1.I MYUB2.I OMTCOM.I PARAME.I SIZCAL.I \
         SIZPAR.I TAMEQ1.I TAMEQ2.I TAMEQ3.I TAMEQ4.I TAMOU2.I TAMOUT.I \
         TAMPB1.I

bd26.o: HDDMAR.I HDDREG.I

bd27.o: ALUHIN.I ALUOUT.I BASEQ3.I CEVBMY.I CEVCOM.I CUMCOM.I FLAGS1.I \
         FLAGS2.I FLAGS3.I FLAGS4.I IM240P.I IMPAR1.I IMPAR2.I MYRCAL.I \
         OFFSET.I PARAME.I PROJEC.I REGISF.I RESUL1.I RESUL2.I RESUL3.I

bd28.o: PARAME.I RLCOM1.I RLCOM2.I RLCOM3.I

bd29.o: HDCCOM.I

bd30.o: CITCIN.I CITPAR.I CITRV1.I CITRV2.I INJECT.I RFORM1.I RFORM2.I

bd31.o: PARAME.I RVPEX1.I RVPEX2.I RVPNAT.I

bd32.o: IMPAR4.I IMPAR5.I IMPAR6.I IMPAR7.I PPEFFD.I

bd33.o: LEVBLK.I REGION.I SCENE1.I SCENE2.I

bd34.o: PARAME.I RESTLS.I RUNLS1.I RUNLS2.I RUNLS3.I

bd35.o: BYMYC1.I BYMYC2.I BYMYC3.I BYMYC4.I BYMYC5.I

bd36.o: DTGGAS.I OXY1.I OXY2.I

bd37.o: FID.I PARAME.I

bd38.o: EVAL.I X1.I X2.I X3.I X4.I X5.I

bd39.o: TTCINF.I

bdsave.o: ETPPAR.I MAXIMA.I MYRSAV.I SAVE01.I SAVE02.I SAVE03.I SAVE04.I \
         SAVE05.I SAVE06.I TAMEQ1.I VMXCOM.I

bef.o: BYMYC1.I BYMYC2.I BYMYC5.I FLAGS1.I FLAGS2.I FLAGS4.I LEVBLK.I \
         MAXIMA.I OFFSET.I OMTCOM.I OXY1.I REGION.I RFORM1.I RFORM2.I \
         YEARS4.I

befun.o: BASE10.I BASE11.I BASEQ1.I BASEQ5.I BASEQ6.I FLAGS2.I IOUCOM.I \
         LEVBLK.I REGION.I TAMEQ2.I TIER1S.I YEARS4.I

bigalh.o: ALUHAC.I ALUHIN.I ALUHLT.I ALUOUT.I MAXIMA.I REGION.I SCENE1.I \
         TEMPS.I

bigcfx.o: FLAGS2.I FLAGS4.I MAXIMA.I MYRCAL.I TEMPC4.I

bigin1.o: FLAGS4.I MAXIMA.I OFFSET.I SCENE1.I SCENE2.I SPEED3.I TEMPC3.I \
         TEMPC5.I TEMPC6.I TEMPS.I

bigin2.o: CUMCOM.I IVPCOM.I MAXIMA.I OPMOD2.I OPMOD3.I PARAME.I SCENE2.I \
         SPEED4.I SPEED6.I TAMEQ2.I TEMPC3.I

bigmet.o: CUMCOM.I INJECT.I NMETH1.I NMETH2.I OFFSET.I REGION.I TEMPC3.I

bigomt.o: OMTCOM.I OPMOD1.I OPMOD2.I OPMOD3.I TAMEQ2.I TEMPC3.I TEMPC4.I \
         TEMPS.I

bigrvp.o: CITRV1.I CITRV2.I FLAGS1.I MYCODE.I RFORM1.I RFORM2.I RVPEX1.I \
         TEMPC3.I TEMPS.I

bigsal.o: ALUOUT.I FLAGS2.I SCENE1.I SPEED4.I SPEED6.I

bigsc1.o: SCENE1.I SPEED3.I SPEED4.I SPEED6.I SPEED7.I

bigsc2.o: SCENE1.I SPEED1.I SPEED3.I SPEED4.I SPEED6.I SPEED8.I

bigsc3.o: SCENE1.I SPEED1.I SPEED3.I SPEED4.I SPEED6.I SPEED8.I

bigtam.o: OMTCOM.I TAMOUT.I TEMPC3.I

bigtcf.o: CITRV2.I FLAGS4.I INJECT.I MYCODE.I OFFSET.I TEMPC1.I TEMPC2.I \
         TEMPC3.I TEMPC4.I TEMPC5.I TEMPC6.I TEMPC7.I TEMPS.I

c1to7.o: LOOKUP.I MYCODE.I SIZCAL.I

c8to12.o: LOOKUP.I MYCODE.I SIZCAL.I SIZPAR.I

calaer.o: ATPAR3.I ATPAR4.I EGSCAL.I EVPDAT.I

calcum.o: FLAGS2.I MAXIMA.I MYRCAL.I MYRSAV.I

caltem.o: FLAGS3.I FLAGS4.I MAXIMA.I TEMPS.I

caludi.o: REGION.I

ccevrt.o: CEVBMY.I CITCIN.I EVAPAR.I EVAPGR.I FLAGS3.I FLAGS4.I INJECT.I \
         LEVBLK.I MAXIMA.I REGION.I TAMOUT.I TPDCOM.I VMXCOM.I

cevmpd.o: CEVBMY.I MAXIMA.I MYRCAL.I TPDCOM.I

cevtpd.o: CEVBMY.I MAXIMA.I TPDCOM.I

chkatp.o: ATPAR1.I ATPAR2.I FLAGS3.I

chkfun.o: PRGCH1.I PRGCH2.I PRSCH1.I PRSCH2.I

chkpgm.o: ATPAR2.I FLAGS2.I FLAGS3.I IMPAR1.I IMPAR6.I PRGCH2.I PRSCH2.I

consec.o: BASEQ9.I FLAGS1.I IOUCOM.I LEVINF.I PROJEC.I

crankc.o: CRACOM.I REGION.I

cvtpdf.o: CEVBMY.I MAXIMA.I MYRCAL.I

disatp.o: MYCODE.I TAMEQ4.I

discal.o: LOOKUP.I MYCODE.I TAMEQ1.I TAMEQ2.I TAMEQ3.I

disimp.o: CUMCOM.I LOOKUP.I MYCODE.I TAMEQ1.I TAMEQ2.I TAMEQ3.I TAMEQ4.I

diurnl.o: CITCIN.I CITRV1.I EVAPAR.I EVAPGR.I FLAGS3.I FLAGS4.I IOUCOM.I \
         MAXIMA.I TAMOUT.I TECHFL.I TEMPS.I

driver.o: COMMTS.I FLAGS3.I IOUCOM.I VDATA.I VNAME.I

efcalx.o: BYMYC1.I BYMYC2.I BYMYC5.I CITRV2.I CUMCOM.I FLAGS1.I FLAGS3.I \
         FLAGS4.I GSFCOM.I MAXIMA.I MYRCAL.I MYRSAV.I REGION.I RESUL1.I \
         RESUL3.I RVPEX2.I SPEED6.I VMXCOM.I

effgrp.o: EGSCAL.I LOOKUP.I MYCODE.I SIZCAL.I

efnftp.o: FLAGS1.I FLAGS4.I NFTPDT.I SCENE1.I

emirat.o: EGSCAL.I EMECAL.I LOOKUP.I MAXIMA.I

etpcal.o: ETPPAR.I FLAGS2.I

etpred.o: ETPPAR.I

ev81.o: EVADU1.I EVAHS1.I EVAPHD.I

ev85.o: EVAPHD.I

evhi.o: CITRV2.I EVADU2.I EVAHIA.I REGION.I

evmain.o: EVAHIA.I REGION.I

evpeme.o: EVPDAT.I IM240P.I LOOKUP.I MYCODE.I TAMEQ4.I TAMOU2.I TAMOUT.I

evpre.o: CITRV2.I EVADU1.I EVADU2.I EVAHS1.I REGION.I

evstd.o: CITRV2.I EVADU1.I EVADU2.I EVAHS1.I REGION.I

evunc.o: CITCIN.I CITRV2.I EVADU2.I REGION.I TECHFL.I TEMPS.I

exmptr.o: NMETH2.I REGION.I

fail.o: CUMCOM.I EVPDAT.I FLAGS3.I IM240P.I IMPAR5.I LOOKUP.I MAXIMA.I \
         MYCODE.I PPEFFD.I PRGCH1.I PRGCH2.I PRSCH1.I PRSCH2.I TAMOUT.I

findpb.o: LOOKUP.I MAXIMA.I MYCODE.I TAMEQ4.I TAMPB1.I

fuel.o: DTGGAS.I EGSCAL.I EVAPAR.I FLAGS1.I FLAGS2.I INJECT.I LOOKUP.I \
         MYCODE.I OFFSET.I OXY1.I OXY2.I RFORM1.I RFORM2.I

fuel2.o: FLAGS1.I FLAGS4.I OXY1.I RFORM1.I

getalh.o: ALUHIN.I FLAGS1.I FLAGS2.I IOUCOM.I MAXIMA.I

getatp.o: ATPAR1.I ATPAR2.I EGSCAL.I EVPDAT.I FLAGS1.I IOUCOM.I YEARS4.I

getbym.o: BYMYC1.I FLAGS1.I IOUCOM.I MAXIMA.I

getcol.o: FLAGS1.I IOUCOM.I TEMPC5.I TEMPC6.I YEARS4.I

getcum.o: CUMCOM.I MAXIMA.I MYRCAL.I

getefr.o: ATPAR4.I FLAGS1.I FLAGS2.I FLAGS3.I IMPAR2.I IMPAR5.I IMPAR6.I \
         IOUCOM.I PRGCH2.I PRSCH2.I

geterq.o: BASEQ3.I BASEQ4.I FLAGS1.I IOUCOM.I YEARS4.I

getetp.o: ETPPAR.I FLAGS1.I IOUCOM.I YEARS4.I

getfyr.o: FLAGS4.I IOUCOM.I NFTPDT.I YEARS4.I

gethor.o: FLAGS1.I FLAGS2.I FLAGS3.I FLAGS4.I IOUCOM.I NFTPDT.I PRGCH2.I \
         PRSCH2.I

getim2.o: FLAGS1.I IMPAR1.I IMPAR2.I IMPAR4.I IMPAR6.I IOUCOM.I TTCINF.I \
         YEARS4.I

getimc.o: FLAGS1.I FLAGS2.I ICR4VA.I ICR4VB.I IMPAR2.I IMPAR4.I IMPAR6.I \
         IMPAR7.I IOUCOM.I RSDINF.I TTCINF.I

getimp.o: FLAGS1.I FLAGS2.I IMPAR1.I IMPAR2.I IMPAR4.I IMPAR6.I IOUCOM.I \
         RSDINF.I TTCINF.I YEARS4.I

getiou.o: FLAGS1.I IOUCOM.I

getlev.o: BASEQ9.I BYMYC2.I FLAGS1.I FLAGS2.I FLAGS4.I IOUCOM.I LEVBLK.I \
         MAXIMA.I MYRCAL.I OFFSET.I OMTCOM.I

getmrn.o: FLAGS1.I FLAGS2.I IOUCOM.I MAXIMA.I MYRSAV.I NAMES1.I

getnat.o: FLAGS1.I FLAGS3.I IOUCOM.I RLCOM2.I YEARS4.I

getppd.o: IOUCOM.I PPEFFD.I PRGCH2.I PRSCH2.I

getprg.o: FLAGS1.I FLAGS2.I IMPAR1.I IMPAR2.I IOUCOM.I PRGCH1.I PRGCH2.I \
         YEARS4.I

getprs.o: FLAGS1.I IOUCOM.I PRSCH1.I PRSCH2.I YEARS4.I

getsc1.o: CITCIN.I EVAL.I FLAGS1.I FLAGS2.I IM12CO.I IM12HC.I IMPAR2.I \
         IMPAR4.I IMPAR6.I IMPAR7.I IOUCOM.I LEVBLK.I MAXIMA.I REGION.I \
         RSDINF.I SCENE1.I TEMPS.I TTCINF.I YEARS4.I

gettam.o: FLAGS1.I FLAGS2.I IOUCOM.I TAMEQ1.I

getttc.o: IOUCOM.I TTCINF.I

getver.o: FLAGS1.I FLAGS2.I FLAGS3.I FLAGS4.I IOUCOM.I NAMES1.I NFTPDT.I \
         PRGCH2.I PRSCH2.I STRING.I

getvmt.o: FLAGS1.I IOUCOM.I MAXIMA.I VMXCOM.I

gsfidx.o: FLAGS1.I GSFCOM.I MAXIMA.I REGISF.I

hccalx.o: BYMYC2.I BYMYC3.I BYMYC5.I CEVBMY.I CITCIN.I CITRV1.I CITRV2.I \
         CUMCOM.I EGSCAL.I EVAPGR.I FLAGS3.I FLAGS4.I GSFCOM.I LEVBLK.I \
         MAXIMA.I MYRCAL.I REGION.I RESUL1.I RESUL2.I RESUL3.I RLCOM1.I \
         RLCOM3.I RVPEX2.I SPEED6.I TEMPS.I VMXCOM.I

hdcgpm.o: BASEQ3.I

hdcpic.o: HDCCOM.I

hddmym.o: HDDMAR.I HDDREG.I MAXIMA.I MYRSAV.I

hdrcsv.o: IOUCOM.I

hotsok.o: CITCIN.I CITRV2.I EVAHS1.I EVAPAR.I EVAPGR.I EVAPHD.I FLAGS3.I \
         FLAGS4.I MAXIMA.I TAMOUT.I TEMPS.I

ibfptr.o: OPMOD2.I

iccopt.o: FLAGS4.I TEMPC5.I TEMPC7.I

idlcal.o: BYMYC2.I FLAGS4.I IVPCOM.I MAXIMA.I MYRCAL.I RESUL2.I RESUL3.I \
         SPEED4.I SPEED6.I VMXCOM.I

ierptr.o: BASEQ2.I BASEQ6.I FLAGS2.I REGION.I

ievptr.o: EVAPAR.I EVASTD.I INJECT.I

igeptr.o: ALUPAC.I REGION.I

iglptr.o: ALUPLT.I REGION.I

igrvpt.o: RVPEX2.I

igsfpt.o: GSFCOM.I REGION.I

igtptr.o: ALUPLT.I REGION.I

iguptr.o: ALUPAC.I REGION.I

imovlp.o: FLAGS2.I IMPAR1.I IMPAR2.I IMPAR6.I

imptr.o: IMPAR5.I

initap.o: EGSCAL.I

irulpt.o: RUNLS3.I

ispptr.o: REGION.I SPEED2.I SPEED4.I SPEED5.I

itampt.o: LOOKUP.I MYCODE.I MYUB1.I MYUB2.I TAMEQ2.I

itcptr.o: TEMPC2.I

levef.o: BASE11.I BASEQ9.I IOUCOM.I LEVBLK.I TAMEQ2.I

local.o: CITCIN.I CITRV1.I CITRV2.I EVAL.I FLAGS1.I FLAGS3.I OXY1.I \
         OXY2.I RFORM1.I RFORM2.I RVPNAT.I TEMPS.I

mobile.o: EVAL.I FLAGS2.I FLAGS3.I FLAGS4.I IOUCOM.I REGION.I SAVE01.I \
         YEARS4.I

myratp.o: FLAGS2.I IMPAR1.I IMPAR2.I

myrimt.o: FLAGS2.I IMPAR1.I IMPAR2.I

nparam.o: ATPAR2.I

obyhdr.o: CITPAR.I PROJEC.I

obyimp.o: BYMYC1.I BYMYC4.I BYMYC5.I IOUCOM.I MAXIMA.I MYRCAL.I

obytcn.o: BYMYC1.I BYMYC2.I BYMYC4.I EVAL.I IOUCOM.I LEVBLK.I MAXIMA.I \
         MYRCAL.I RESUL1.I SPEED6.I

onesec.o: FLAGS1.I FLAGS2.I FLAGS3.I FLAGS4.I IMPAR6.I IOUCOM.I NFTPDT.I \
         PRGCH2.I PRSCH2.I TAMPB1.I TTCINF.I

otcalc.o: EVPDAT.I TAMEQ4.I

outatp.o: ATPAR1.I ATPAR2.I FLAGS3.I IMPAR5.I IOUCOM.I PPEFFD.I PRGCH1.I \
         PRGCH2.I PRSCH1.I PRSCH2.I STRING.I

outby1.o: BYMYC1.I BYMYC2.I BYMYC3.I BYMYC4.I EVAL.I IOUCOM.I LEVBLK.I \
         MAXIMA.I MYRCAL.I RESUL1.I RESUL2.I SPEED6.I

outby2.o: BYMYC1.I BYMYC4.I IOUCOM.I MAXIMA.I

outby3.o: BYMYC1.I BYMYC4.I IOUCOM.I MAXIMA.I

outby4.o: BYMYC1.I BYMYC2.I BYMYC3.I BYMYC4.I EVAL.I IOUCOM.I LEVBLK.I \
         MAXIMA.I MYRCAL.I RESUL1.I RESUL2.I SPEED6.I

outcsv.o: ATPAR1.I ATPAR2.I BASE10.I CITPAR.I CITRV1.I ETPPAR.I EVAL.I \
         FLAGS1.I FLAGS2.I FLAGS3.I FLAGS4.I IMPAR1.I IMPAR2.I IMPAR4.I \
         IMPAR6.I IOUCOM.I LEVBLK.I MAXIMA.I OXY1.I PRGCH1.I PRGCH2.I \
         PROJEC.I PRSCH1.I PRSCH2.I REGION.I RESUL1.I RESUL2.I RESUL3.I \
         RFORM1.I RLCOM2.I RUNLS1.I SCENE1.I STRING.I TEMPS.I VMXCOM.I \
         YEARS4.I

outdt1.o: CITPAR.I FLAGS3.I FLAGS4.I IOUCOM.I MAXIMA.I REGION.I RESUL1.I \
         RESUL2.I RESUL3.I SCENE1.I TEMPS.I VMXCOM.I

outdt2.o: FLAGS4.I IOUCOM.I MAXIMA.I REGION.I RESUL1.I RESUL2.I RESUL3.I \
         SCENE1.I TEMPS.I VMXCOM.I

outdt3.o: ALUHIN.I BASE10.I CITCIN.I CITRV1.I EVAPGR.I FLAGS1.I FLAGS2.I \
         FLAGS3.I FLAGS4.I IOUCOM.I LEVBLK.I MAXIMA.I OXY1.I OXY2.I \
         REGION.I RESUL1.I RESUL2.I RESUL3.I RFORM1.I SCENE1.I STRING.I \
         TEMPC5.I TEMPC6.I TEMPS.I USDATA.I VMXCOM.I

outdt4.o: ALUHIN.I BASE10.I CITCIN.I FLAGS1.I FLAGS2.I FLAGS3.I FLAGS4.I \
         IOUCOM.I LEVBLK.I MAXIMA.I OXY1.I OXY2.I REGION.I RESUL1.I \
         RESUL2.I RESUL3.I RFORM1.I SCENE1.I STRING.I TEMPS.I USDATA.I \
         VMXCOM.I

outdt5.o: EVAL.I FLAGS4.I IOUCOM.I

outdtg.o: DTGGAS.I IOUCOM.I

outetp.o: ETPPAR.I IOUCOM.I

outgsf.o: IOUCOM.I MAXIMA.I REGISF.I

outhd1.o: CITRV1.I FLAGS2.I FLAGS3.I IOUCOM.I MAXIMA.I OXY1.I REGION.I \
         RFORM1.I STRING.I TEMPS.I

outhd2.o: IOUCOM.I MAXIMA.I REGION.I

outhd3.o: ETPPAR.I FLAGS1.I FLAGS2.I FLAGS3.I IOUCOM.I PRGCH2.I PRSCH2.I \
         REGION.I

outhd4.o: ETPPAR.I FLAGS1.I FLAGS2.I FLAGS3.I IOUCOM.I PRGCH2.I PRSCH2.I \
         REGION.I

outhdr.o: PROJEC.I VDATA.I

outimp.o: FLAGS2.I IMPAR1.I IMPAR2.I IMPAR4.I IMPAR5.I IMPAR6.I IOUCOM.I \
         STRING.I TTCINF.I

outlap.o: CITPAR.I CITRV1.I FLAGS3.I IOUCOM.I RFORM1.I TEMPS.I

outmth.o: FLAGS3.I FLAGS4.I IOUCOM.I

outnew.o: BASEQ3.I BASEQ4.I IOUCOM.I STRING.I

outpol.o: FLAGS4.I MAXIMA.I

output.o: EVAL.I FLAGS1.I FLAGS3.I IOUCOM.I LEVBLK.I LEVINF.I MYRCAL.I \
         REGION.I YEARS4.I

outstg.o: IOUCOM.I RLCOM2.I

outtam.o: FLAGS2.I FLAGS3.I IOUCOM.I TAMEQ1.I

outzev.o: BASE10.I IOUCOM.I LEVBLK.I MAXIMA.I YEARS4.I

outztt.o: BASE10.I IOUCOM.I LEVBLK.I MAXIMA.I MYRCAL.I

ovlp2.o: FLAGS2.I IMPAR1.I IMPAR2.I IMPAR6.I

p3wpox.o: EGSCAL.I LOOKUP.I

parsec.o: FLAGS1.I FLAGS4.I

pcleft.o: FLAGS2.I ICR4VA.I ICR4VB.I IM12CO.I IM12HC.I IMPAR1.I IMPAR2.I \
         IMPAR5.I IMPAR6.I MAXIMA.I TTCINF.I

pctlev.o: BASEQ9.I LEVBLK.I

piv.o: TEMPC5.I TEMPC6.I TEMPC7.I

quiter.o: ATPAR2.I BASEQ4.I BASEQ6.I FLAGS2.I FLAGS3.I IMPAR1.I IMPAR6.I \
         IOUCOM.I NAMES1.I PRGCH2.I PRSCH2.I QUITXQ.I YEARS4.I

rdcity.o: CITPAR.I CITRV1.I DTGGAS.I FLAGS1.I FLAGS2.I FLAGS3.I IOUCOM.I \
         MAXIMA.I OXY1.I REGISF.I RFORM1.I RFORM2.I RVPNAT.I TEMPS.I \
         YEARS4.I

refuel.o: CITRV1.I FLAGS2.I FLAGS3.I MAXIMA.I RLCOM1.I RLCOM2.I RLCOM3.I \
         TAMOU2.I TEMPS.I

regmod.o: FLAGS2.I FLAGS3.I GSFCOM.I MAXIMA.I MYRCAL.I MYRSAV.I REGION.I \
         REGISF.I

reinit.o: ATPAR1.I ATPAR2.I BASEQ6.I BYMYC2.I BYMYC3.I BYMYC5.I ETPPAR.I \
         EVPDAT.I FLAGS2.I IM240P.I IMPAR1.I IMPAR2.I IMPAR5.I IMPAR6.I \
         MAXIMA.I MYRSAV.I OFFSET.I PARAME.I PPEFFD.I QUITXQ.I REGION.I \
         RUNLS1.I TAMOU2.I TAMOUT.I TAMPB1.I

restor.o: ETPPAR.I IMPAR4.I MAXIMA.I MYRSAV.I SAVE01.I SAVE04.I SAVE05.I \
         SAVE06.I TAMEQ1.I VMXCOM.I

rstlos.o: FLAGS3.I RESTLS.I TEMPS.I

ruloss.o: CITCIN.I FLAGS3.I RUNLS1.I RUNLS2.I RUNLS3.I SCENE1.I TAMOUT.I \
         TEMPS.I

rvpexh.o: CITRV2.I MAXIMA.I RVPEX1.I RVPEX2.I TEMPS.I

savepb.o: ATPAR1.I IMPAR1.I IMPAR2.I LOOKUP.I MAXIMA.I MYCODE.I TAMEQ4.I \
         TAMPB1.I YEARS4.I

saver.o: EVAPGR.I MAXIMA.I RESUL1.I RESUL2.I RESUL3.I VMXCOM.I X1.I X2.I \
         X3.I X4.I X5.I

scfsub.o: MAXIMA.I SPEED1.I SPEED3.I SPEED4.I

setegs.o: EGSCAL.I LOOKUP.I SIZCAL.I

stoerq.o: BASEQ4.I BASEQ5.I BASEQ6.I

tamper.o: ATPAR1.I ATPAR2.I EGSCAL.I FLAGS3.I FLAGS4.I IM240P.I IMPAR1.I \
         IMPAR2.I IMPAR6.I LOOKUP.I MAXIMA.I MYCODE.I TAMEQ4.I TAMOUT.I \
         TAMPB1.I YEARS4.I

tech0.o: LOOKUP.I MYCODE.I SIZCAL.I TAMEQ4.I

tfcalx.o: FLAGS1.I MAXIMA.I MYRCAL.I REGISF.I VMXCOM.I

togfid.o: FID.I OXY1.I OXY2.I

usdptr.o: FLAGS1.I FLAGS2.I USDATA.I

vocfid.o: FID.I OXY1.I OXY2.I

wtlrul.o: FLAGS1.I IOUCOM.I RUNLS1.I

yrtest.o: YEARS4.I

