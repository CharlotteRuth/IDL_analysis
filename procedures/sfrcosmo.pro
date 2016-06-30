PRO SFRCOSMO, s,meansfr, BINSIZE=binsize, OVERPLOT=overplot,move=move,NOPLOT=noplot,YMAX=ymax,MINT=mint,MASSUNIT=massunit,MAXTIME=maxtime,HALFTIME=halftime,AVETIME=avetime,SANDTAU=sandtau,TOFFSET=toffset,BIRTHRATE=birthrate,_EXTRA=_extra
;plots the star formation rate for a simulation

timeunit=4.0435875e10
IF (keyword_set(massunit) EQ 0) THEN massunit=3.171e15

IF (keyword_set(binsize) EQ 0) THEN binsize=1.e7

ind=WHERE(s.tform GT 0.0)
IF ( keyword_set(move) EQ 0 ) THEN tform=s[ind].tform*timeunit ELSE tform=s[ind].tform*timeunit + move
mass=s[ind].mass*massunit
IF (keyword_set(mint) EQ 0) THEN mintime=MIN(tform) ELSE mintime=mint
maxtime=MAX(tform)
nbins=FIX((maxtime-mintime)/binsize)+1
timearray=FINDGEN(nbins)*binsize+mintime
midtimearray=FINDGEN(nbins)*binsize+mintime+binsize*0.5
minz=floor(min(time2z(midtimearray)))
if (minz LT 0) then minz = 0
if (minz GT 5) then minz = 5
zs=findgen(6-minz)+minz
ztimevals=1.345773e10-lookback(zs)
sfr=FLTARR(nbins)

FOR i=0,nbins-1 DO BEGIN
    inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize),ninbin)
    IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize
ENDFOR

IF (KEYWORD_SET(ymax) EQ 0) THEN ymax=MAX(sfr)
 
IF KEYWORD_SET(overplot) THEN oplot, midtimearray/1e9,sfr,psym=10,_EXTRA=_extra ELSE BEGIN
IF KEYWORD_SET(birthrate) THEN BEGIN
  meansfr=mean(sfr)
  ymax = ymax/meansfr
   plot,midtimearray/1e9,sfr/meansfr,psym=10,yrange=[0,ymax],xtitle="Time [Gyrs]",ytitle='b [SFR/<SFR>]',xstyle=8,_EXTRA=_extra 
ENDIF ELSE  plot,midtimearray/1e9,sfr,psym=10,yrange=[0,ymax],xtitle="Time [Gyrs]",ytitle='SFR [M!lsolar!n/yr]',xstyle=8,_EXTRA=_extra

   axis,/xaxis,xtickv=ztimevals/1e9,xtickname=strtrim(string(zs,FORMAT='(I)'),2),xtit='z',xticks=n_elements(zs)-1,xticklayout=0
ENDELSE
;IF KEYWORD_SET(maxtime) THEN BEGIN
  firstguess=timearray[where(sfr EQ ymax)]
  maxtime=firstguess[0]
  plots,[maxtime,maxtime],[0,MIN([MAX(sfr),ymax])],linestyle=2
;ENDIF
;IF KEYWORD_SET(halftime) THEN BEGIN
  halfmass = 0.5*total(mass)
  i=0
  while(total(sfr[0:i]*binsize) LT halfmass) DO i=i+1
  halftime=timearray[i]
  plots,[halftime,halftime],[0,MIN([MAX(sfr),ymax])],linestyle=2
  avetime=mean(tform)
;ENDIF
;IF KEYWORD_SET(sandagetau) THEN BEGIN
  ;A = [1.0e4,2.0e9,1e9]
if (n_elements(sfr) GT 2)  then begin
  A = [1.0e4,2.0e9]
  weights=fltarr(n_elements(midtimearray))+1.0
  sandfit = curvefit(midtimearray,sfr,weights,A,SIGMA,FUNCTION_NAME='sandagetau')
;  oplot,midtimearray/1e9,sandfit
  sandtau=A[1]
endif
  ;toffset=A[2]
;ENDIF

END
