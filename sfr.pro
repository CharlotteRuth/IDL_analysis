PRO SFR, s, meansfr,binsize=binsize, OVERPLOT=overplot,move=move,NOPLOT=noplot,YMAX=ymax,LINESTYLE=linestyle,MINT=mint,maxt=maxtime,MASSUNIT=massunit,TIMEUNIT=timeunit,MASSFORM=massform,sfh=sfr,_EXTRA=_extra, gMass=gMass,xrange=xrange
;plots the star formation rate for a simulation

IF (keyword_set (timeunit) EQ 0) THEN timeunit=4.0435875e10
IF (keyword_set(massunit) EQ 0) THEN massunit=13.6e16
IF (keyword_set(gMass) EQ 0) THEN gMass=1

IF (keyword_set(binsize) EQ 0) THEN binsize=1.e7

ind=WHERE(s.tform GT 0.0)
if (keyword_set(move) EQ 0 ) then tform=s[ind].tform*timeunit else tform=s[ind].tform*timeunit + move
IF (keyword_set(massform) EQ 0) THEN mass=s[ind].mass*massunit $
ELSE mass=massform[ind]*massunit
IF (keyword_set(mint) EQ 0) THEN mintime=MIN(tform) ELSE mintime=mint
maxtime=MAX(tform)
nbins=FIX((maxtime-mintime)/binsize)+1
timearray=FINDGEN(nbins)*binsize+mintime
midtimearray=FINDGEN(nbins)*binsize+mintime+binsize*0.5
sfr=FLTARR(nbins)

FOR i=0,nbins-1 DO BEGIN
    inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize),ninbin)
    IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize
ENDFOR
IF (KEYWORD_SET(ymax) EQ 0) THEN ymax=MAX(sfr)
 
IF (KEYWORD_SET(noplot) EQ 0) THEN BEGIN
ENDIF

stop
;IF (KEYWORD_SET(overplot)EQ 0) THEN BEGIN plot, midtimearray/1e9,sfr/gMass,psym=10,color=color,linestyle=linestyle,_EXTRA=_extra ELSE plot,midtimearray/1e9,sfr,psym=10,color=color,yrange=[0,ymax],xtitle="Time [Gyr]",ytitle='SFR [M'+sunsymbol()+' yr!u-1!n]',ytitle='SFR/Gass Mass 10e'+gMass ' M'+sunsymbol(), linestyle=linestyle,_EXTRA=_extra 
   
 
;IF KEYWORD_SET(overplot) THEN oplot, midtimearray/1e9,sfr/10.^gMass*10.,psym=10,color=color,linestyle=linestyle,_EXTRA=_extra ELSE plot,midtimearray/1e9,sfr/10.^gMass*10.,psym=10,color=color,yrange=[0,ymax],xtitle="Time [Gyr]",ytitle='SFR/Gas Mass',linestyle=linestyle,_EXTRA=_extra, xrange=xrange,title='SFR/Gas Mass -- Dark Matter Mass: 10^'+trim(gMass)+sunsymbol()

;print,ymax
;ytitle='SFR [M]+sunsymbol()+' yr!u-1!n]

IF N_PARAMS() EQ 2 THEN BEGIN
    mint=mintime+5.e8
    maxt=mintime+15.e8
    ind=WHERE(midtimearray GE mint AND midtimearray LE maxt)
    meansfr=MEAN(sfr[ind])
    IF (KEYWORD_SET(noplot) EQ 0) THEN BEGIN
       ;if ( ymax GT max(sfr)) then dashmax=ymax
       plots,[mint/1e9,mint/1e9],[0,ymax],linestyle=2
       plots,[maxt/1e9,maxt/1e9],[0,ymax],linestyle=2
    ENDIF
ENDIF

END
