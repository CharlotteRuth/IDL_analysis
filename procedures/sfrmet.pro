PRO SFRMET, s,meansfr, BINTSIZE=tbinsize, BINMSIZE=mbinsize,move=move,NOPLOT=noplot,YMAX=ymax,MINT=mint,MINMET=minmet,MAXMET=maxmet,MASSUNIT=massunit,_EXTRA=_extra
;plots the star formation rate for a simulation

timeunit=4.0435875e10
IF (keyword_set(massunit) EQ 0) THEN massunit=3.171e15

IF (keyword_set(tbinsize) EQ 0) THEN tbinsize=0.1
IF (keyword_set(mbinsize) EQ 0) THEN mbinsize=1.e-3

ind=WHERE(s.tform GT 0.0)
if ( keyword_set(move) EQ 0 ) then tform=s[ind].tform*timeunit/1.e9 else tform=s[ind].tform*timeunit + move
mass=s[ind].mass*massunit
IF (keyword_set(mint) EQ 0) THEN mintime=MIN(tform) ELSE mintime=mint
maxtime=MAX(tform)
;ntbins=FIX((maxtime-mintime)/tbinsize)+1
ntbins=FIX((maxtime)/tbinsize)+1
timearray=FINDGEN(ntbins)*tbinsize+mintime
midtimearray=FINDGEN(ntbins)*tbinsize+tbinsize*0.5;+mintime

IF (keyword_set(minmet) EQ 0) THEN minmet=MIN(s.metals)
IF (keyword_set(maxmet) EQ 0) THEN maxmet=MAX(s.metals)
nzbins=FIX((maxmet-minmet)/mbinsize)+1
metarray=FINDGEN(nzbins)*mbinsize+minmet
midmetarray=FINDGEN(nzbins)*mbinsize+mbinsize*0.5;+minmet

;minz=floor(min(time2z(midtimearray)))
;if (minz LT 0) then minz = 0
;zs=findgen(10-minz)+minz
;ztimevals=1.345773e10-lookback(zs)
;sfr=FLTARR(ntbins,nzbins)
sfr=hist_2d(s.metals,tform,bin2=tbinsize,bin1=mbinsize,max1=maxmet,min1=minmet)

;FOR i=0,ntbins-1 DO BEGIN
;  FOR j=0,nzbins-1 DO BEGIN
;    inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+tbinsize) AND s.metals GE metarray[j] AND s.metals LT (metarray[j]+mbinsize),ninbin)
;    IF (ninbin EQ 0) THEN sfr[i,j]=0 ELSE sfr[i,j]=TOTAL(mass[inbin])/tbinsize
;  ENDFOR
;ENDFOR

IF (KEYWORD_SET(ymax) EQ 0) THEN ymax=MAX(sfr)
 
surface,sfr,midmetarray,midtimearray,zrange=[0,ymax],ytitle="Time [Gyrs]",xtit="Metal Fraction",ztitle='SFR [M!lsolar!n/yr]',_EXTRA=_extra,az=110,zaxis=2
;axis,/xaxis,xtickv=ztimevals,xtickname=strtrim(string(zs,FORMAT='(I)'),2),xtit='z',xticks=n_elements(zs)-1,xticklayout=0

END
