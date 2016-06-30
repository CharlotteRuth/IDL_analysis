PRO SFRCOSMOBOX,s, binsize, meansfr,OVERPLOT=overplot,COLOR=color,NODATA=nodata,move=move,NOPLOT=noplot,YMAX=ymax,YMIN=ymin,LINESTYLE=linestyle, FILENAME=filename,MASSUNIT=massunit,BOXSIZE=boxsize
;plots the star formation rate for a simulation

timeunit=4.0435875e10
IF (keyword_set(massunit) EQ 0) THEN massunit=4.95469e13
IF (keyword_set(boxsize) EQ 0) THEN boxsize=7.14  ; =5/h Mpc
zmax=10

IF (keyword_set(s) EQ 0) THEN rtipsy,filename,h,g,d,s

IF (keyword_set(binsize) EQ 0) THEN binsize=1.e7

ind=WHERE(s.tform GT 0.0)
if ( keyword_set(move) EQ 0 ) then tform=s[ind].tform*timeunit else tform=s[ind].tform*timeunit + move
mass=s[ind].mass*massunit
mintime=MIN(tform)
maxtime=MAX(tform)
nbins=FIX((maxtime-mintime)/binsize)+1
timearray=FINDGEN(nbins)*binsize+mintime
midtimearray=FINDGEN(nbins)*binsize+mintime+binsize*0.5
sfr=FLTARR(nbins)

FOR i=0,nbins-1 DO BEGIN
    inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize),ninbin)
    IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/(boxsize^3)
ENDFOR

midtimearray=time2z(midtimearray)

IF (KEYWORD_SET(ymax) EQ 0) THEN ymax=MAX(sfr)
IF (KEYWORD_SET(ymin) EQ 0) THEN ymin=sfr[where(midtimearray EQ zmax)]
 
IF (KEYWORD_SET(noplot) EQ 0) THEN BEGIN
      IF KEYWORD_SET(nodata) THEN plot, midtimearray,sfr,psym=10,/nodata,yrange=[ymin,ymax],xrange=[min(midtimearray),6],xtitle="Time [yrs]",ytitle='SFR [M!lsolar!n/yr]',/ylog,xstyle=1,ystyle=1
ENDIF
    
giavz=[3,3.75,4.03,4.8,5.7]
giavalisco=[10^(-0.75),10^(-0.8),10^(-0.85),10^(-0.95),10^(-0.9)]
giaverr=[10^(-0.7)-10^(-0.75),10^(-0.7)-10^(-0.8),10^(-0.75)-10^(-0.85),10^(-0.85)-10^(-0.95),10^(-0.7)-10^(-0.9)]
IF KEYWORD_SET(overplot) THEN oplot, midtimearray,sfr,psym=10,color=color,thick=2,linestyle=linestyle ELSE BEGIN
  plot,midtimearray,sfr,psym=10,color=color,thick=2,yrange=[ymin,ymax],xrange=[min(midtimearray)-0.2,zmax],xtitle="z",ytitle='SFR [M!lsolar!n/yr/Mpc!u3!n]',linestyle=linestyle,/ylog,xstyle=1,ystyle=1
  oploterr,giavz,giavalisco,giaverr,6
ENDELSE

END
