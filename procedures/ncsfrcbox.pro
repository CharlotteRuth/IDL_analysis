PRO NCSFRCBOX,s, meansfr, BINSIZE=binsize, OVERPLOT=overplot,COLOR=color,NODATA=nodata,move=move,NOPLOT=noplot,YMAX=ymax,YMIN=ymin,LINESTYLE=linestyle, FILENAME=filename,MASSUNIT=massunit,BOXSIZE=boxsize,_EXTRA=_extra
;plots the star formation rate for a simulation

timeunit=3.8e10
IF (keyword_set(massunit) EQ 0) THEN massunit=4.7526e16
IF (keyword_set(boxsize) EQ 0) THEN boxsize=68.493  ; =50/h Mpc
zmax=10

IF (keyword_set(s) EQ 0) THEN BEGIN
  rncarray,'timeform',h,tform
  s = replicate({tform: 1., mass: 1.},h.n)
  s.tform=tform
  rncarray,'massform',h,m
  s.mass=m
ENDIF

IF (keyword_set(binsize) EQ 0) THEN binsize=5.e7

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
IF (KEYWORD_SET(ymin) EQ 0) THEN ymin=sfr[where(midtimearray GT zmax)]
 
IF (KEYWORD_SET(noplot) EQ 0) THEN BEGIN
      IF KEYWORD_SET(nodata) THEN plot, midtimearray,sfr,psym=10,/nodata,yrange=[ymin,ymax],xrange=[min(midtimearray),6],xtitle="Time [yrs]",ytitle='SFR [M!lsolar!n/yr]',/ylog,xstyle=1,ystyle=1,_EXTRA=_extra
ENDIF
    
;Bouwens 2006
giavz=[0.25,0.5,0.75,1,2,3,3.75,4,5,6,7.5,10]
giavalisco=[10^(-3.45),10^(-2.85),10^(-2.55),10^(-2.5),10^(-1.95),10^(-1.75),10^(-1.8),10^(-1.8),10^(-1.95),10^(-2.25),10^(-3.),10^(-4.)]
giaverr=[10^(-0.7)-10^(-0.75),10^(-0.7)-10^(-0.8),10^(-0.75)-10^(-0.85),10^(-0.85)-10^(-0.95),10^(-0.7)-10^(-0.9)]
IF KEYWORD_SET(overplot) THEN oplot, midtimearray,sfr,psym=10,color=color,thick=2,linestyle=linestyle ELSE BEGIN
  plot,midtimearray,sfr,psym=10,color=color,thick=2,yrange=[ymin,ymax],xrange=[min(midtimearray)-1.0,zmax],xtitle="z",ytitle='SFR [M!lsolar!n/yr/Mpc!u3!n]',linestyle=linestyle,/ylog,xstyle=1,ystyle=1,_EXTRA=_extra
  ;oploterr,giavz,giavalisco,giaverr,6
  oplot,giavz,giavalisco,psym=6
ENDELSE

END
