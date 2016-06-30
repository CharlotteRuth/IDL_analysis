pro SFR, s, meansfr,binsize=binsize, OVERPLOT=overplot,move=move,NOPLOT=noplot,YMAX=ymax,LINESTYLE=linestyle,MINTime=mintime,maxt=maxtime,MASSUNIT=massunit,TIMEUNIT=timeunit,MASSFORM=massform,sfh=sfr,_EXTRA=_extra, gMass=gMass,xrange=xrange,ntrial = ntrial,color=color,title=title,tarray = tarray,sarray = sarray,STARLOG = starlog,cumlative = cumlative,redshift = redshift
;plots the star formation rate for a simulation

IF KEYWORD_SET(STARLOG) THEN BEGIN
    s1 = s
    s2_base = {mass:0.0, x:0.0, y:0.0, z:0.0, vx:0.0, vy:0.0, vz:0.0, metals:0.0, tform:0.0, eps:0.0, phi:0.0}
    s2 = replicate(s2_base,N_ELEMENTS(s))
    FOR i = 0L, N_ELEMENTS(s1) - 1 DO BEGIN
       s2[i].mass=s1[i].massform
       s2[i].x=s1[i].x
       s2[i].y=s1[i].y
       s2[i].z=s1[i].z
       s2[i].vx=s1[i].vx
       s2[i].vy=s1[i].vy
       s2[i].vz=s1[i].vz
       s2[i].tform=s1[i].timeform   
   ENDFOR
   s = s2
ENDIF

IF (max(s.mass) eq 0) then BEGIN
    IF (keyword_set(ymax) EQ 0) THEN ymax = 1
    IF KEYWORD_SET(overplot) THEN oplot,[0,0],[0,0],psym=10,color=color,linestyle=linestyle,_EXTRA=_extra ELSE plot,[0,0],[0,0] ,psym=10,color=color,xtitle="Time [Gyr]",ytitle='SFR [M'+sunsymbol()+' yr!u-1!n]',linestyle=linestyle,_EXTRA=_extra, xrange=xrange,yrange=[0,ymax],title=title
    return
ENDIF

IF (keyword_set (timeunit) EQ 0) THEN timeunit=4.0435875e10
IF (keyword_set(massunit) EQ 0) THEN massunit=2.31000e+15 ;13.6e16
IF (keyword_set(gMass) EQ 0) THEN gMass=1
IF (keyword_set(ntrial) EQ 0) THEN ntrial = 1.

IF (keyword_set(binsize) EQ 0) THEN binsize=1.e7

ind=WHERE(s.tform GT 0.0)
if (keyword_set(move) EQ 0 ) then tform=s[ind].tform*timeunit else tform=s[ind].tform*timeunit + move
IF (keyword_set(massform) EQ 0) THEN mass=s[ind].mass*massunit $
ELSE mass=massform[ind]*massunit
IF (keyword_set(mintime) EQ 0) THEN mintime=MIN(tform)*timeunit
IF (keyword_set(maxtime) EQ 0) THEN maxtime=MAX(tform)*timeunit
IF NOT KEYWORD_SET(xrange) THEN xrange = [mintime,maxtime]/1e9

nbins=FIX((maxtime-mintime)/binsize)+1
timearray=FINDGEN(nbins)*binsize+mintime
midtimearray=FINDGEN(nbins)*binsize+mintime+binsize*0.5
sfr=FLTARR(nbins)
stop
FOR i=0,nbins-1 DO BEGIN
    IF KEYWORD_SET(cumlative) THEN inbin=WHERE(tform LE timearray[i]+binsize,ninbin) ELSE inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize),ninbin)
    IF (ninbin EQ 0) THEN sfr[i]=0 ELSE begin
                                ;print,TOTAL(mass[inbin])
        IF keyword_set(cumlative) THEN sfr[i] = TOTAL(mass[inbin])/ntrial ELSE sfr[i]=TOTAL(mass[inbin])/binsize/ntrial
    endelse
ENDFOR

IF keyword_set(cumlative) THEN ytitle = 'Stellar Mass [M'+sunsymbol() + ']' ELSE ytitle = 'SFR [M'+sunsymbol()+' yr!u-1!n]'

IF (KEYWORD_SET(ymax) EQ 0) THEN ymax=MAX(sfr)/gmass
 
IF (KEYWORD_SET(noplot) EQ 0) THEN BEGIN
ENDIF

;IF (KEYWORD_SET(overplot)EQ 0) THEN BEGIN plot, midtimearray/1e9,sfr/gMass,psym=10,color=color,linestyle=linestyle,_EXTRA=_extra ELSE plot,midtimearray/1e9,sfr,psym=10,color=color,yrange=[0,ymax],xtitle="Time [Gyr]",ytitle='SFR [M'+sunsymbol()+' yr!u-1!n]',ytitle='SFR/Gasxs Mass 10e'+gMass ' M'+sunsymbol(), linestyle=linestyle,_EXTRA=_extra 
   
stop 
IF keyword_set(cumlative) THEN xarray = timearray+binsize ELSE xarray = midtimearray
IF KEYWORD_SET(overplot) THEN oplot, xarray/1e9,sfr/gMass,psym=10,color=color,linestyle=linestyle,_EXTRA=_extra $
ELSE BEGIN
    IF NOT KEYWORD_SET(redshift) THEN plot,xarray/1e9,sfr/gmass,psym=10,color=color,xtitle="Time [Gyr]",ytitle=ytitle,linestyle=linestyle,_EXTRA=_extra, xrange=xrange,yrange=[0,ymax],title=title $
    ELSE BEGIN
        ageUniverse = 13.7346*1e9 ;wmap3_lookback(100)
        tickred_in_t = REVERSE(findgen(5))
        ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e9
        plot,xarray/1e9,sfr/gmass,psym=10,color=color,xtitle="Time [Gyr]",ytitle=ytitle,linestyle=linestyle,yrange=[0,ymax],title=title,xstyle = 9,xrange = xrange,_EXTRA=_extra
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
     ENDELSE
 ENDELSE

;print,MAX(sfr/gmass)
;print,gmass
;print,ymax
;ytitle='SFR [M]+sunsymbol()+' yr!u-1!n]
IF N_PARAMS() EQ 2 THEN BEGIN
    mint=mintime+5.e8
    maxt=mintime+15.e8
    ind=WHERE(midtimearray GE mint AND midtimearray LE maxt)
    meansfr=MEAN(sfr[ind])
    print,'MeanSFR: ',meansfr
    IF (KEYWORD_SET(noplot) EQ 0) THEN BEGIN
        if ( ymax GT max(sfr)) then dashmax=ymax
       plots,[mint/1e9,mint/1e9],[0,ymax],linestyle=2
       plots,[maxt/1e9,maxt/1e9],[0,ymax],linestyle=2
   ENDIF
ENDIF
sarray = sfr/gMass
tarray = xarray
return

END
