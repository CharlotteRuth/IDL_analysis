
pro compSFR_twins,files,time,m,colors,linestyle = linestyle, key = key,xtitle = xtitle, ytitle = ytitle,print=print, outplot = outplot


files =['/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.00512/h603.cosmo50cmb.2304g2bwK.00512.1.std','/astro/net/scratch1/abrooks/FABIO/h516.cosmo25cmb.2304g2bwK/h516.cosmo25cmb.2304g2bwK.00512/h516.cosmo25cmb.2304g2bwK.00512.1.std']
time = [15,15]
colors = [50,240]
m = [1.,1.]
key = ['h603','h516']
outroot = 'twinsSFH2'
massunit = [1.84793e16,2.310e15]
timeunit = [3.8785614e+10,3.8784879e+10]
linestyle = [0,0]
ytitle = textoidl('SFR/M_b [yr^(-1)]')

;This graphs the (star fromation rate)/(gass mass) vs time in
;different ways and is taken from compSFR_spec.pro

nfiles = N_ELEMENTS(files)

binsize=5.e7
time = time*timeunit
integrate = dblarr(nfiles)
lines = linestyle
ymax=0 
FOR fct=0,nfiles -1 DO BEGIN
    first = 0
    start = 0

    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
    file = files[fct]
    rtipsy,file,h,g,d,s
    IF (s[0].mass NE 0)THEN BEGIN
        ind=WHERE(s.tform GT 0.0)
        tform=s[ind].tform*timeunit[fct]
        mass=s[ind].mass*massunit[fct]
        mintime=MIN(tform)
        maxtime=MAX(tform)
        nbins=FIX((maxtime-mintime)/binsize)+1
        timearray=FINDGEN(nbins)*binsize+mintime
        sfr=FLTARR(nbins)
        m[fct] = (total(g.mass) + total(s.mass))*massunit[fct]
        FOR i=0,nbins-1 DO BEGIN
            inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize),ninbin)
            IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/m[fct] ;massin bin/time for each bin/total mass of galaxy/fraction that is gas
        ENDFOR
        ymaxtemp=MAX(sfr)
        IF (ymaxtemp GT ymax)THEN ymax=ymaxtemp
        integrate[fct]=TOTAL(mass)/time[fct]
        print,files[fct]
        print,'Average SFR [yr^-1]: ',integrate[fct]
        print,'Average SFR [yr^-1]/Mbaryon: ',integrate[fct]/m[fct]
        print,'Average SFR [yr^-1]/Mstar: ',integrate[fct]/total(s.mass)/massunit[fct]        
        print,'Range of SFR: ',MINMAX(sfr)*m[fct]
        print,'Range of SFR: ',MINMAX(sfr)
        print,'Baryonic Mass',m[fct]
;        print,'Baryonic Mass',(total(g.mass) + total(s.mass))*massunit[fct]
        print,'Total Mass   ',(total(g.mass) + total(s.mass) + total(d.mass))*massunit[fct]
    ENDIF ELSE integrate[fct]=0.0000
ENDFOR
ymax=1.01*ymax
set_plot,'x'
loadct,39
set_plot,'ps'
device,filename='/astro/net/scratch1/christensen/Twins/'+outroot+'.eps',/color,bits_per_pixel=8,/times
FOR fct=0,nfiles - 1 DO BEGIN
    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
    file = files[fct]
    rtipsy,file,h,g,d,s
    IF (s[0].mass NE 0)THEN BEGIN
        IF first EQ 0 THEN BEGIN
            first = 1
            sfr,s,massu=massunit[fct],time=timeunit[fct],OVERPLOT=0,xrange=[0,MAX(time/timeunit[fct])],binsize=binsize,gmass=m[fct],xtitle = xtitle, ytitle = ytitle, linestyle = lines[fct],yrange = [0,ymax]
        ENDIF ELSE BEGIN
            first = first+1
        ENDELSE
        sfr,s,massu=massunit[fct],time=timeunit[fct],COLOR=COLORS[fct],OVERPLOT=1,binsize=binsize,gmass=m[fct],linestyle = lines[fct]
    ENDIF ELSE start = start+1
ENDFOR
IF (KEYWORD_SET(key)) THEN legend,key,linestyle = lines, color = colors,/right ELSE BEGIN
    IF (nfiles eq 1) THEN legend,[key[0]+' '+trim(integrate[0])],linestyle=[0],color=[colors[0]],/right
    IF (nfiles eq 2) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1])],linestyle=[0,0],color=[colors[0],colors[1]],/right
    IF (nfiles eq 3) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2])],linestyle=[0,0,0],color=[colors[0],colors[1],colors[2]],/right  
ENDELSE
device,/close

END
