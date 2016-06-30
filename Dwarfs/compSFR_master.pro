PRO compSFR_master
loadct,39
!P.CHARSIZE = 1.25
!P.thick = 1.5
!X.Charsize = 1
!Y.Charsize = 1
!Y.Style = 1
!X.style = 1

base = ['/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/','/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/']
imf = ['k'];,'ms']
legend_name0 = [textoidl('10^6'),textoidl('10^5'),textoidl('10^4'),textoidl('10^3')]
res_str0 = ['1E6R','1E5R','1E4R','1E3R']
colors0 = [0,0,0,0]
linestyle0 = [0,0,2,1]
legend_name1 = [textoidl('2.5\times10^{-4} R_{vir}'),textoidl('10^{-3} R_{vir} '),textoidl('10^{-2} R_{vir} '),textoidl('4\times10^{-2} R_{vir}')]
res_str1 = ['50pc','200pc','2kpc','8kpc']
colors1 = [0,0,0,0]
linestyle1 = [0,0,2,1]
thickness = [3,1,3,1]
mass_str = ['12','10']
m = [12,10]

FOR imfct = 0, N_ELEMENTS(imf) -1 DO BEGIN
    set_plot,'x'
    set_plot,'ps'
    device,filename = '/astro/net/scratch1/christensen/DwarfResearch/results/sfr_grid_'+imf[imfct]+'.eps',/color,bits_per_pixel=8,/times
    !P.MULTI = [0]
    !P.MULTI = [0,2,2,0,0]
    !P.FONT = 0
    FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
        FOR ploti = 0, N_ELEMENTS(base) - 1 DO BEGIN
            IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00300' ELSE files = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00300'
            xtitle = 'Time [Gyr]'
            multiplot
            time = fltarr(N_ELEMENTS(files))+3
            mass = fltarr(N_ELEMENTS(files))+m[mct]
            ytitle = textoidl('SFR/M_B [yr^{-1}]')
            IF(mct eq 0) THEN BEGIN
                IF (ploti eq 0) THEN subcompSFR,files,time, mass,colors0, linestyle = linestyle0,ytitle = ytitle, key = legend_name0,xtitle = '',thickness=thickness,xmax=2.99,ytickname = ['0',textoidl('10^{-10}'),textoidl('2\times10^{-10}')] ELSE subcompSFR,files,time, mass,colors1,linestyle = linestyle1, xtitle = ' ', ytitle = ' ',key=legend_name1,thickness=thickness,xmax=3.0
            ENDIF ELSE BEGIN
               IF (ploti eq 0) THEN subcompSFR,files, time,mass,colors0,linestyle = linestyle0, ytitle = ytitle, xtitle = xtitle,thickness=thickness,xmax=2.99, ytickname = ['0',textoidl('10^{-10}'),textoidl('2\times10^{-10}')] ELSE subcompSFR,files,time,mass,colors1,linestyle = linestyle1, xtitle = xtitle,ytitle = '',thickness=thickness,xmax = 3.0
            ENDELSE
        ENDFOR
    ENDFOR
    xyouts,[3000,10100,8100,15100],[11000,11000,6000,6000],[textoidl('10^{12}M')+sunsymbol(),textoidl('10^{12}M')+sunsymbol(),textoidl('10^{10}M')+sunsymbol(),textoidl('10^{10}M')+sunsymbol()],charsize = 1.0,/device
     device,/close
ENDFOR
END

pro subcompSFR,files,time,m,colors,linestyle = linestyle, key = key,xtitle = xtitle, ytitle = ytitle,print=print, outplot = outplot,thickness=thickness,xmax=xmax,ytickname = ytickname
;generic program to plot sfr of a given galaxy or two
;files =['ResMassTests/1E5R/10M_nofeedback/o10M_1.00300','ResMassTests/1E4R/10M_nofeedback/o10M_1.00300','ResMassTests/1E3R/10M_nofeedback/o10M_1.00300','ResMassTests/1E2R/10M_nofeedback/o10M_1.00300']
nfiles = N_ELEMENTS(files)

binsize=5.e7
timeunit=1e9
massunit=2.325e5
time = time*timeunit
integrate = dblarr(nfiles)
gmass = dblarr(nfiles)
lines = linestyle

base = '/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/1E5R/'
;addpath,'/astro/users/stinson/idl/trace'
FOR fct=0,nfiles -1 DO BEGIN
    first = 0
    start = 0
    ymax=0   
    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
;    file=base + files[fct]
  ;  stop
    file = files[fct]
    rtipsy,file,h,g,d,s
    m[fct] = ALOG10((total(g.mass) + total(d.mass) + total(s.mass))*massunit)
    gmass[fct] = (total(g.mass) + total(s.mass))*massunit
    distance_s = SQRT(s.x^2 + s.y^2 + s.z^2)
 ;   s = s[where(distance_s le 1)]
    IF (s[0].mass NE 0)THEN BEGIN
        ind=WHERE(s.tform GT 0.0)
        tform=s[ind].tform*timeunit
        mass=s[ind].mass*massunit
        mintime=MIN(tform)
        maxtime=MAX(tform)
        nbins=FIX((maxtime-mintime)/binsize)+1
        timearray=FINDGEN(nbins)*binsize+mintime
        maxtime = 3e9
        mintime = 0
        nbins = (maxtime-mintime)/binsize
        timearray=FINDGEN(nbins)*binsize+mintime
        sfr=FLTARR(nbins)
        FOR i=0,nbins-1 DO BEGIN
            inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize),ninbin)
            IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/gmass[fct] ;massin bin/time for each bin/total mass of galaxy/fraction that is gas
        ENDFOR
        print,MAX(sfr),MAX(sfr)*gmass[fct],TOTAL(mass)/time[fct]
        ymaxtemp=MAX(sfr)
        IF (ymaxtemp GT ymax)THEN ymax=ymaxtemp
        integrate[fct]=TOTAL(mass)/time[fct]
    ENDIF ELSE integrate[fct]=0.0000
ENDFOR
ymax = 2.9e-10
;set_plot,'ps'
;device,filename='/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/results/'+outroot+'.eps',/color,bits_per_pixel=8
;stop
FOR fct=0,nfiles - 1 DO BEGIN
    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
    file = base + files[fct]
    file = files[fct]
    rtipsy,file,h,g,d,s
    gmass[fct] = (total(g.mass) + total(s.mass))*massunit
    IF (s[0].mass NE 0)THEN BEGIN
        IF first EQ 0 THEN BEGIN
            first = 1
            print,'gmass: ',gmass[fct]
            sfr,s,massu=2.325e5,time=timeunit,OVERPLOT=0,xrange=[0,xmax],binsize=binsize,gmass=gmass[fct],xtitle = xtitle, ytitle = ytitle, linestyle = lines[fct],yrange = [0,ymax],thick=thickness[fct],ytickinterval=1e-10,ytickname = ytickname;,yrange = [0,4e-10]
        ENDIF ELSE BEGIN
            first = first+1
        ENDELSE
        sfr,s,massu=2.325e5,time=timeunit,COLOR=COLORS[fct],OVERPLOT=1,binsize=binsize,gmass=gmass[fct] ,linestyle = lines[fct],thick=thickness[fct]
    ENDIF ELSE start = start+1
ENDFOR
IF (KEYWORD_SET(key)) THEN legend,key,linestyle = lines, color = colors,thick=thickness,/right,charsize=1.0,pspacing = 1.5*1.25
END
