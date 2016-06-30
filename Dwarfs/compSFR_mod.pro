PRO compSFR_master
type = ['s','d','g']
type_str = ['Stellar','Dark','Gas']

loadct,39
!P.CHARSIZE = 1.25
base = ['/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/','/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/']
imf = ['k'];,'ms']
legend_name0 = ['1000K','100K','10K','1K']
res_str0 = ['1E6R','1E5R','1E4R','1E3R']
colors0 = [0,0,0,0]
linestyle0 = [0,0,2,2]
legend_name1 = [textoidl('R_{vir}/4000'),textoidl('R_{vir}/1000'),textoidl('R_{vir}/100'),textoidl('R_{vir}/25')]
res_str1 = ['50pc','200pc','2kpc','8kpc']
colors1 = [0,0,0,0]
linestyle1 = [0,0,2,2]
thickness = [3,1,3,1]
mass_str = ['12','10']
m = [12,10]

FOR imfct = 0, N_ELEMENTS(imf) -1 DO BEGIN
;    set_plot,'x'
    !P.MULTI = [0]
    multiplot,/reset
    set_plot,'x'
    set_plot,'ps'
    device,filename = '/astro/net/scratch1/christensen/DwarfResearch/results/sfr_grid_'+imf[imfct]+'_inner1kpc.eps',/color,bits_per_pixel=8
    !P.MULTI = [0,2,2,0,0]
    FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
        FOR ploti = 0, N_ELEMENTS(base) - 1 DO BEGIN
            IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00300' ELSE files = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00300'
            xtitle = 'Time [Gyr]'
            multiplot
            time = fltarr(N_ELEMENTS(files))+3
            mass = fltarr(N_ELEMENTS(files))+m[mct]
            ytitle = 'SFR/Initial Gas Mass'
            IF(mct eq 0) THEN BEGIN
                IF (ploti eq 0) THEN compSFR_mod,files,time, mass,colors0, linestyle = linestyle0,ytitle = ytitle, key = legend_name0,xtitle = '',thickness=thickness ELSE compSFR_mod,files,time, mass,colors1,linestyle = linestyle1, xtitle = ' ', ytitle = ' ',key=legend_name1,thickness=thickness
            ENDIF ELSE BEGIN
               IF (ploti eq 0) THEN compSFR_mod,files, time,mass,colors0,linestyle = linestyle0, ytitle = ytitle, xtitle = xtitle,thickness=thickness ELSE compSFR_mod,files,time,mass,colors1,linestyle = linestyle1, xtitle = xtitle,ytitle = '',thickness=thickness
            ENDELSE
        ENDFOR
    ENDFOR
    xyouts,[3000,10000,3000,10000],[11000,11000,6000,6000],[textoidl('10^{12}'),textoidl('10^{12}'),textoidl('10^{10}'),textoidl('10^{10}')],/device
     device,/close
ENDFOR
END

pro compSFR_mod,files,time,m,colors,linestyle = linestyle, key = key,xtitle = xtitle, ytitle = ytitle,print=print, outplot = outplot,thickness=thickness
;generic program to plot sfr of a given galaxy or two
;files =['ResMassTests/1E5R/10M_nofeedback/o10M_1.00300','ResMassTests/1E4R/10M_nofeedback/o10M_1.00300','ResMassTests/1E3R/10M_nofeedback/o10M_1.00300','ResMassTests/1E2R/10M_nofeedback/o10M_1.00300']
;time = [3.0,3.0,3.0,3.0]
;colors = [50,100,175,240]
;m = [10,10,10,10]
;key = ['1E5R','1E4R','1E3R','1E2R']
;outroot = 'nofeedback'

;files =['ResMassTests/1E5R/10M_original/o10M_1.00300','ResMassTests/1E4R/10M/o10M_1.00300','ResMassTests/1E5R/lowResGas/o10M.00300','ResMassTests/1E5R/lowResDark/o10M.00300']
;time = [3.0,3.0,3.0,3.0]
;colors = [30,70,130,240]
;m = [10,10,10,10]
;key = ['1E5 G, 1E5 D','1E4 G, 1E4 D','1E4 G, 1E5 D','1E4 D, 1E5 G']
;outroot = 'd_g_res'
;files =['ResMassTests/1E5R/10M_kroupa_robert/o10M_1.00300','ResMassTests/1E5R/10M_ms_robert/o10M_1.00300','ResMassTests/1E4R/10M_kroupa_robert/o10M.00300','ResMassTests/1E4R/10M_ms_robert/o10M.00300']
;time = [3.0,3.0,3.0,3.0]
;colors = [50,50,240,240]
;m = [10,10,10,10]
;key = ['1E5R Kroupa','1E5R Miller & Scalo','1E4R Kroupa','1E4R Miller & Scalo']
;outroot = 'diffIMF'
;lines = [0,2,0,2]

;files =['13M/o13M_1.00300','13M_ms/o13M_1.00300','12M/o12M_1.00300','12M_ms/o12M_1.00300','11M/o11M_1.00300','11M_ms/o11M_1.00300','10M_kroupa_robert/o10M_1.00300','10M_ms_robert/o10M_1.00300','9M/o9M_1.00300','9M_ms/o9M_1.00300']
;time = [3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]
;colors = [20,20,80,80,140,140,200,200,240,240]
;m = [13,13,12,12,11,11,10,10,9,9]
;key = ['13 Kroupa','13 Miller & Scalo','12 Kroupa','12 Miller & Scalo','11 Kroupa','11 Miller & Scalo','10 Kroupa','10 Miller & Scalo','9 Kroupa','9 Miller & Scalo']
;outroot = 'diffIMF_mass'
;lines = [0,2,0,2,0,2,0,2,0,2]

;files =['/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.2304g2bwK/h986.cosmo50cmb.2304g2bwK.00512/h986.cosmo50cmb.2304g2bwK.00512.1.std','/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.00512/h603.cosmo50cmb.2304g2bwK.00512.1.std']
;time = [15,15]
;colors = [50,240]
;m = [10,10,10,10]
;key = ['h986','h603']
;outroot = 'twins'
;massunit = '1.84793e16'
;timeunit = 3.8785614e+10
;linestyle = [0,0]

;This graphs the (star fromation rate)/(gass mass) vs time in
;different ways and is taken from compSFR_spec.pro

nfiles = N_ELEMENTS(files)

binsize=5.e7
timeunit=1e9
massunit=2.325e5
;loadct,39
;set_plot,'ps'
;set_plot,'x'
;colors = [50,100,150,200,240]
time = time*timeunit
integrate = dblarr(nfiles)
gmass = dblarr(nfiles)
;lines = fltarr(N_ELEMENTS(files))
lines = linestyle

base = '/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/1E5R/'
addpath,'/astro/users/stinson/idl/trace'
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
    s = s[where(distance_s le 1)]
    IF (s[0].mass NE 0)THEN BEGIN
        ind=WHERE(s.tform GT 0.0)
        tform=s[ind].tform*timeunit
        mass=s[ind].mass*massunit
        mintime=MIN(tform)
        maxtime=MAX(tform)
        nbins=FIX((maxtime-mintime)/binsize)+1
        timearray=FINDGEN(nbins)*binsize+mintime
        sfr=FLTARR(nbins)
        FOR i=0,nbins-1 DO BEGIN
            inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize),ninbin)
            IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/10.^m[fct]*10 ;massin bin/time for each bin/total mass of galaxy/fraction that is gas
        ENDFOR
        ymaxtemp=MAX(sfr)
;        print,ymaxtemp
        IF (ymaxtemp GT ymax)THEN ymax=ymaxtemp
        integrate[fct]=TOTAL(mass)/time[fct]
    ENDIF ELSE integrate[fct]=0.0000
ENDFOR
ymax=1.01*ymax
ymax = 4.5e-10
;set_plot,'ps'
;device,filename='/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/results/'+outroot+'.eps',/color,bits_per_pixel=8
;stop
FOR fct=0,nfiles - 1 DO BEGIN
    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
    file = base + files[fct]
    file = files[fct]
   rtipsy,file,h,g,d,s
    distance_s = SQRT(s.x^2 + s.y^2 + s.z^2)
    s = s[where(distance_s le 1)]
gmass[fct] = (total(g.mass) + total(s.mass))*massunit
    IF (s[0].mass NE 0)THEN BEGIN
        print,MAX(time/timeunit)
        IF first EQ 0 THEN BEGIN
            first = 1
            print,'gmass: ',m[fct]
            sfr,s,massu=2.325e5,time=timeunit,OVERPLOT=0,xrange=[0,MAX(time/timeunit)],binsize=binsize,gmass=gmass[fct],xtitle = xtitle, ytitle = ytitle, linestyle = lines[fct],yrange = [0,ymax],thick=thickness[fct];,yrange = [0,4e-10]
        ENDIF ELSE BEGIN
            first = first+1
        ENDELSE
        sfr,s,massu=2.325e5,time=timeunit,COLOR=COLORS[fct],OVERPLOT=1,binsize=binsize,gmass=gmass[fct] ,linestyle = lines[fct],thick=thickness[fct]
    ENDIF ELSE start = start+1
ENDFOR
IF (KEYWORD_SET(key)) THEN legend,key,linestyle = lines, color = colors,/right
;    IF (nfiles eq 1) THEN legend,[key[0]+' '+trim(integrate[0])],linestyle=[0],color=[colors[0]],/right
;    IF (nfiles eq 2) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1])],linestyle=[0,0],color=[colors[0],colors[1]],/right
;    IF (nfiles eq 3) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2])],linestyle=[0,0,0],color=[colors[0],colors[1],colors[2]],/right
;    IF (nfiles eq 4) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3])],linestyle=[lines[0],lines[1],lines[2],lines[3]],color=[colors[0],colors[1],colors[2],colors[3]],/right
;    IF (nfiles eq 5) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4])],linestyle=[0,0,0,0,0],color=[colors[0],colors[1],colors[2],colors[3],colors[4]],/right
;    IF (nfiles eq 6) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4]),key[5]+' '+trim(integrate[5])],linestyle=[0,0,0,0,0,0],color=[colors[0],colors[1],colors[2],colors[3],colors[4],colors[5]],/right
;    IF (nfiles eq 7) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4]),key[5]+' '+trim(integrate[5]),key[6]+' '+trim(integrate[6])],linestyle=[0,0,0,0,0,0,0],color=[colors[0],colors[1],colors[2],colors[3],colors[4],colors[5],colors[6]],/left,/top
;    IF (nfiles eq 8) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4]),key[5]+' '+trim(integrate[5]),key[6]+' '+trim(integrate[6]),key[7]+' '+trim(integrate[7])],linestyle=[0,0,2,3,2,0,0,0],color=[colors[0],colors[1],colors[2],colors[3],colors[4],colors[5],colors[6],colors[7]],/left,/top
;    IF (nfiles eq 10) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4]),key[5]+' '+trim(integrate[5]),key[6]+' '+trim(integrate[6]),key[7]+' '+trim(integrate[7]),key[8]+' '+trim(integrate[8]),key[9]+' '+trim(integrate[9])],linestyle=[lines[0],lines[1],lines[2],lines[3],lines[4],lines[5],lines[6],lines[7],lines[8],lines[9]],color=[colors[0],colors[1],colors[2],colors[3],colors[4],colors[5],colors[6],colors[7],colors[8],colors[9]],/right
;ENDIF
;device,/close
;set_plot,'x'
;first = 0
;FOR fct=0,nfiles - 1 DO BEGIN
 ;   s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
  ;  file = base + files[fct]
   ; rtipsy,file,h,g,d,s
    ;IF (s[0].mass NE 0)THEN BEGIN
     ;   print,MAX(time/timeunit)
      ;  IF first EQ 0 THEN BEGIN
       ;     first = 1
        ;    sfr,s,massu=2.325e5,time=timeunit,OVERPLOT=0,xrange=[0,MAX(time/timeunit)],binsize=binsize,gmass=m[fct],title = 'SFR/Gas Mass',yrange = [0,ymax],linestyle = lines[fct]
;        ENDIF ELSE BEGIN
 ;           first = first+1
  ;      ENDELSE
   ;     sfr,s,massu=2.325e5,time=timeunit,COLOR=COLORS[fct],OVERPLOT=1,binsize=binsize,gmass=m[fct] ,linestyle = lines[fct]
   ; ENDIF ELSE start = start+1
;ENDFOR
;IF (nfiles eq 1) THEN legend,[key[0]+' '+trim(integrate[0])],linestyle=[0],color=[colors[0]],/right
;IF (nfiles eq 2) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1])],linestyle=[0,0],color=[colors[0],colors[1]],/right
;IF (nfiles eq 3) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2])],linestyle=[0,0,0],color=[colors[0],colors[1],colors[2]],/right
;IF (nfiles eq 4) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3])],linestyle=[lines[0],lines[1],lines[2],lines[3]],color=[colors[0],colors[1],colors[2],colors[3]],/right
;IF (nfiles eq 5) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4])],linestyle=[0,0,0,0,0],color=[colors[0],colors[1],colors[2],colors[3],colors[4]],/right
;IF (nfiles eq 6) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4]),key[5]+' '+trim(integrate[5])],linestyle=[0,0,0,0,0,0],color=[colors[0],colors[1],colors[2],colors[3],colors[4],colors[5]],/right
;IF (nfiles eq 7) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4]),key[5]+' '+trim(integrate[5]),key[6]+' '+trim(integrate[6])],linestyle=[0,0,0,0,0,0,0],color=[colors[0],colors[1],colors[2],colors[3],colors[4],colors[5],colors[6]],/left,/top
;IF (nfiles eq 8) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4]),key[5]+' '+trim(integrate[5]),key[6]+' '+trim(integrate[6]),key[7]+' '+trim(integrate[7])],linestyle=[0,0,0,0,0,0,0,0],color=[colors[0],colors[1],colors[2],colors[3],colors[4],colors[5],colors[6],colors[7]],/left,/top
;IF (nfiles eq 10) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4]),key[5]+' '+trim(integrate[5]),key[6]+' '+trim(integrate[6]),key[7]+' '+trim(integrate[7]),key[8]+' '+trim(integrate[8]),key[9]+' '+trim(integrate[9])],linestyle=[lines[0],lines[1],lines[2],lines[3],lines[4],lines[5],lines[6],lines[7],lines[8],lines[9]],color=[colors[0],colors[1],colors[2],colors[3],colors[4],colors[5],colors[6],colors[7],colors[8],colors[9]],/right


END
