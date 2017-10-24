pro compSFR_mod,files,time,colors,key,outroot
;generic program to plot sfr of a given galaxy or two
;This graphs the (star fromation rate)/(gass mass) vs time in
;different ways and is taken from compSFR_spec.pro

files =['/astro/net/scratch2/fabio/Xruns/h1201-h1755-X2X2g1bwK/h1201-h1755-X2X2g1bwK.00512.2686.00.8.rot','/astro/net/scratch2/fabio/Xruns/h1201-h1755-X5X5g1bwK/h1201-h1755-X5X5g1bwK.00512.5330.00.1.rot']
time = [0.333091,0.333091]
colors = [50,240]
key = ['X2X2 8','X5X5 1']
outroot = 'sfr_1_8'

nfiles = N_ELEMENTS(files)

msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 1.69875e16 ; Solar mass in system units
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 50000.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70

binsize=1.e8
timeunit=4.0452811e10
massunit=dMsolUnit;2.325e5
loadct,39
set_plot,'ps'
set_plot,'x'
colors = [50,100,150,200,240]
time = time*timeunit
m = Fltarr(N_ELEMENTS(files))

integrate = dblarr(nfiles)

base = ''
addpath,'/astro/users/stinson/idl/trace'
;device,filename=''+outroot+'.eps',/color,bits_per_pixel=8
FOR fct=0,nfiles -1 DO BEGIN
    first = 0
    start = 0
    ymax=0   
    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
    file=base + files[fct]
    rtipsy,file,h,g,d,s
    totalmass = TOTAL(g.mass + d.mass + s.mass) * dMsolUnit
    m[fct] = totalmass
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
            IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/m[fct] ;massin bin/time for each bin/total mass of galaxy/fraction that is gas
        ENDFOR
        ymaxtemp=MAX(sfr)
        IF (ymaxtemp GT ymax)THEN ymax=ymaxtemp
        integrate[fct]=TOTAL(mass)/time[fct]
    ENDIF ELSE integrate[fct]=0.0000
ENDFOR
stop
ymax=1.01*ymax

FOR fct=0,nfiles - 1 DO BEGIN
    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
    file = base + files[fct]
    rtipsy,file,h,g,d,s
    IF (s[0].mass NE 0)THEN BEGIN
        print,MAX(time/timeunit)
        IF first EQ 0 THEN BEGIN
            first = 1
            sfr,s,massu=massunit,time=timeunit,OVERPLOT = 0,binsize=binsize,gmass=ALOG10(m[fct]),title = 'SFR/Gas Mass'
        ENDIF ELSE BEGIN
            first = first+1
        ENDELSE
        sfr,s,massu=massunit,time=timeunit,COLOR=COLORS[fct],OVERPLOT=1,binsize=binsize,gmass=ALOG10(m[fct]) 
    ENDIF ELSE start = start+1
ENDFOR
IF (nfiles eq 1) THEN legend,[key[0]+' '+trim(integrate[0])],linestyle=[0],color=[colors[0]],/right
IF (nfiles eq 2) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1])],linestyle=[0,0],color=[colors[0],colors[1]],/right
IF (nfiles eq 3) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2])],linestyle=[0,0,0],color=[colors[0],colors[1],colors[2]],/right
IF (nfiles eq 4) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3])],linestyle=[0,0,0,0],color=[colors[0],colors[1],colors[2],colors[3]],/right
IF (nfiles eq 5) THEN legend,[key[0]+' '+trim(integrate[0]),key[1]+' '+trim(integrate[1]),key[2]+' '+trim(integrate[2]),key[3]+' '+trim(integrate[3]),key[4]+' '+trim(integrate[4])],linestyle=[0,0,0,0,0],color=[colors[0],colors[1],colors[2],colors[3],colors[4]],/right
;device,/close
stop
END
