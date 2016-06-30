pro compSFR
;This graphs the (star fromation rate)/(gass mass) over time for galaxies with
;different resolution.  One graph per galactic mass is produced.

;********IF this code doesn't seem to work, make sure that my version
;(in the same directory) of SFR.pro is compiled


;********This code is specifically to produce sfr plots to explain the
;wierd 10M sfr

files=['1E5R/10M','1E5R/10M_stable','1E5R/10M_jeans','1E5R/10M_sfr2','1E4R/10M', '1E3R/10M', '1E2R/10M', '5E1R/10M']
n_files = N_ELEMENTS(files)

integrate=fltarr(n_files)
binsize=1.e8
timeunit=1e9
massunit=2.325e5
loadct,39
set_plot,'ps'
;set_plot,'x'
colors = [50,50,50,50,100,150,200,240]
styles = [0,1,2,3,0,0,0,0]

avesfr = dblarr(n_files) ;An array that will hold the ave sfr/mass of galaxy for each different file

base = '/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/'
addpath,'/astro/users/stinson/idl/trace'

first = 0
start = 0
ymax=0   
device,filename='/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/results/10M__massRes_sfr2.eps',/color,bits_per_pixel=8
FOR fct=0, n_files - 1 DO BEGIN
    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
    if (fct eq 3)then file=base+files[fct]+'/o10M_1.00200' else file=base+files[fct]+'/o10M_1.00300'
    rtipsy,file,h,g,d,s
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
            IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/10.^10.*10. ;massin bin/time for each bin/total mass of galaxy/fraction that is gas
        ENDFOR
        ymaxtemp=MAX(sfr)
        IF (ymaxtemp GT ymax)then ymax=ymaxtemp
        if (fct eq 3) then integrate[fct]=TOTAL(mass)/2.0e9 else integrate[fct]=TOTAL(mass)/3.0e9
    ENDIF ELSE integrate[fct]=0.0000
    avesfr[fct] = integrate[fct] ;Integrated SFR/gass mass divided by time
    print,file,avesfr[fct]
ENDFOR

ymax=1.01*ymax

FOR fct=0,n_files - 1 DO BEGIN
    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
;    file=base+files[fct]+'/o10M_1.00300'
     if (fct eq 3)then file=base+files[fct]+'/o10M_1.00200' else file=base+files[fct]+'/o10M_1.00300'
rtipsy,file,h,g,d,s
    IF (s[0].mass NE 0)THEN BEGIN
        IF first EQ 0 THEN BEGIN
            first = 1
            sfr,s,massu=2.325e5,time=1e9,OVERPLOT=0,xrange=[0,3],binsize=binsize,gamss='10M',ymax=ymax,title = 'SFR/Gas Mass -- Total Mass: 1e10 -- 200pc'
        ENDIF ELSE BEGIN
            first = first+1
        ENDELSE
        sfr,s,massu=2.325e5,time=1e9,OVERPLOT=1,binsize=binsize,gmass='10M',COLOR=COLORS[first-1], linestyle=styles[first-1]
    ENDIF ELSE start = start+1
ENDFOR
legend,['100K '+trim(integrate[0]),'100K fixed '+trim(integrate[1]),'100K jeans1 '+trim(integrate[2]),'100K sfr'+trim(integrate[3]),'10K'+trim(integrate[4]),'1K '+trim(integrate[5]),'100 '+trim(integrate[6]),'50 '+trim(integrate[7])],linestyle=[0,1,2,3,0,0,0,0], color=[colors[0],colors[1],colors[2],colors[3],colors[4],colors[5],colors[6],colors[7]],/right
device,/close
END
