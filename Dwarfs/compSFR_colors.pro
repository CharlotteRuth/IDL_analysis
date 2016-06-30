pro compSFR_color
;This finds the amount of star formation in the last 400 Million years (about
;the lifetime of an A star) as this
;should relate to to the color of the galaxy 

;********IF this code doesn't seem to work, make sure that my version
;(in the same directory) of SFR.pro is compiled

r=['1E5R','1E4R', '1E3R', '1E2R', '5E1R']
;m=['9M', '10M', '11M','12M','13M']
m=['10M','12M']
sp=['8kpc', '2kpc', '1kpc','500pc','200pc','50pc']
n_mass = N_ELEMENTS(m)
n_res = N_ELEMENTS(r)
n_spec = N_ELEMENTS(sp)

;integrate=fltarr(n_res)
integrate=fltarr(n_spec)
binsize=1.e8
timeunit=1e9
massunit=2.325e5
loadct,39
;set_plot,'ps'
set_plot,'x'
;colors = [50,100,150,200,240]
colors = [25,75,100,150,200,240]

;avesfr = dblarr(n_res,n_mass) 
;A 2D array that will hold the ave sfr/mass of galaxy for each
;different spaciatial resolution and each different mass
avesfr = dblarr(n_res,n_spec) 

;base = '/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/'
base = '/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/'
addpath,'/astro/users/stinson/idl/trace'

FOR mct=0,n_mass-1 DO BEGIN
FOR sct = 0, n_spec-1 DO BEGIN
    FOR rct=0,n_res-1 DO BEGIN
        s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
 ;       cd,base+r[rct]+'/'+m[mct]
        cd,base+r[rct]+'/'+m[mct]+'/'+sp[sct]
 ;       file='o'+m[mct]+'_1.00300'
        file=sp[sct]+'.00300'
        rtipsy,file,h,g,d,s
        ind=WHERE(s.tform GT 2.6)
        IF (ind[0] NE -1)THEN BEGIN
            ind=WHERE(s.tform GT 2.6)
            tform=s.tform*timeunit
            mass=s[ind].mass*massunit
            mintime=MIN(tform)
            maxtime=MAX(tform)
            nbins=FIX((maxtime-mintime)/binsize)+1
            timearray=FINDGEN(nbins)*binsize+mintime
            sfr=FLTARR(nbins)
            FOR i=0,nbins-1 DO BEGIN
                inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize),ninbin)
                IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/10.^m[mct]*10. ;massin bin/time for each bin/total mass of galaxy/fraction that is gas
            ENDFOR
            integrate[sct]=TOTAL(mass)
        ENDIF ELSE integrate[sct]=0.0000
        avesfr[rct,sct] = integrate[sct] ;Integrated SFR/gass mass divided by time
        cd,'../..'
    ENDFOR
ENDFOR

set_plot,'ps'
;set_plot,'x'
;device,filename='/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/results/massRes_recentSFR.eps',/color,bits_per_pixel=8
device,filename='/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/results/specRes_recentSFR'+trim(m[mct])+'.eps',/color,bits_per_pixel=8
!y.style=0
if (mct eq 0) then plot,[1e5,1e4,1e3,1e2,50],[1e5,1e5,1e5,1e5,1e5],xtitle = 'Number of Particles',ytitle='M'+sunsymbol()+' Formed in past 400 Myrs',/xlog,/ylog,xrange=[10,1e6],yrange=[1e6,1e8],title='Effect of Resolution on Recent SF -- 10^10M'+sunsymbol()$
else plot,[1e5,1e4,1e3,1e2,50],[1e7,1e7,1e7,1e7,1e7],xtitle = 'Number of Particles',ytitle='M'+sunsymbol()+' Formed in past 400 Myrs',/xlog,/ylog,xrange=[10,1e6],yrange=[1e8,1e10],title='Effect of Resolution on Recent SF -- 10^12M'+sunsymbol()

oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,0],linestyle = 0,psym=-1,color=colors[0]
oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,1],linestyle = 0,psym=-2,color=colors[1]
oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,2],linestyle = 0,psym=-4,color=colors[2]
oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,3],linestyle = 0,psym=-5,color=colors[3]
oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,4],linestyle = 0,psym=-6,color=colors[4]
oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,5],linestyle = 0,psym=-7,color=colors[5]

;legend,['10^'+trim(m[4])+sunsymbol(),'10^'+trim(m[3])+sunsymbol(),'10^'+trim(m[2])+sunsymbol(),'10^'+trim(m[1])+sunsymbol(),'10^'+trim(m[0])+sunsymbol()],linestyle=[0,0,0,0,0],psym=[-6,-5,-4,-2,-1]
legend,['8kpc','2kpc','1kpc','500pc','200pc','50pc'],linestyle=[0,0,0,0,0,0],psym=[-1,-2,-4,-5,-6,-7],color = [colors[0],colors[1],colors[2],colors[3],colors[4],colors[5]]
device,/close
ENDFOR
cd,'/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/procedures'
END
