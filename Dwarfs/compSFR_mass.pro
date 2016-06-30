pro compSFR_mass
;This graphs the (star fromation rate)/(gass mass) over time for galaxies with
;different resolution.  One graph per galactic mass is produced.

;********IF this code doesn't seem to work, make sure that my version
;(in the same directory) of SFR.pro is compiled

r=['1E6R','1E5R','1E4R', '1E3R', '1E2R', '5E1R']
m=['9M','10M','11M','12M','13M']
imf=['k']
gmass = 10.^m/10. ;Somehow this works . . .
n_mass = N_ELEMENTS(m)
n_res = N_ELEMENTS(r)
n_imf = N_ELEMENTS(imf)

binsize=1.e8
timeunit=1e9
massunit=2.325e5
loadct,39
;set_plot,'ps'
!P.thick = 1.5
!X.Charsize = 1.35
!Y.Charsize = 1.35
!P.CHARSIZE = 1.5
!X.STYLE = 1
!P.FONT = 0
set_plot,'x'
colors = [20,20,70,70,100,100,150,150,240,240]

avesfr = dblarr(n_res,n_mass,n_imf) ;A 3D array that will hold the ave sfr/mass of galaxy for each different spaciatial resolution, imf, and each different mass

base = '/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/'
addpath,'/astro/users/stinson/idl/trace'


for mct=0,n_mass-1 do begin
    integrate=fltarr(n_res*n_imf)    
    first = 0
    start = 0
    ymax=0   
;    device,filename='/astro/net/scratch1/christensen/DwarfResearch/results/'+m[mct]+'_massRes_sfr.eps',/color,bits_per_pixel=8
                                ;  device,filename='/astro/net/scratch1/christensen/DwarfResearch/results/'+m[mct]+'_massRes_sfr_'+imf[0]+'.eps',/color,bits_per_pixel=8
    IF (mct ne 1 AND mct ne 3) then rct_min = 1 else rct_min = 0
    FOR rct=rct_min,n_res - 1 DO BEGIN
        FOR imfct = 0,n_imf-1 DO BEGIN
            s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
            file=base+r[rct]+'/' + m[mct] + '_'+imf[imfct]+'/'+'o'+m[mct]+'_1.00300'
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
                    IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/gmass[mct] ;massin bin/time for each bin/total mass of galaxy/fraction that is gas
                ENDFOR
                ymaxtemp=MAX(sfr)
                IF (ymaxtemp GT ymax)then ymax=ymaxtemp
                integrate[rct*n_imf+imfct]=TOTAL(mass)/3.0e9
               ; if(mct gt 2 AND imfct eq 1) THEN stop
            ENDIF ELSE integrate[rct]=0.0000
            avesfr[rct,mct,imfct] = integrate[rct*n_imf+imfct];Integrated SFR/gass mass divided by time
        ENDFOR
    ENDFOR

    ymax=1.01*ymax

    FOR rct=rct_min,n_res-1 DO BEGIN
        FOR imfct = 0,n_imf-1 DO BEGIN
            s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
            file=base+r[rct]+'/' + m[mct] + '_'+imf[imfct]+'/'+'o'+m[mct]+'_1.00300'
            rtipsy,file,h,g,d,s
            IF (s[0].mass NE 0)THEN BEGIN
                IF first EQ 0 THEN BEGIN
                    first = 1
                    sfr,s,massu=2.325e5,time=1e9,OVERPLOT=0,xrange=[0,3],binsize=binsize,gamss=gmass[mct],ymax=ymax,title = 'SFR/Gas Mass -- Total Mass: 10^'+trim(m[mct])+sunsymbol()+' -- 200pc',linestyle = imfct*2,xtitle = 'Time (Gyr)', ytitle='Average SFR/(Galactic Initial Gas Mass) (yr^-1)'
                ENDIF ELSE BEGIN
                    first = first+1
                ENDELSE
                sfr,s,massu=2.325e5,time=1e9,OVERPLOT=1,binsize=binsize,gmass=gmass[mct],COLOR=COLORS[(first-1)*2],linestyle = imfct*2
            ENDIF ELSE start = start+1
            print,file,' ',mct,' ',rct,' ',imfct
        ENDFOR
    ENDFOR

;    legend,['1M '+trim(integrate[0]),'100K '+trim(integrate[1]),'10K '+trim(integrate[2]),'1K '+trim(integrate[3])],linestyle=[0,0,0,0],color=[colors[0],colors[1],colors[2],colors[3]],/right
    ;legend,['100K '+trim(integrate[0]),'10K '+trim(integrate[1]),'1K '+trim(integrate[2]),'100 '+trim(integrate[3]),'50 '+trim(integrate[4])],  linestyle=[0,0,0,0,0],    color=[colors[0],colors[1],colors[2],colors[3],colors[4]],/right
    legend,['100K '+trim(integrate[0]),'10K '+trim(integrate[1]),'1K '+trim(integrate[2]),'100 '+trim(integrate[3]),'50 '+trim(integrate[4])],  linestyle=[0,0,0,0,0],    color=[colors[0],colors[2],colors[4],colors[6],colors[8]],/right
;    legend,['100K, K'+trim(integrate[0]),'100K, MS '+trim(integrate[1]),'10K, K '+trim(integrate[2]),'10K, MS '+trim(integrate[3]),'1K, K '+trim(integrate[4]),'1K, MS '+trim(integrate[5]),'100, K '+trim(integrate[6]),'100, MS '+trim(integrate[7]),'50, K '+trim(integrate[8]),'50, MS '+trim(integrate[9])],  linestyle=[0,2,0,2,0,2,0,2,0,2], color=[colors[0],colors[1],colors[2],colors[3],colors[4],colors[5],colors[6],colors[7],colors[8],colors[9]],/right
;    stop
;    device,/close
endfor

;set_plot,'ps'
;device,filename='/astro/net/scratch1/christensen/DwarfResearch/results/massRes_sfr.eps',/color,bits_per_pixel=8
;!y.style=0
;plot,[1e5,1e4,1e3,1e2,50],avesfr[*,0,0],xtitle = 'Number of Particles',ytitle='Average SFR/(Galactic Initial Gass Mass) (yr^-1)',xrange=[10,1e6],/xlog,yrange=[1e-4,1e4],/ylog,linestyle = 0,psym=-1,title='Effect of Resolution and Mass on SFR'
;oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,0,0],linestyle = 0,psym=-1,color = colors[8]
;oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,0,1],linestyle = 1,psym=-1,color = colors[8]
;oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,1,0],linestyle = 0,psym=-2,color = colors[6]
;oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,1,1],linestyle = 1,psym=-2,color = colors[6]
;oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,2,0],linestyle = 0,psym=-4,color = colors[4]
;oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,2,1],linestyle = 1,psym=-4,color = colors[4]
;oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,3,0],linestyle = 0,psym=-5,color = colors[2]
;oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,3,1],linestyle = 1,psym=-5,color = colors[2]
;oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,4,0],linestyle = 0,psym=-6,color = colors[0]
;oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,4,1],linestyle = 1,psym=-6,color = colors[0]
;legend,['K,  10^'+trim(m[4])+sunsymbol(),'MS, 10^'+trim(m[4])+sunsymbol(),'K,  10^'+trim(m[3])+sunsymbol(),'MS, 10^'+trim(m[3])+sunsymbol(),'K,  10^'+trim(m[2])+sunsymbol(),'MS, 10^'+trim(m[2])+sunsymbol(),'K,  10^'+trim(m[1])+sunsymbol(),'MS, 10^'+trim(m[1])+sunsymbol(),'K,  10^'+trim(m[0])+sunsymbol(),'MS, 10^'+trim(m[0])+sunsymbol()],linestyle=[0,1,0,1,0,1,0,1,0,1],psym=[-6,-6,-5,-5,-4,-4,-2,-2,-1,-1],color=[colors[0],colors[0],colors[2],colors[2],colors[4],colors[4],colors[6],colors[6],colors[8],colors[8]]
;device,/close


set_plot,'ps'
device,filename='/astro/net/scratch1/christensen/DwarfResearch/results/massRes_norm1e4_sfr_'+imf[0]+'.eps',/color,bits_per_pixel=8,/TIMES
!x.style = 1
!y.style=1
plot,[1e5,1e4,1e3,1e2,50],avesfr[1:n_res - 1,0,0]/avesfr[3,0,0],xtitle = 'Number of DM Particles',ytitle=textoidl('<SFR>/<SFR>_{N = 1000}'),xrange=[50,2e7],/xlog,linestyle = 1,psym=-7,yrange = [0,1.7];,ytickname = [textoidl('10^{-4}'),textoidl('10^{-3}'),textoidl('10^{-2}'),textoidl('10^{-1}'),'10',textoidl('10^{1}'),textoidl('10^{2}')];format = '(g10.1)';,title='Effect of Resolution and Mass on SFR'
oplot,[1e5,1e4,1e3,1e2,50],avesfr[1:n_res - 1,4,0]/avesfr[3,4,0],linestyle = 0,psym=-2,thick=3
oplot,[1e6,1e5,1e4,1e3,1e2,50],avesfr[*,3,0]/avesfr[3,3,0],linestyle = 0,psym=-4,thick=1
oplot,[1e5,1e4,1e3,1e2,50],avesfr[1:n_res -1,2,0]/avesfr[3,2,0],linestyle = 2,psym=-5,thick=3
oplot,[1e6,1e5,1e4,1e3,1e2,50],avesfr[*,1,0]/avesfr[3,1,0],linestyle = 2,psym=-6,thick=1
oplot,[1e5,1e4,1e3,1e2,50],avesfr[1:n_res - 1,0,0]/avesfr[3,0,0],linestyle = 3,psym=-7,thick=3
legend,[textoidl('10^{13}')+'M'+sunsymbol(),textoidl('10^{12}')+'M'+sunsymbol(),textoidl('10^{11}')+'M'+sunsymbol(),textoidl('10^{10}')+'M'+sunsymbol(),textoidl('10^{9}')+'M'+sunsymbol()],linestyle=[0,0,2,2,3],psym=[-2,-4,-5,-6,-7],thick=[3,1,3,1,3],charsize = 1.5, number = 0.5,/bottom,/right,pspacing = 1.5*1.5;,color=[colors[0],colors[2],colors[4],colors[6],colors[8]]
device,/close

set_plot,'ps'
device,filename='/astro/net/scratch1/christensen/DwarfResearch/results/massRes_sfr_'+imf[0]+'.eps',/color,bits_per_pixel=8,/TIMES
!x.style = 1
!y.style=1
plot,[1e5,1e4,1e3,1e2,50],avesfr[1:n_res - 1,0,0],xtitle = 'Number of DM Particles',ytitle=textoidl('<SFR> [M')+sunsymbol()+textoidl(' yr^{-1}]'),xrange=[50,2e7],/xlog,yrange=[1e-4,2e2],/ylog,linestyle = 1,psym=-7,ytickname = [textoidl('10^{-4}'),textoidl('10^{-3}'),textoidl('10^{-2}'),textoidl('10^{-1}'),'10',textoidl('10^{1}'),textoidl('10^{2}')];format = '(g10.1)';,title='Effect of Resolution and Mass on SFR'
oplot,[1e5,1e4,1e3,1e2,50],avesfr[1:n_res - 1,4,0],linestyle = 0,psym=-2,thick=3
oplot,[1e6,1e5,1e4,1e3,1e2,50],avesfr[*,3,0],linestyle = 0,psym=-4,thick=1
oplot,[1e5,1e4,1e3,1e2,50],avesfr[1:n_res -1,2,0],linestyle = 2,psym=-5,thick=3
oplot,[1e6,1e5,1e4,1e3,1e2,50],avesfr[*,1,0],linestyle = 2,psym=-6,thick=1
oplot,[1e5,1e4,1e3,1e2,50],avesfr[1:n_res - 1,0,0],linestyle = 3,psym=-7,thick=3
legend,[textoidl('10^{13}')+'M'+sunsymbol(),textoidl('10^{12}')+'M'+sunsymbol(),textoidl('10^{11}')+'M'+sunsymbol(),textoidl('10^{10}')+'M'+sunsymbol(),textoidl('10^{9}')+'M'+sunsymbol()],linestyle=[0,0,2,2,3],psym=[-2,-4,-5,-6,-7],thick=[3,1,3,1,3],charsize = 1.5, number = 0.5,/bottom,/right,pspacing = 1.5*1.5;,color=[colors[0],colors[2],colors[4],colors[6],colors[8]]
device,/close

set_plot,'ps'
device,filename='/astro/net/scratch1/christensen/DwarfResearch/results/massRes_normsfr_'+imf[0]+'.eps',/color,bits_per_pixel=8,/TIMES
device,/times
!y.style=1
plot,[1e5,1e4,1e3,1e2,50],avesfr[*,0,0]/gmass[0],xtitle = 'Number of DM Particles',ytitle=textoidl('<SFR>/M_{g} [yr^{-1}]'),xrange=[1,2e5],/xlog,/ylog,linestyle = 1,psym=-7,yrange=[7e-13,3e-10];format = '(g10.1)';,title='Effect of Resolution and Mass on SFR'
oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,4,0]/gmass[4],linestyle = 0,psym=-2,color = 0
oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,3,0]/gmass[3],linestyle = 2,psym=-4,color = 0
oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,2,0]/gmass[2],linestyle = 0,psym=-5,color = 100
oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,1,0]/gmass[1],linestyle = 2,psym=-6,color = 100
oplot,[1e5,1e4,1e3,1e2,50],avesfr[*,0,0]/gmass[0],linestyle = 1,psym=-7,color =0
legend,[textoidl('10^{13}')+'M'+sunsymbol(),textoidl('10^{12}')+'M'+sunsymbol(),textoidl('10^{11}')+'M'+sunsymbol(),textoidl('10^{10}')+'M'+sunsymbol(),textoidl('10^{9}')+'M'+sunsymbol()],linestyle=[0,2,0,2,1],psym=[-2,-4,-5,-6,-7],color = [0,0,100,100,0],charsize = 1.5, number = 2,pspacing = 1.5*1.5;,color=[colors[0],colors[2],colors[4],colors[6],colors[8]]
device,/close
;stop

END
