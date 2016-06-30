pro aveSFR_spec
;This graphs the (star fromation rate)/(gass mass) over time for galaxies with
;different spatial resolution.  One graph per galactic mass is produced.
;********IF this code doesn't seem to work, make sure that my version
;(in the same directory) of SFR.pro is compiled

p=['1E5R','1E4R','1E3R','1E2R'];,'5E1R']
r=['8kpc','2kpc', '1kpc', '500pc','200pc','50pc']
resvalues = [2e-2,1e-2,5e-3,2.5e-3,1e-3,2.5e-4]
m=['10M','12M']
imf=['k']
n_mass = N_ELEMENTS(m)
n_res = N_ELEMENTS(r)
n_par = N_ELEMENTS(p)
n_imf = N_ELEMENTS(imf)

integrate=fltarr(n_res * n_imf)
binsize=1.e8
timeunit=1e9
massunit=2.325e5
loadct,39
;set_plot,'ps'
!P.thick = 1.5
!X.Charsize = 1.25
!Y.Charsize = 1.25
!P.CHARSIZE = 1.5
set_plot,'x'
  
avesfr = dblarr(n_res,n_par,n_imf) ;A 2D array that will hold the ave sfr/mass of galaxy for each different spaciatial resolution and each different mass
base = '/astro/net/scratch1/christensen/DwarfResearch/'

FOR mct = 0,n_mass-1 DO BEGIN
    FOR pct=0,n_par -1 DO BEGIN
        FOR rct=0,n_res-1 DO BEGIN
            FOR imfct = 0, n_imf-1 DO BEGIN
                s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
                file=base+'ResSpecTests/'+p[pct]+'/'+m[mct]+'/'+r[rct]+'_'+imf[imfct]+'/'+r[rct]+'.00300'
                rtipsy,file,h,g,d,s
                IF (s[0].mass NE 0)THEN BEGIN
                    ind=WHERE(s.tform GT 0.0)
                    tform=s[ind].tform*timeunit
                    mass=s[ind].mass*massunit
                    integrate[rct]=TOTAL(mass)/3.0e9
                ENDIF ELSE integrate[rct]=0.0000
                avesfr[rct,pct,imfct] = integrate[rct]
                print,file,integrate[rct]

            ENDFOR
        ENDFOR
    ENDFOR

    set_plot,'ps'
    device,filename=base+'results/massSpec_sfr_'+m[mct]+'_'+imf[0]+'.eps',/color,bits_per_pixel=8
        !Y.style = 1
        if (mct eq 0) then plot,resvalues,[10,10,10,10,10],xtitle = 'Softening Length ' + textoidl('[R_{vir}]'),ytitle=textoidl('<SFR> [M')+sunsymbol()+textoidl(' yr^{-1}]'),/xlog,/ylog,xrange = [1e-1,1e-4],yrange = [4e-4,0.25],XTICKNAME = [textoidl('10^{-1}'),textoidl('10^{-2}'),textoidl('10^{-3}'),textoidl('10^{-4}'),'1'],YTICKNAME=[textoidl('10^{-3}'),textoidl('10^{-2}'),textoidl('10^{-1}')]$
         else plot,resvalues,[1,1,1,1,1],xtitle = 'Softening Length ' + textoidl('[R_{vir}]'),ytitle=textoidl('<SFR> [M_')+sunsymbol()+textoidl(' yr^{-1}]'),/xlog,/ylog,xrange = [0.1,0.0001],yrange = [6,25],XTICKNAME = [textoidl('10^{-1}'),textoidl('10^{-2}'),textoidl('10^{-3}'),textoidl('10^{-4}'),'1']
        oplot,resvalues,avesfr[*,0,0],linestyle = 0,psym=-2,thick=3
        oplot,resvalues,avesfr[*,1,0],linestyle = 0,psym=-4,thick = 1
        oplot,resvalues,avesfr[*,2,0],linestyle = 2,psym=-5,thick = 3
        oplot,resvalues,avesfr[*,3,0],linestyle = 2,psym=-6,thick = 1
        if (mct eq 0) then begin
            legend,[textoidl('10^5'),textoidl('10^4'),textoidl('10^3'),textoidl('10^2')],linestyle=[0,0,2,2],psym=[-2,-4,-5,-6],thick=[2,1,2,1],/bottom
            xyouts,[4000],[10500],[textoidl('10^{10}M')+sunsymbol()],charsize = 2,/device
        endif else begin
            legend,[textoidl('10^5'),textoidl('10^4'),textoidl('10^3'),textoidl('10^2')],linestyle=[0,0,2,2],psym=[-2,-4,-5,-6],thick=[2,1,2,1],/right
            xyouts,[4000],[10500],[textoidl('10^{10}M')+sunsymbol()],charsize = 2,/device
        endelse
        device,/close
ENDFOR
END
