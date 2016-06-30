PRO starform_analysis_master

;dir = '/astro/net/nbody1/jillian/h277/evenlower/'
;file = 'h277.cosmo50cmb.1536g1MBWKBH.00512/h277.cosmo50cmb.1536g1MBWKBH.00512'
;sfile = 'h277.cosmo50cmb.1536g1MBWKBH.starlog'
;pfile = 'h277.cosmo50cmb.1536g1MBWKBH.param'
;outfile = '/astro/net/scratch2/christensen/SFresults/h277.1536g1MBWKBH.'

;dir = '/net/nbody1/fabio/RUNS/h277.cosmo50cmb.3072gs1MbwK/'
;file = 'h277.cosmo50cmb.3072gs1MbwK.00128'
;sfile = 'h277.cosmo50cmb.3072gs1MbwK.starlog'
;pfile = 'h277.cosmo50cmb.3072gs1MbwK.param'
;outfile = '/astro/net/scratch2/christensen/SFresults/h277.3072gs1MbwK.'

;dir = '/astro/net/nbody1/abrooks/h603.cosmo50cmb.3072gs1MbwK/'
;file = 'h603.cosmo50cmb.3072gs1MbwK.00512/h603.cosmo50cmb.3072gs1MbwK.00512'
;sfile = 'h603.cosmo50cmb.3072gs1MbwK.starlog'
;pfile = 'h603.cosmo50cmb.3072gs1MbwK.param'
;outfile = '/astro/net/scratch2/christensen/SFresults/h603.3072gs1MbwK.'

;dir = '/astro/net/nbody1/abrooks/h516.cosmo25cmb.3072g1MBWK/'
;file = 'h516.cosmo25cmb.3072g1MBWK.00512/h516.cosmo25cmb.3072g1MBWK.00512'
;sfile = 'h516.cosmo25cmb.3072g1MBWK.starlog'
;pfile = 'h516.cosmo25cmb.3072g1MBWK.param'
;outfile = '/astro/net/scratch2/christensen/SFresults/h516.3072g1MBWK'

;dir = '/astro/net/scratch2/christensen/SFresults/h285.cosmo50cmb.3072gs1MbwK.00216.dir/'
;file = 'h285.cosmo50cmb.3072gs1MbwK.00216'
;sfile = 'h285.cosmo50cmb.3072gs1MbwK.starlog'
;pfile = 'h285.cosmo50cmb.3072gs1MbwK.param'
;outfile = '/astro/net/scratch2/christensen/SFresults/h285.3072gs1MbwK'

dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/'
file = 'steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512'
sfile = 'h516.cosmo25cmb.1536g3HBWK.starlog'
pfile = 'h516.cosmo25cmb.1536g3HBWK.param'
outfile = '/astro/net/scratch2/christensen/SFresults/h516.1536g3HBWK'


starform_analysis,dir,file,sfile,pfile,outfile = outfile
END


PRO starform_analysis,dir,file,sfile,pfile,outfile = outfile
loadct,39
!P.CHARTHICK=1.5
!X.THICK=1.5
!Y.THICK=1.5
!p.charsize=1.0
!x.charsize=1.5
!y.charsize=1.5  
set_plot,'x'

units = tipsyunits(dir + pfile)
rtipsy,dir + file,h,g,d,s
slog = rstarlog(dir + sfile)
readarr,dir+file+'.iord',h,siord,part = 'star',/ascii
;stop
;ind_arr = [0]
;FOR is = 0L, N_ELEMENTS(siord) - 1 DO BEGIN
;    ind = WHERE(siord[is] eq slog.iorderstar)
;    if ind then ind_arr = [ind_arr,ind]
;ENDFOR
;slog = slog[ind_arr[1:N_ELEMENTS(ind_arr) - 1]]
slog = slog[siord - min(siord)]

window,0
plot,g.dens*units.rhounit/h.time/h.time/h.time,g.tempg,psym = 3,/xlog,/ylog,xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]',title = file,xrange = [1e-5,1e5],yrange = [10,1e7]
oplot,slog.rhoform*units.rhounit,slog.tempform,psym = 3, color = 240

window,1
nstars = histogram(alog10(slog.tempform),nbins = 100,locations = tempform, min = 1, max = 4.5)
plot,tempform,nstars,psym = 10,xrange = [1,4.5],xtitle = 'SF Temperature [K]',ytitle = 'Number of Stars' 
dtime = 0.05
FOR it = 0, CEIL(MAX(s.tform)/dtime - 1) DO BEGIN
    indstars = where(s.tform ge it*dtime and s.tform lt (it + 1)*dtime)
    IF (indstars[0] ne -1) THEN BEGIN
        nstars = histogram(alog10(slog[indstars].tempform),nbins = 100,locations = tempform, min = 1, max = 4.5)
        oplot,tempform,nstars,psym = 10,color = 240/CEIL(MAX(s.tform)/dtime)*(it + 1)      
    ENDIF
ENDFOR

window,2
plot,s.metals,slog.tempform,psym = 3,/xlog,/ylog,xtitle = 'Stellar Metallicity',ytitle = 'SF Temperature [K]',title = file

window,3
sfr,s,title = file,timeunit = units.timeunit,massunit = units.massunit ;300
IF ((where(slog.tempform gt 6000 and slog.tempform lt 10000))[0] ne -1) THEN $
  sfr,s[where(slog.tempform gt 6000 and slog.tempform lt 10000)],/overplot,color = 240,timeunit = units.timeunit,massunit = units.massunit
IF ((where(slog.tempform  gt 10000))[0] ne -1 ) THEN $
  sfr,s[where(slog.tempform  gt 10000)],/overplot,color = 50,timeunit = units.timeunit,massunit = units.massunit
legend,['SF Temperature > 60000 and < 10^4','SF Temperatre > 10^4'],color = [240,50],linestyle = [0,0],/right

stop

;********************************************
IF KEYWORD_SET(outfile) THEN BEGIN
    set_plot,'ps' 
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0

;    outfile = '/astro/net/scratch2/christensen/SFresults/'

    device,filename=outfile+'SF_RhoTemp.eps',/color,bits_per_pixel= 8,/times
    plot,g.dens*units.rhounit/h.time/h.time/h.time,g.tempg,psym = 3,/xlog,/ylog,xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]',title = file,xrange = [1e-5,1e5],yrange = [10,1e7]
    oplot,slog.rhoform*units.rhounit,slog.tempform,psym = 3, color = 240
    device,/close

    device,filename=outfile+'SF_TempDist.eps',/color,bits_per_pixel= 8
    nstars = histogram(alog10(slog.tempform),nbins = 100,locations = tempform, min = 1, max = 4.5)
    plot,tempform,nstars,psym = 10,xrange = [1,4.5],xtitle = 'LOG(SF Temperature) [K]',ytitle = 'Number of Stars',thick = 2,xstyle = 1,title = file
    oplot,[alog10(6000),alog10(6000)],[0,MAX(nstars)]
    FOR it = 0, CEIL(MAX(s.tform)/dtime - 1) DO BEGIN
        indstars = where(s.tform ge it*dtime and s.tform lt (it + 1)*dtime)
        IF (indstars[0] ne -1) THEN BEGIN
            nstars = histogram(alog10(slog[indstars].tempform),nbins = 100,locations = tempform, min = 1, max = 4.5)
            oplot,tempform,nstars,psym = 10,color = 240/CEIL(MAX(s.tform)/dtime)*(it + 1),thick = 2
        ENDIF
    ENDFOR
    device,/close

    device,filename=outfile+'SF_ZTemp.eps',/color,bits_per_pixel= 8
    plot,s.metals,slog.tempform,psym = 3,/xlog,/ylog,xtitle = 'Stellar Metallicity',ytitle = 'SF Temperature [K]',title = file
    device,/close

    device,filename=outfile+'SFH_temp.eps',/color,bits_per_pixel= 8
    sfr,s,title = file,timeunit = units.timeunit,massunit = units.massunit ;300
IF ((where(slog.tempform gt 6000 and slog.tempform lt 10000))[0] ne -1) THEN $
  sfr,s[where(slog.tempform gt 6000 and slog.tempform lt 10000)],/overplot,color = 240,timeunit = units.timeunit,massunit = units.massunit
IF ((where(slog.tempform  gt 10000))[0] ne -1 ) THEN $
  sfr,s[where(slog.tempform  gt 10000)],/overplot,color = 50,timeunit = units.timeunit,massunit = units.massunit
legend,['SF Temperature > 60000 and < 10^4','SF Temperatre > 10^4'],color = [240,50],linestyle = [0,0],/right
    set_plot,'x'
ENDIF

;********************************************************************

;readcol,dir+file+'.iord',array
tempformarray = fltarr(h.n)
tempformarray[h.n - h.nstar:h.n - 1] = slog.tempform

openw,1,dir+file+'.tempform'
printf,1,LONG(h.n)
close,1
writecol,dir+file+'.tempform',format = '',tempformarray,/update

IF (1) THEN BEGIN
s2 = s
s2.x = 2.0*((alog10(slog.rhoform) - MIN(alog10(slog.rhoform)))/(MAX(alog10(slog.rhoform)) - MIN(alog10(slog.rhoform)))) - 1.0
s2.y = 2.0*((alog10(slog.tempform) - MIN(alog10(slog.tempform)))/(MAX(alog10(slog.tempform)) - MIN(alog10(slog.tempform)))) - 1.0
s2.z = 2.0*((s.tform)/(MAX(s.tform))) - 1.0
IF ((where(s.tform lt 0))[0] ne -1) THEN s2[where(s.tform lt 0)].z = 1.0

wtipsy,dir+file+'_starform',h,g,d,s2,/standard
ENDIF

END
