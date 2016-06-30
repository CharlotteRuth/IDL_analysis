PRO SN_FB_H2_master
;dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g8HBWK/steps/h516.cosmo25cmb.1536g9HBWK.00192.dir/'
dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g8HBWK/steps.11/h516.cosmo25cmb.1536g11HBWK.00408.dir/'
;file = 'h516.cosmo25cmb.1536g9HBWK.00192.halo.1'
file = 'h516.cosmo25cmb.1536g11HBWK.00408.halo.1'
;pfile = '../../h516.cosmo25cmb.1536g8HBWK.param'
pfile = '../../h516.cosmo25cmb.1536g9HBWK.param'
;sfile = '../../h516.cosmo25cmb.1536g9HBWK.starlog'
sfile = '../../h516.cosmo25cmb.1536g11HBWK.starlog'
;snaptime = 1.343671e-1
SN_FB_H2,dir,file,pfile,sfile;,snaptime = snaptime
END

PRO SN_FB_H2,dir,file,pfile,sfile,snaptime = snaptime
loadct,39
cm_in_kpc = 3.08568025e21

cd,dir
units = tipsyunits(pfile)
rtipsy,file,h,g,d,s
IF NOT KEYWORD_SET(snaptime) then snaptime = MAX(s.tform)
snaptime = snaptime*units.timeunit
readarr,file+'.H2',h,H2,part = 'gas',/ascii
readarr,file+'.HI',h,HI,part = 'gas',/ascii
Hall = 2.0*H2 + HI
;readarr,file+'.smoothlength',h,smooth,part = 'gas',/ascii
readarr,file+'.correL',h,correL,/ascii,part = 'gas'
correL = correL*units.lengthunit*cm_in_kpc
readarr,file+'.coolontime',h,coolontime,part = 'gas',/ascii
coolontime = coolontime*units.TIMEUNIT
slog = rstarlog(sfile,/MOLECULARH)
readarr,file+'.iord',h,siord,part = 'star',/ascii
readarr,file+'.iord',h,giord,part = 'gas',/ascii
a3 = h.time*h.time*h.time

ind_arr = [0]
FOR is = 0L, N_ELEMENTS(siord) - 1 DO BEGIN
    ind = WHERE(siord[is] eq slog.iorderstar)
    if (ind[0] ne -1) then ind_arr = [ind_arr,ind]
ENDFOR
slog = slog[ind_arr[1:N_ELEMENTS(ind_arr) - 1]]

ind_arr_g = [0]
ind_arr_slog = [-1]
FOR ig = 0L, N_ELEMENTS(slog) - 1 DO BEGIN
    ind = WHERE(slog[ig].iordergas eq giord)
    if (ind[0] ne -1) then begin
        if N_ELEMENTS(giord[ind] eq slog[ig].iordergas) gt 1 then stop
        ind_arr_g = [ind_arr_g,ind]
        ind_arr_slog = [ind_arr_slog,ig]
    endif
ENDFOR
gstar = g[ind_arr_g[1:N_ELEMENTS(ind_arr_g) - 1]]
ind_arr_slog = ind_arr_slog[1:N_ELEMENTS(ind_arr_slog) - 1]

!p.multi = 0
window,3
indco = where(coolontime ne 0, complement = nindco)
ycoolontime = histogram(coolontime[indco],nbins = 100,locations = xcoolontime,min = 0,max = 5.5e9)
plot,xcoolontime,ycoolontime,psym = 10
g_coolontime = coolontime[ind_arr[1:N_ELEMENTS(ind_arr) - 1]]
indco_g = where(g_coolontime ne 0, complement = nindco_g)
ycoolontime = histogram(g_coolontime[indco_g],nbins = 100,locations = xcoolontime,min = 0,max = 5.5e9)
oplot,xcoolontime,ycoolontime,psym = 10,color = 240
oplot,[snaptime,snaptime],[0,6000]

mincool = 0 ;4.7e9;MIN(coolontime[ind])
maxcool = MAX(coolontime[indco] - snaptime)
color = 240*(coolontime - snaptime - mincool)/(maxcool - mincool) + 14
;mincool = MIN(alog10(coolontime[indco]))
;maxcool = MAX(alog10(coolontime[indco]))
;color = 240*(alog10(coolontime) - mincool)/(maxcool - mincool)
base = 255
color[where(color lt 14)] = base
coolon = where(coolontime - snaptime lt 0, complement = cooloff)

mincoolH2 = -6 ;4.7e9;MIN(coolontime[ind])
maxcoolH2 = 0
colorH2 = 240*(alog10(2.0*H2/Hall) - mincoolH2)/(maxcoolH2 - mincoolH2) + 14
colorH2[where(colorH2 lt 14)] = base

window,0,xsize = 712,ysize = 712
!p.multi =[0,2,2]
plot,g.dens*units.rhounit/a3,2.0*H2/Hall,psym = 3,/xlog,/ylog,xrange = [1e-6,1e4],yrange = [1e-10,1]
FOR i = 0L, N_ELEMENTS(g) - 1 DO $
  if(color[i] ne base) then oplot,[g[i].dens,g[i].dens]*units.rhounit/a3,2.0*[H2[i]/Hall[i],H2[i]/Hall[i]],psym = 2,color = color[i],symsize = 0.5

plot,g.dens*units.rhounit/a3,2.0*H2/Hall,psym = 3,/xlog,/ylog,xrange = [1e-6,1e4],yrange = [1e-6,1]
oplot,slog.rhoform*units.rhounit,slog.H2form,color = 240, psym = 3

plot,g[coolon].dens*units.rhounit/a3*correL[coolon],2.0*H2[coolon]/Hall[coolon],psym = 3,/xlog,/ylog,xrange = [1e18,1e24],yrange = [1e-6,1]
oplot,g[cooloff].dens*units.rhounit/a3*correL[cooloff],2.0*H2[cooloff]/Hall[cooloff],psym = 3,color = 50

plot,g[coolon].dens*units.rhounit/a3,2.0*H2[coolon]/Hall[coolon],psym = 3,/xlog,/ylog,xrange = [1e-6,1e4],yrange = [1e-6,1]
oplot,g[cooloff].dens*units.rhounit/a3,2.0*H2[cooloff]/Hall[cooloff],color = 50,psym = 3
oplot,slog.rhoform*units.rhounit*correL,slog.H2form,color = 240, psym = 3
;plot,g.dens*units.rhounit/a3,g.tempg,psym = 3,/xlog,/ylog,xrange = [1e-6,1e4],yrange = [10,1e6]
;FOR i = 0L, N_ELEMENTS(g) - 1 DO $
;  if(color[i] ne base) then oplot,[g[i].dens,g[i].dens]*units.rhounit/a3,[g[i].tempg,g[i].tempg],psym = 2,color = color[i],symsize = 0.5

;plot,g.dens*units.rhounit/a3,g.tempg,psym = 3,/xlog,/ylog,xrange = [1e-6,1e4],yrange = [10,1e6]
;oplot,slog.rhoform*units.rhounit,slog.tempform,color = 240, psym = 3


window,5,xsize = 712,ysize = 712
!p.multi = [0,2,2]
Hall = 2.0*H2 + HI
indh = where(ABS(g.z*units.lengthunit*h.time) lt 0.75)
;plot,g.x*units.lengthunit*h.time,g.y*units.lengthunit*h.time,psym = 3,xrange = [-10,10],yrange = [-10,10]
;oplot,g[indco].x*units.lengthunit*h.time,g[indco].y*units.lengthunit*h.time,psym = 3,color = 240
;plot,g[indco].x*units.lengthunit*h.time,g[indco].y*units.lengthunit*h.time,psym
;= 3,xrange = [-3,3],yrange = [-3,3]
plot,g[indh].x*units.lengthunit*h.time,g[indh].y*units.lengthunit*h.time,psym = 3,xrange = [-4,4],yrange = [-4,4]
FOR i = 0L, N_ELEMENTS(g) - 1 DO $
  if(color[i] ne base AND (ABS(g[i].z*units.lengthunit) lt 0.75)) then oplot,[g[i].x,g[i].x]*units.lengthunit*h.time,[g[i].y,g[i].y]*units.lengthunit*h.time,psym = 2,color = color[i],symsize = 0.5
;
plot,g[indco].y*units.lengthunit*h.time,g[indco].z*units.lengthunit*h.time,psym = 3,xrange = [-3,3],yrange = [-3,3]
FOR i = 0L, N_ELEMENTS(g) - 1 DO $
  if(color[i] ne base) then oplot,[g[i].y,g[i].y]*units.lengthunit*h.time,[g[i].z,g[i].z]*units.lengthunit*h.time,psym = 2,color = color[i],symsize = 0.5

plot,g[indh].x*units.lengthunit*h.time,g[indh].y*units.lengthunit*h.time,psym = 3,xrange = [-4,4],yrange = [-4,4]
FOR i = 0L, N_ELEMENTS(g) - 1 DO $
  if(colorH2[i] ne base AND (ABS(g[i].z*units.lengthunit) lt 0.75)) then oplot,[g[i].x,g[i].x]*units.lengthunit*h.time,[g[i].y,g[i].y]*units.lengthunit*h.time,psym = 3,color = colorH2[i]
;
plot,g[indco].y*units.lengthunit*h.time,g[indco].z*units.lengthunit*h.time,psym = 3,xrange = [-3,3],yrange = [-3,3]
FOR i = 0L, N_ELEMENTS(g) - 1 DO $
  if(colorH2[i] ne base) then oplot,[g[i].y,g[i].y]*units.lengthunit*h.time,[g[i].z,g[i].z]*units.lengthunit*h.time,psym = 3,color = colorH2[i]


!p.multi = 0

window,6
plot,coolontime,g.tempg,psym = 3,/ylog
oplot,coolontime[ind_arr_g],gstar.tempg,psym = 3,color = 240
oplot,coolontime[ind_arr_g],slog[ind_arr_slog].tempform,psym = 3,color = 50

plot,coolontime[ind_arr_g],slog[ind_arr_slog].timeform*units.timeunit,psym = 3
plot,slog[ind_arr_slog].timeform*units.timeunit,coolontime[ind_arr_g]-slog[ind_arr_slog].timeform*units.timeunit,psym = 3

stop
END
