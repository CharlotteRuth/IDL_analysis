pro testSplatter


cd,'/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/1E5R/12M_k'
file = 'o12M_1.00300'
dMsolUnit       = 2.362e5
dKpcUnit        = 1
rtipsy,file,h,g1,d1,s1
g = g1
;g.x = -1.0*g.x
g.mass = g1.mass*dMsolUnit

data = read_cube_fits('data.fits',header,kpcunit = dKpcUnit,munit = dMsolUnit)
data1 = data[*,*,1]
xaxes = (findgen(header.NAXIS1) - header.NAXIS1/2.0)*header.CDELT1
yaxes = (findgen(header.NAXIS2) - header.NAXIS1/2.0)*header.CDELT2

ncell = header.naxis1
minrange = header.crval1
res = header.cdelt1
gas_sd = fltarr(ncell,ncell)
print,'IDL:   ',total(g.mass)
print,'Tipsy: ',total(data1)
;FOR ix = 0, ncell - 1 DO BEGIN
;    FOR iy = 0,ncell - 1 DO BEGIN 
;        indg = where((g.x ge (res*ix + minrange)) AND (g.x lt (res*(ix+1) + minrange)) AND (-1.0*g.y ge (res*iy + minrange)) AND (-1.0*g.y lt (res*(iy+1) + minrange)))
;        IF indg[0] NE -1 THEN gas_sd[ix,iy] = TOTAL(g[indg].mass) ELSE gas_sd[ix,iy] = 0
;        stop
;     ENDFOR
; ENDFOR
FOR ix = 0, ncell - 1 DO BEGIN
    FOR iy = 0,ncell - 1 DO BEGIN 
        indg = where((g.x ge (res*ix + minrange)) AND (g.x lt (res*(ix+1) + minrange))$
                 AND (g.y ge (res*iy + minrange)) AND (g.y lt (res*(iy+1) + minrange)))
        IF indg[0] NE -1 THEN gas_sd[ix,iy] = TOTAL(g[indg].mass) ELSE gas_sd[ix,iy] = 0
;        IF ix eq 8 AND iy eq 8 THEN STOP
 ;      IF ix eq 9 AND iy eq 10 THEN STOP
 ;       IF ix eq 8 AND iy eq 10 THEN STOP
 ;       IF ix eq 10 AND iy eq 10 THEN STOP
;       IF (WHERE(indg eq 0))[0] ne -1 then stop
;        stop
     ENDFOR
 ENDFOR
data2 = rotate(data1,2)


window,0,xsize = 400,ysize = 400
contour,alog10(gas_sd),xaxes,yaxes,/fill,nlevels = 60,yrange = [-50,50],xrange = [-50,50],min_value = 8,max_value = 10,xstyle = 1,ystyle = 1
;oplot,g.x,-1.0*g.y,psym = 2
window,1,xsize = 400,ysize = 400
contour,alog10(data1),xaxes,yaxes,/fill,nlevels = 60,yrange = [-50,50],xrange = [-50,50],min_value = 8,max_value = 10,xstyle = 1,ystyle = 1
;oplot,g.x,-1.0*g.y,psym = 2
window,2
plot,data1,gas_sd/data1,psym = 2,/ylog,/xlog;,xrange = [1e4,1e10]
stop

end
