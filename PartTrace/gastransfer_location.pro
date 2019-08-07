;

PRO gastransfer_location_MASTER

dir = '/nobackupp8/crchrist/MolecH/cptmarvel.cosmo25.4096g/cptmarvel.cosmo25cmb.4096g5HbwK1BH'
filebase = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
step = '004096'
haloids = ['1','2','4','5','6','7','10','11','13','14']

dir = '/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK'
filebase = 'h986.cosmo50cmb.3072g14HBWK'
step = '00512'
;haloids = ['1','2','3','8', '15','16'] ;11, 14
haloids = ['11','27']

;dir = '/nobackupp8/crchrist/MolecH/elektra.cosmo25cmb.4096g/elektra.cosmo25cmb.4096g5HbwK1BH/'
;filebase='elektra.cosmo25cmb.4096g5HbwK1BH'
;haloids=['1','2','3','4','5','8','9','10','11','12','17','18']
;haloids_int = ['1','5','11']
;step = '004096'

;dir='/nobackupp8/crchrist/MolecH/rogue.cosmo25cmb.4096g/rogue.cosmo25cmb.4096g5HbwK1BH/'
;filebase='rogue.cosmo25cmb.4096g5HbwK1BH'
;haloids = ['1','3','7','8','10','11','12','16','17','18','30','31','32','34','36']
;haloids_int = ['1','3','7','8','10','32','34']
;step = '004096'

;dir = '/nobackupp8/crchrist/MolecH/storm.cosmo25cmb.4096g/storm.cosmo25cmb.4096g5HbwK1BH/' 
;filebase='storm.cosmo25cmb.4096g5HbwK1BH' 
;haloids= ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60']
;haloids_int = ['1','4','8','35','43','49']
;step = '004096'
gastransfer_location,dir,filebase,step,haloids
END

PRO gastransfer_location,dir,filebase,step,haloids
formatplot

cd,dir
colors = (indgen(n_elements(haloids)) + 1)*(254/(n_elements(haloids)))
n_steps = n_elements(files1)

rtipsy,filebase+'.'+step+'/'+filebase+'.'+step,h,g,d,s
units = tipsyunits(filebase+'.param')
g.x = g.x*h.time*units.lengthunit
g.y = g.y*h.time*units.lengthunit
g.z = g.z*h.time*units.lengthunit
s.x = s.x*h.time*units.lengthunit
s.y = s.y*h.time*units.lengthunit
s.z = s.z*h.time*units.lengthunit
loadct,39
xminmax = [-0.1,0.08]*h.time*units.lengthunit
yminmax = [-0.12,0.02]*h.time*units.lengthunit
zminmax = [-0.02,0.1]*h.time*units.lengthunit
print,xminmax
print,yminmax
print,zminmax
iord = read_lon_array(filebase+'.'+step+'/'+filebase+'.'+step+'.iord')
amigastat = read_stat_struc_amiga(filebase+'.'+step+'/'+filebase+'.'+step+'.amiga.stat')

window,2,xsize = 1200,ysize = 600
multiplot,[2,1],/square
plot,s.x,s.y,psym = 3;,xrange = xminmax,yrange = yminmax,xstyle = 1,ystyle = 1,/nodata
oplot,s.x,s.y,psym = 3,color = 190
FOR ihalos = 0, n_elements(haloids) - 1 DO BEGIN
    allgas1 = mrdfits('grp'+haloids[ihalos]+'.allgas.iord.fits')
    match,iord,allgas1,ind1_tip,ind1
    oplot,g[ind1_tip].x,g[ind1_tip].y,psym = 3,color = colors[ihalos]
    vir = circle(amigastat[where(amigastat.group EQ haloids[ihalos])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[ihalos])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[ihalos])].rvir)
    oplot,vir[0,*],vir[1,*]
ENDFOR
oplot,s.x,s.y,psym = 3,color = 190
multiplot
plot,s.x,s.z,psym = 3;,xrange = xminmax,yrange = zminmax,xstyle = 1,ystyle = 1,/nodata
oplot,s.x,s.z,psym = 3,color = 190
FOR ihalos = 0, n_elements(haloids) - 1 DO BEGIN
    allgas1 = mrdfits('grp'+haloids[ihalos]+'.allgas.iord.fits')
    match,iord,allgas1,ind1_tip,ind1
    oplot,g[ind1_tip].x,g[ind1_tip].z,psym = 3,color = colors[ihalos]
    vir = circle(amigastat[where(amigastat.group EQ haloids[ihalos])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[ihalos])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[ihalos])].rvir)
    oplot,vir[0,*],vir[1,*]
ENDFOR
oplot,s.x,s.z,psym = 3,color = 190
multiplot,/reset

END
