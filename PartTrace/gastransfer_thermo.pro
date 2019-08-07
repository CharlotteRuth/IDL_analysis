PRO MASTER

dir = '/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK'
filebase = 'h986.cosmo50cmb.3072g14HBWK'
step = '00512'
haloida = '11'
haloidb = '27'

;dir='/nobackupp8/crchrist/MolecH/rogue.cosmo25cmb.4096g/rogue.cosmo25cmb.4096g5HbwK1BH/'
;filebase='rogue.cosmo25cmb.4096g5HbwK1BH'
;haloida = '1'
;haloidb = '32' ;'34'
;step = '004096'

;dir = '/nobackupp8/crchrist/MolecH/elektra.cosmo25cmb.4096g/elektra.cosmo25cmb.4096g5HbwK1BH/'
;filebase='elektra.cosmo25cmb.4096g5HbwK1BH'
;haloida = '1'
;haloidb = '11'
;step = '004096'

;dir = '/nobackupp8/crchrist/MolecH/storm.cosmo25cmb.4096g/storm.cosmo25cmb.4096g5HbwK1BH/' 
;filebase='storm.cosmo25cmb.4096g5HbwK1BH' 
;haloida = '1'
;haloidb = '35'
;step = '004096'

gastransfer_thermo,dir,filebase,step,haloida,haloidb

END

PRO gastransfer_thermo,dir,filebase,step,haloida,haloidb
cd,dir
loadct,39
rtipsy,filebase+'.'+step+'/'+filebase+'.'+step,h,g,d,s,/justhead
iord = read_lon_array(filebase+'.'+step+'/'+filebase+'.'+step+'.iord')
amigastat = read_stat_struc_amiga(filebase+'.'+step+'/'+filebase+'.'+step+'.amiga.stat')
units = tipsyunits(filebase+'.param')

print,'1, grp'+haloida
allgas1 = mrdfits('grp'+haloida+'.allgas.iord.fits')
match,iord,allgas1,ind1_tip,ind1
readcol,filebase + '.grp'+haloida+'.haloid.dat',files1,pasthalos1,format='(A,A)'
n_steps = n_elements(files1)
dtime = 13.7/fix((strsplit(files1[0],'.',/extract))[-1])
time1 = fltarr(n_elements(files1))
FOR is = 0, n_elements(time1) - 1 DO time1[is] = dtime*fix((strsplit(files1[is],'.',/extract))[-1])

allgas2 = mrdfits('grp'+haloidb+'.allgas.iord.fits') &  $
print,'2, grp'+haloidb & $
match,iord,allgas2,ind2_tip,ind2 
match,allgas2,allgas1,temp1,temp2 

rtipsy,filebase+'.'+step+'/'+filebase+'.'+step,h,g,d,s
match,iord,allgas1,iordind1,temp
match,iord,allgas2,iordind2,temp

IF temp1[0] NE -1 THEN print,'2-1',n_elements(temp1),n_elements(allgas2),n_elements(allgas1) ELSE RETURN
readcol,filebase + '.grp'+haloidb+'.haloid.dat',files2,pasthalos2,format='(A,A)'
time2 = fltarr(n_elements(files1))
FOR is = 0, n_elements(time2) - 1 DO time2[is] = dtime*fix((strsplit(files2[is],'.',/extract))[-1])

data1 = mrdfits('grp'+haloida+'.allgas.entropy.fits',1)
data2 = mrdfits('grp'+haloidb+'.allgas.entropy.fits',1)
match,allgas1,allgas2,temp1,temp2
data12 = data1[temp1]
FOR is = n_steps - min([n_elements(pasthalos1),n_elements(pasthalos2)]), n_steps - 1 DO BEGIN &  $
    IF (where(data12[*].grp[is] EQ pasthalos1[n_elements(pasthalos1) - 1 - is]))[0] NE -1 THEN data12[where(data12[*].grp[is] EQ pasthalos1[n_elements(pasthalos1) - 1 - is])].grp[is] = long(haloida) &  $
    IF (where(data12[*].grp[is] EQ pasthalos2[n_elements(pasthalos2) - 1 - is]))[0] NE -1 THEN data12[where(data12[*].grp[is] EQ pasthalos2[n_elements(pasthalos2) - 1 - is])].grp[is] = long(haloidb) &  $
    IF (where((data12[*].grp[is] NE long(haloida)) AND $
              (data12[*].grp[is] NE long(haloidb)) AND $
              (data12[*].grp[is] NE 0)))[0] NE -1 THEN $
              data12[where( $
              data12[*].grp[is] NE long(haloida) AND $
              data12[*].grp[is] NE long(haloidb) AND $
              data12[*].grp[is] NE 0)].grp[is] = -1 &  $
  ENDFOR

data12colors = data12.grp
data12colors[where(data12.grp EQ -1)] = 60 ;Different halo
data12colors[where(data12.grp EQ 0)] = 20 ;Not in halo
data12colors[where(data12.grp EQ long(haloida))] = 254 ;halo 1
data12colors[where(data12.grp EQ long(haloidb))] = 200 ;halo 2

window,0
plot,[0,0],[0,0],xrange = minmax(data12.x),yrange = minmax(data12.y),/nodata
FOR ip = 0, n_elements(data12.iord) - 1 DO $ ;oplot,data12[ip].x,data12[ip].y,color = ip
;  FOR is = 0, n_steps -2 DO $
  FOR is = n_steps - 12, n_steps -2 DO $
    oplot,[data12[ip].x[is],data12[ip].x[is+1]],[data12[ip].y[is],data12[ip].y[is+1]],color = data12colors[is,ip],psym = 3 ;ip
FOR ip = 0, n_elements(data12.iord) - 1 DO $
  oplot,[(data12.x)[n_elements(data12[0].x) - 1,ip],(data12.x)[n_elements(data12[0].x) - 1,ip]],[(data12.y)[n_elements(data12[0].x) - 1,ip],(data12.y)[n_elements(data12[0].x) - 1,ip]],psym = 4,color = data12colors[n_elements(data12[0].x) - 1,ip]
vir = circle(amigastat[where(amigastat.group EQ haloidb)].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloidb)].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloidb)].rvir)
oplot,vir[0,*],vir[1,*],color = 120
vir = circle(amigastat[where(amigastat.group EQ haloida)].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloida)].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloida)].rvir)
oplot,vir[0,*],vir[1,*],color = 120
legend,[haloida,haloidb],color = [254,200],psym = [4,4]
stop

window,1
plot,[0,0],[0,0],xrange = [9,max(time1)],yrange = [1e-6,1e4],/ylog,xtitle = 'Time [Gyr]',ytitle = 'Density [amu/cc]'
FOR ip = 0, n_elements(data12.iord) - 1 DO $
;  FOR is = 0, n_steps -2 DO $
  FOR is = n_steps - 12, n_steps -2 DO $
    oplot,[(reverse(time1))[is],(reverse(time1))[is + 1]],[data12[ip].rho[is],data12[ip].rho[is+1]],color = data12colors[is,ip] ;ip
stop

;for Anna
startstep = 60
ind = where((data12.grp)[startstep:n_steps - 1,*] eq haloidb)
s = size((data12.grp)[startstep:n_steps - 1,*])
ncol = s(1)
ind_step = (ind MOD ncol) + startstep
ind_part = ind / ncol
window,1
plot,[0,0],[0,0],xrange = [9,max(time1)],yrange = [1e-6,1e4],/ylog,xtitle = 'Time [Gyr]',ytitle = 'Density [amu/cc]'
FOR ip = 0, n_elements(ind_part) - 1 DO $
;  FOR is = 0, n_steps -2 DO $
  FOR is = startstep, n_steps - 2 DO $
    oplot,[(reverse(time1))[is],(reverse(time1))[is + 1]],[data12[ind_part[ip]].rho[is],data12[ind_part[ip]].rho[is+1]],color = data12colors[is,ind_part[ip]] ;ip
END
