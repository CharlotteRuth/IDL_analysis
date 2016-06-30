PRO align_over_time
device,decomposed = 0
step = '00120'

dir = '/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK'
filename = 'h516.cosmo25cmb.3072g14HBWK'
haloid = '2';1,2

dir = '/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK'
filename = 'h799.cosmo25cmb.3072g14HBWK'
haloid = '6';1,4,6

dir = '/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK'
filename = 'h986.cosmo50cmb.3072g14HBWK'
haloid = '8';1,2,3,8

dir = '/nobackupp8/crchrist/MolecH/h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK'
filename = 'h285.cosmo50cmb.3072g14HMbwK'
haloid = '1';1,4,9

dir = '/nobackupp8/crchrist/MolecH/h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK'
filename = 'h258.cosmo50cmb.3072g14HMbwK'
haloid = '1';1,4,8

dir = '/nobackupp8/crchrist/MolecH/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK'
filename = 'h239.cosmo50cmb.3072g14HMbwK'
haloid = '1'

cd,dir
units = tipsyunits(filename + '.param')

sl = rstarlog(filename + '.starlog',/molecularH)
;sl = mrdfits(filename + '_2merge.starlog.fits',1)
halodat = mrdfits('grp'+haloid+'.alignment.fits',1)

;readcol,filename + '.log',dTime,z,E,T,U,Eth,Lx,Ly,Lz,WallTime,dWMax,dImax,dEMax,dMultiEff,/silent
;readcol,'/nobackupp8/crchrist/MolecH/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/h239.cosmo50cmb.3072g14HMbwK.log',dTime,z,E,T,U,Eth,Lx,Ly,Lz,WallTime,dWMax,dImax,dEMax,dMultiEff,/silent
;zfix = fltarr(n_elements(halodat))
;tfix = fltarr(n_elements(halodat))
;FOR i = 0,n_elements(halodat) -1 DO BEGIN
;   temp = min(abs(z - halodat[i].z),indz)
;   zfix[i] = z[indz]
;   tfix[i] = dtime[indz]
;   print,zfix[i],halodat[i].z,tfix[i],halodat[i].time*1e9/units.timeunit
;END
;halodat.time = tfix;halodat.time*1e9/units.timeunit
;halodat.z = zfix
;writecol,'grp'+haloid+'.alignment.txt',halodat.file,halodat.haloid,tfix,zfix,halodat.xc,halodat.yc,halodat.zc,halodat.vxc,halodat.vyc,halodat.vzc,halodat.xa,halodat.ya,halodat.za,format = '(A,I,F,F,F,F,F,F,F,F,F,F,F)'
halodat.time = halodat.time*1e9/units.timeunit
a = 1/(1 + halodat.z)

rtipsy,filename + '.' + step + '/' + filename + '.' + step,h,g,d,s
read_tipsy_arr,filename + '.' + step + '/' + filename + '.' + step + '.iord',h,siord,part = 'star',type = 'long'
ind_tcurr = where(halodat.file EQ filename + '.' + step + '/' + filename + '.' + step)
tcurr = halodat[ind_tcurr].time;max(s.tform)
;srecent = s[where(s.tform GT tcurr - 0.1*1e9/units.timeunit)]
siord_recent = siord[where(s.tform GT tcurr - 0.1*1e9/units.timeunit)]

timeform_uniq = sl[uniq(sl.timeform,sort(sl.timeform))].timeform

xcfit = spline(halodat.time,halodat.xc,timeform_uniq) ;comoving
ycfit = spline(halodat.time,halodat.yc,timeform_uniq) ;comoving
zcfit = spline(halodat.time,halodat.zc,timeform_uniq) ;comoving
afit = spline(halodat.time,a,timeform_uniq)

match2,sl.timeform,timeform_uniq,ind1,ind2
slx = (sl.x - xcfit[ind1]);*afit[ind1]);comoving
sly = (sl.y - ycfit[ind1]);*afit[ind1]);comoving
slz = (sl.z - zcfit[ind1]);*afit[ind1]);comoving
afit = afit[ind1]

rmax = 10/units.lengthunit
indgal = where(sqrt(slx*slx + sly*sly + slz*slz)*afit LT rmax)
sl = sl[indgal]
slx = slx[indgal]
sly = sly[indgal]
slz = slz[indgal]
;indgal = where(sqrt((srecent.x - (halodat[ind_tcurr].xc)[0])*(srecent.x - (halodat[ind_tcurr].xc)[0]) + (srecent.y -  (halodat[ind_tcurr].yc)[0])*(srecent.y - (halodat[ind_tcurr].yc)[0]) + (srecent.z -  (halodat[ind_tcurr].zc)[0])*(srecent.z -  (halodat[ind_tcurr].zc)[0]))*(a[ind_tcurr])[0] LT rmax)
slrecent = sl[where(sl.timeform GT tcurr - 0.1*1e9/units.timeunit AND sl.timeform LT tcurr)]
match,slrecent.iorderstar,siord,ind1,ind2
slrecent = slrecent[ind1]
srecent = s[ind2]

srecent_r = sqrt(srecent.x*srecent.x + srecent.y*srecent.y + srecent.z*srecent.z)
slrecent_r = sqrt(slrecent.x*slrecent.x + slrecent.y*slrecent.y + slrecent.z*slrecent.z)
srecent_v = sqrt(srecent.vx*srecent.vx + srecent.vy*srecent.vy + srecent.vz*srecent.vz)
slrecent_v = sqrt(slrecent.vx*slrecent.vx + slrecent.vy*slrecent.vy + slrecent.vz*slrecent.vz)
window,0
plot,srecent.tform,srecent_v*h.time*h.time/slrecent_v,psym = 3;,yrange = [-5,5]
oplot,[0,1],[1,1],color = 100
stop
plot,srecent.tform,srecent_r/slrecent_r,psym=3;,yrange = [-5,5]
oplot,[0,1],[1,1],color = 100
stop

IF 0 THEN BEGIN
loadct,39
;window,0
;plot,sl.timeform*units.timeunit/1e9,sl.x*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_x [comoving]',title = filename + ' ' + haloid;,xrange = [3,3.5];0.5,14];,yrange = [-1e4,1e4]
;oplot,halodat.time/1e9*units.timeunit,halodat.xc*units.lengthunit,psym= 4,color = 254
;oplot,timeform_uniq*units.timeunit/1e9,xcfit*units.lengthunit,color = 254
;oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
;oplot,srecent.tform*units.timeunit/1e9,srecent.x*units.lengthunit,psym = 3,color = 100
;oplot,sl.timeform*units.timeunit/1e9,sl.x*units.lengthunit,psym = 3

window,2
plot,sl.timeform*units.timeunit/1e9,slx*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_x [comoving]',title = filename + ' ' + haloid,xrange = [1,14],xstyle = 1
oplot,halodat.time/1e9*units.timeunit,0*halodat.time,psym = 4,color = 254
oplot,[0,14],[0,0],color = 254
oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
oplot,srecent.tform*units.timeunit/1e9,(srecent.x - (halodat[ind_tcurr].xc)[0])*units.lengthunit,psym = 3,color = 100
stop

;window,0
;plot,sl.timeform*units.timeunit/1e9,sl.y*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_y [comoving]',title = filename + ' ' + haloid;,xrange = [3,3.5];0.5,14];,yrange = [-1e4,1e4]
;oplot,halodat.time/1e9*units.timeunit,halodat.yc*units.lengthunit,psym= 4,color = 254
;oplot,timeform_uniq*units.timeunit/1e9,ycfit*units.lengthunit,color = 254
;oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
;oplot,srecent.tform*units.timeunit/1e9,srecent.y*units.lengthunit,psym = 3,color = 100
;oplot,sl.timeform*units.timeunit/1e9,sl.y*units.lengthunit,psym = 3

window,2
plot,sl.timeform*units.timeunit/1e9,sly*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_y [comoving]',title = filename + ' ' + haloid,xrange = [1,14],xstyle = 1
oplot,halodat.time/1e9*units.timeunit,0*halodat.time,psym = 4,color = 254
oplot,[0,14],[0,0],color = 254
oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
oplot,srecent.tform*units.timeunit/1e9,(srecent.y - (halodat[ind_tcurr].yc)[0])*units.lengthunit,psym = 3,color = 100
stop

;window,0
;plot,sl.timeform*units.timeunit/1e9,sl.z*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_z [comoving]',title = filename + ' ' + haloid;,xrange = [3,3.5];0.5,14];,yrange = [-1e4,1e4]
;oplot,halodat.time/1e9*units.timeunit,halodat.zc*units.lengthunit,psym= 4,color = 254
;oplot,timeform_uniq*units.timeunit/1e9,zcfit*units.lengthunit,color = 254
;oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
;oplot,srecent.tform*units.timeunit/1e9,srecent.z*units.lengthunit,psym = 3,color = 100
;oplot,sl.timeform*units.timeunit/1e9,sl.z*units.lengthunit,psym = 3

window,2
plot,sl.timeform*units.timeunit/1e9,slz*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_z [comoving]',title = filename + ' ' + haloid,xrange = [1,14],xstyle = 1
oplot,halodat.time/1e9*units.timeunit,0*halodat.time,psym = 4,color = 254
oplot,[0,14],[0,0],color = 254
oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
oplot,srecent.tform*units.timeunit/1e9,(srecent.z - (halodat[ind_tcurr].zc)[0])*units.lengthunit,psym = 3,color = 100
stop
ENDIF
;------------------------------------------ Velocity --------------------------------
halodat.vxc = halodat.vxc;a;a ;comoving
halodat.vyc = halodat.vyc;a;a ;comoving
halodat.vzc = halodat.vzc;a;a ;comoving
vxcfit = spline(halodat.time,halodat.vxc,timeform_uniq)
vycfit = spline(halodat.time,halodat.vyc,timeform_uniq)
vzcfit = spline(halodat.time,halodat.vzc,timeform_uniq)
afit = spline(halodat.time,a,timeform_uniq) 
match2,sl.timeform,timeform_uniq,ind1,ind2

slvx = sl.vx/afit[ind1] - vxcfit[ind1]*afit[ind1]
slvy = sl.vy/afit[ind1] - vycfit[ind1]*afit[ind1]
slvz = sl.vz/afit[ind1] - vzcfit[ind1]*afit[ind1]

srecent_r = sqrt((srecent.x - (halodat[ind_tcurr].xc)[0])*(srecent.x - (halodat[ind_tcurr].xc)[0]) + $
                 (srecent.y - (halodat[ind_tcurr].yc)[0])*(srecent.y - (halodat[ind_tcurr].yc)[0]) + $
                 (srecent.z - (halodat[ind_tcurr].zc)[0])*(srecent.z - (halodat[ind_tcurr].zc)[0]))
srecent_v = sqrt((srecent.vx - (halodat[ind_tcurr].vxc)[0])*(srecent.vx - (halodat[ind_tcurr].vxc)[0]) + $
                 (srecent.vy - (halodat[ind_tcurr].vyc)[0])*(srecent.vy - (halodat[ind_tcurr].vyc)[0]) + $
                 (srecent.vz - (halodat[ind_tcurr].vzc)[0])*(srecent.vz - (halodat[ind_tcurr].vzc)[0]))
sr = sqrt((s.x - (halodat[ind_tcurr].xc)[0])*(s.x - (halodat[ind_tcurr].xc)[0]) + $
          (s.y - (halodat[ind_tcurr].yc)[0])*(s.y - (halodat[ind_tcurr].yc)[0]) + $
          (s.z - (halodat[ind_tcurr].zc)[0])*(s.z - (halodat[ind_tcurr].zc)[0]))
sv = sqrt((s.vx - (halodat[ind_tcurr].vxc)[0])*(s.vx - (halodat[ind_tcurr].vxc)[0]) + $
          (s.vy - (halodat[ind_tcurr].vyc)[0])*(s.vy - (halodat[ind_tcurr].vyc)[0]) + $
          (s.vz - (halodat[ind_tcurr].vzc)[0])*(s.vz - (halodat[ind_tcurr].vzc)[0]))
dr = sqrt((d.x - (halodat[ind_tcurr].xc)[0])*(d.x - (halodat[ind_tcurr].xc)[0]) + $
          (d.y - (halodat[ind_tcurr].yc)[0])*(d.y - (halodat[ind_tcurr].yc)[0]) + $
          (d.z - (halodat[ind_tcurr].zc)[0])*(d.z - (halodat[ind_tcurr].zc)[0]))
dv = sqrt((d.vx - (halodat[ind_tcurr].vxc)[0])*(d.vx - (halodat[ind_tcurr].vxc)[0]) + $
          (d.vy - (halodat[ind_tcurr].vyc)[0])*(d.vy - (halodat[ind_tcurr].vyc)[0]) + $
          (d.vz - (halodat[ind_tcurr].vzc)[0])*(d.vz - (halodat[ind_tcurr].vzc)[0]))


window,0
plot,sl.timeform*units.timeunit/1e9,sl.vx*units.vunit/afit[ind1],psym = 3,xtitle = 'Time [Gyr]',ytitle = 'v_x',title = filename + ' ' + haloid,yrange = [-500,500];,xrange = [3,3.5];[tcurr - (0.05*1e9/units.timeunit),tcurr]/1e9*units.timeunit
oplot,halodat.time/1e9*units.timeunit,halodat.vxc*units.vunit*a,psym= 4,color = 254
oplot,sl.timeform/1e9*units.timeunit,vxcfit[ind1]*units.vunit*afit[ind1],color = 254
oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
;oplot,srecent.tform*units.timeunit/1e9,srecent.vx*units.vunit*h.time*h.time,psym = 1,color = 100
;oplot,sl.timeform*units.timeunit/1e9,sl.vx*units.vunit,psym = 3

window,2
plot,sl.timeform*units.timeunit/1e9,slvx*units.vunit/afit[ind1],psym = 3,xtitle = 'Time [Gyr]',ytitle = 'v_x',title = filename + ' ' + haloid ,xrange = [3,3.5];[tcurr - (0.05*1e9/units.timeunit),tcurr]/1e9*units.timeunit
oplot,[0,14],[0,0],color = 254
oplot,halodat.time/1e9*units.timeunit,0*halodat.time,psym = 4,color = 254
oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
oplot,srecent.tform*units.timeunit/1e9,(srecent.vx - (halodat[ind_tcurr].vxc)[0])*units.vunit*h.time,psym = 3,color = 100
 stop

window,0
plot,sl.timeform*units.timeunit/1e9,sl.vy*units.vunit/afit[ind1],psym = 3,xtitle = 'Time [Gyr]',ytitle = 'v_y',title = filename + ' ' + haloid;,xrange = [3,3.5];[tcurr - (0.05*1e9/units.timeunit),tcurr]/1e9*units.timeunit
oplot,halodat.time/1e9*units.timeunit,halodat.vyc*units.vunit*a,psym= 4,color = 254
oplot,sl.timeform/1e9*units.timeunit,vycfit[ind1]*units.vunit*afit[ind1],color = 254
oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
;oplot,srecent.tform*units.timeunit/1e9,srecent.vy*units.vunit*h.time*h.time,psym = 1,color = 100;*(a[ind_tcurr])[0],psym = 1,color = 100
;oplot,sl.timeform*units.timeunit/1e9,sl.vy*units.vunit,psym = 3

window,2
plot,sl.timeform*units.timeunit/1e9,slvy*units.vunit/afit[ind1],psym = 3,xtitle = 'Time [Gyr]',ytitle = 'v_y',title = filename + ' ' + haloid;,xrange = [3,3.5]
oplot,[0,14],[0,0],color = 254
oplot,halodat.time/1e9*units.timeunit,0*halodat.time,psym = 4,color = 254
oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
oplot,srecent.tform*units.timeunit/1e9,(srecent.vy - (halodat[ind_tcurr].vyc)[0])*units.vunit*h.time,psym = 1,color = 100
;stop

window,0
plot,sl.timeform*units.timeunit/1e9,sl.vz*units.vunit/afit[ind1],psym = 3,xtitle = 'Time [Gyr]',ytitle = 'v_z',title = filename + ' ' + haloid,yrange = [min([sl.vz*units.vunit,halodat.vzc*units.vunit]),max([sl.vz*units.vunit,halodat.vzc*units.vunit])];[tcurr - (0.05*1e9/units.timeunit),tcurr]/1e9*units.timeunit
oplot,halodat.time/1e9*units.timeunit,halodat.vzc*units.vunit*a,psym= 4,color = 254
oplot,sl.timeform/1e9*units.timeunit,vzcfit[ind1]*units.vunit*afit[ind1],color = 254
oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
;oplot,srecent.tform*units.timeunit/1e9,srecent.vz*units.vunit*h.time*h.time,psym = 1,color = 100;*(a[ind_tcurr])[0],psym = 1,color = 100
;oplot,sl.timeform*units.timeunit/1e9,sl.vz*units.vunit,psym = 3

window,2
plot,sl.timeform*units.timeunit/1e9,slvz*units.vunit/afit[ind1],psym = 3,xtitle = 'Time [Gyr]',ytitle = 'v_z',title = filename + ' ' + haloid;,xrange = [3,3.5]
oplot,[0,14],[0,0],color = 254
oplot,halodat.time/1e9*units.timeunit,0*halodat.time,psym = 4,color = 254
oplot,[tcurr*units.timeunit/1e9,tcurr*units.timeunit/1e9],[-1e8,1e8]
oplot,srecent.tform*units.timeunit/1e9,(srecent.vz - (halodat[ind_tcurr].vzc)[0])*units.vunit*h.time,psym = 1,color = 100
;stop

;------------------------------ Change of Axis ------------------------
timeform_uniq = sl[uniq(sl.timeform,sort(sl.timeform))].timeform
pos = [[slx],[sly],[slz]]
vel = [[slvx],[slvy],[slvz]]
axfit = spline(halodat.time,halodat.xa,timeform_uniq)
ayfit = spline(halodat.time,halodat.ya,timeform_uniq)
azfit = spline(halodat.time,halodat.za,timeform_uniq)
match2,sl.timeform,timeform_uniq,ind1,ind2
;axfit = axfit[ind1]
;ayfit = ayfit[ind1]
;azfit = azfit[ind1]
pos_prime = pos
vel_prime = vel
FOR i = 0, n_elements(timeform_uniq) - 1 DO BEGIN
   az_prime = [axfit[i],ayfit[i],azfit[i]]
   mag = sqrt(total(az_prime*az_prime))
   az_prime = az_prime/mag
   ax_prime = [az_prime[2]/sqrt(az_prime[0]*az_prime[0] + az_prime[2]*az_prime[2]),0,-1.0*az_prime[0]/sqrt(az_prime[0]*az_prime[0] + az_prime[2]*az_prime[2])]
   ay_prime = crossp(az_prime,ax_prime) 
   basis = [[ax_prime],[ay_prime],[az_prime]]
   ind_time = where(timeform_uniq[i] EQ sl.timeform)
   pos_prime[ind_time,*] = pos[ind_time,*]#basis
   vel_prime[ind_time,*] = vel[ind_time,*]#basis
;   stop

;   az_prime = transpose([[axfit],[ayfit],[azfit]])
;   mag = sqrt((az_prime*az_prime)[0,*] + (az_prime*az_prime)[1,*] + (az_prime*az_prime)[2,*])
;   mag = transpose([[reform(mag)],[reform(mag)],[reform(mag)]])
;   az_prime = az_prime/mag
;   ax_prime = transpose([[reform(az_prime[2,*]/sqrt(az_prime[0,*]*az_prime[0,*] + az_prime[2,*]*az_prime[2,*]))],$
;               [fltarr(n_elements(axfit))],$
;               [reform(-1.0*az_prime[0,*]/sqrt(az_prime[0,*]*az_prime[0,*] + az_prime[2,*]*az_prime[2,*]))]])
;   ay_prime = transpose([[reform(-ax_prime[2,*]*az_prime[1,*])],[reform(ax_prime[2,*]*az_prime[0,*] - ax_prime[0,*]*az_prime[2,*])],[reform(ax_prime[0,*]*az_prime[1,*])]])
ENDFOR
plot,sl.timeform*units.timeunit,vel_prime[*,2]*units.vunit,psym = 3
stop


;IDL> sl2.x = pos_prime[*,0]
;IDL> sl2.y = pos_prime[*,1]
;IDL> sl2.z = pos_prime[*,2]
;IDL> sl2.vx = vel_prime[*,0]
;IDL> sl2.vy = vel_prime[*,1]
;IDL> sl2.vz = vel_prime[*,2]
END
