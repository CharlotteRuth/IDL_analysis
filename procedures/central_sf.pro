;7/1/2015
;Charlotte Christensen
;We see that there is a pile up of gas inside the smoothing length.
;Could that gas be forming most of the central stars?
;Could this be why our bulges are too massive and overly concentrated?
;Here, we investigate


PRO central_sf,galaxy,finalstep,halo

units = tipsyunits(galaxy + '.param')

IF file_test(galaxy + '_2merge.starlog.fits') THEN BEGIN
   print,galaxy + '_2merge.starlog.fits'
   sl = mrdfits(galaxy + '_2merge.starlog.fits',1) 
ENDIF ELSE BEGIN
   IF file_test(galaxy + '.starlog.fits') THEN BEGIN
      print,galaxy + '.starlog.fits'
      sl = mrdfits(galaxy + '.starlog.fits',1)        
   ENDIF ELSE BEGIN
      print,'Starlog: ',galaxy + '.starlog'
      sl = rstarlog(galaxy + '.starlog',/molecularH)
   ENDELSE
ENDELSE
;sl = rstarlog(galaxy + '.starlog',/molecularH)
stop

sl.timeform = sl.timeform*units.timeunit/1e9
sl.x = sl.x*units.lengthunit
sl.y = sl.y*units.lengthunit
sl.z = sl.z*units.lengthunit
sl.vx = sl.vx*units.vunit
sl.vy = sl.vy*units.vunit
sl.vz = sl.vz*units.vunit

rtipsy,'steps/' + galaxy + '.' + finalstep + '.dir/' + galaxy + '.' + finalstep + '.halo.' + halo,h,g,d,s
readarr,'steps/' + galaxy + '.' + finalstep + '.dir/' + galaxy + '.' + finalstep + '.halo.' + halo + '.iord',h,iord,type = 'long',/ascii
siord = iord[h.ngas + h.ndark:h.n - 1]
s.tform = s.tform*units.timeunit/1e9
s.x = s.x*units.lengthunit
s.y = s.y*units.lengthunit
s.z = s.z*units.lengthunit
s.vx = s.vx*units.vunit
s.vy = s.vy*units.vunit
s.vz = s.vz*units.vunit

match2,siord,sl.iorderstar,ind_t,ind_sl
slhalo = sl[where(ind_sl NE -1)]
siord = siord[where(ind_t NE -1)]
s = s[where(ind_t NE -1)]

;
plot,slhalo.timeform,slhalo.x,psym = 3
oplot,s.tform,s.x,psym = 3,color = 100
;stop

;siord_all = iordall[hall.ngas + hall.ndark:hall.n - 1]
;match,siord_all,sl.iorderstar,temp,ind
;sall = sall[temp]
;siord_all = siord_all[temp]
;slcosmo = sl[ind]
;match,siord_all,siord,haloind,temp

align = mrdfits('grp' + halo + '.alignment.fits',1)
timeform_uniq = slhalo[uniq(slhalo.timeform)].timeform
xc_extrap_clip = spline(align.time,align.xc*units.lengthunit,timeform_uniq)
yc_extrap_clip = spline(align.time,align.yc*units.lengthunit,timeform_uniq)
zc_extrap_clip = spline(align.time,align.zc*units.lengthunit,timeform_uniq)
oplot,timeform_uniq,xc_extrap_clip,psym = 3, color = 200

match2,slhalo.timeform,timeform_uniq,ind,temp

xc_extrap = xc_extrap_clip[ind]
yc_extrap = yc_extrap_clip[ind]
zc_extrap = zc_extrap_clip[ind]

slhalo.x = slhalo.x - xc_extrap
slhalo.y = slhalo.y - yc_extrap
slhalo.z = slhalo.z - zc_extrap

;plot,slhalo.timeform,slhalo.x,psym = 3
;oplot,s.tform,s.x,psym = 3,color = 100
;stop
;slhalo.x = slhalo.x

eps = max(s.eps*units.lengthunit)
histogramp,slhalo.x,nbins = 100,min = -10,max = 10,/normalize
histogramp,slhalo.y,nbins = 100,min = -10,max = 10,/overplot,color = 100,/normalize
histogramp,slhalo.z,nbins = 100,min = -10,max = 10,/overplot,color = 60,/normalize
oplot,[-1*eps,-1*eps],[0,1],linestyle = 2
oplot,[eps,eps],[0,1],linestyle = 2
stop
END


PRO central_sf_master,outplot = outplot,thick = thick
formatplot,outplot = outplot,thick = thick

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' $
ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN prefix = '/nobackupp8/crchrist/MolecH/' $
ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
file986 = 'h986.cosmo50cmb.3072g14HBWK'
key986 = 'h986'
dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
file603 = 'h603.cosmo50cmb.3072g14HBWK'
key603 = 'h603'
dir277 = prefix + 'h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/'
file277 = 'h277.cosmo50cmb.3072g14HMbwK'
key277 = 'h277'
dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
file258 = 'h258.cosmo50cmb.3072g14HMbwK'
key258 = 'h258'
dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
file285 = 'h285.cosmo50cmb.3072g14HMbwK'
key285 = 'h285'
dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
file239 = 'h239.cosmo50cmb.3072g14HMbwK'
key239 = 'h239'
;dirs = [dir986,dir603,dir277,dir258,dir285,dir239]
;files = [file986,file603,file277,file258,file285,file239]
dirs = [dir986,dir258,dir285,dir239]
files = [file986,file258,file285,file239]
finalstep = '00512'
outfiles = dirs + files + '.' + finalstep + '/' + files + '.' + finalstep
haloid = ['1','1','1','1','1']
outhalos = outfiles + '.halo.' + haloid 
key = [key986,key258,key285,key239] + ', ' + haloid

FOR i = 0, n_elements(files) - 1 DO BEGIN
    cd,dirs[i]
    central_sf,files[i],finalstep,haloid[i]
ENDFOR

END
