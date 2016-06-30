;tfile = 'h603.cosmo50cmb.3072g14HBWK.00444'
;dist_units = 50000
;outarray = ['HI','H2','iord','smoothlength','lw']


;cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00372.dir'
;tfile = 'h986.cosmo50cmb.3072g14HBWK.00372'

PRO mark2tipsy,tfile,dist_units,outarray = outarray
;tfile = 'h603.cosmo50cmb.3072g14HBWK.00444'
;dist_units = 50000
;outarray = ['HI','H2','iord','smoothlength','lw']

mfile = tfile + '.mark'
rtipsy,tfile,h,g,d,s
openr,1,mfile
readf,1,n,ngas,nstar
ndark = n - ngas - nstar
close,1
readcol,mfile,iord,format='(d)'
iord = LONG(iord) - 1 ;Mark file iorders are +1
ngashalo = 1

gas_iord = iord[where(iord LT h.ngas)]
;dark_iord = iord[where(iord LT h.ndark + h.ngas AND iord GE h.ngas)] - h.ngas
;star_iord = iord[where(iord LT h.ndark + h.ngas + nstar AND iord GE h.ngas + h.ndark)] - (h.ngas + h.ndark)

gas = g[gas_iord]
gx = MEAN(gas.x)
gy = MEAN(gas.y)
gz = MEAN(gas.z)
gas.x = gas.x - gx
gas.y = gas.y - gy
gas.z = gas.z - gz
radii = max(sqrt(gas.x*gas.x + gas.y*gas.y + gas.z*gas.z))

s.x =  s.x - gx
s.y =  s.y - gy
s.z =  s.z - gz
rads = sqrt(s.x*s.x + s.y*s.y + s.z*s.z)
stars_iord = where(rads le radii)
stars = s[stars_iord]

d.x = d.x - gx
d.y = d.y - gy
d.z = d.z - gz
radd = sqrt(d.x*d.x + d.y*d.y + d.z*d.z)
dark_iord = where(radd le radii)
dark = d[dark_iord]

;loadct,39
;plot,dark.x,dark.y,psym = 3
;oplot,gas.x,gas.y,psym = 3,color = 240
;oplot,stars.x,stars.y,psym = 3,color = 50

;dark = d[dark_iord]
;star = s[star_iord]

;************** Center Halo ***************
;foo = min(stars.phi, cmindex)
;cx = stars[cmindex].x
;cy = stars[cmindex].y
;cz = stars[cmindex].z

;stars.x =  stars.x - cx
;stars.y =  stars.y - cy
;stars.z =  stars.z - cz

;dark.x = dark.x - cx
;dark.y = dark.y - cy
;dark.z = dark.z - cz

;gas.x = gas.x - cx
;gas.y = gas.y - cy
;gas.z = gas.z - cz

propx=mean(stars.vx)
propy=mean(stars.vy)
propz=mean(stars.vz)
    
stars.vx=stars.vx-propx
stars.vy=stars.vy-propy
stars.vz=stars.vz-propz
gas.vx=gas.vx-propx
gas.vy=gas.vy-propy
gas.vz=gas.vz-propz
dark.vx=dark.vx-propx
dark.vy=dark.vy-propy
dark.vz=dark.vz-propz


;********************* Align ***********
read_tipsy_arr,tfile+'.HI',h,HI,part = 'gas'
indg = gas_iord
indd = dark_iord + h.ngas
inds = stars_iord + h.ngas + h.ndark 
ind = [indg,indd,inds]
hi = hi[ind]
limit=5*h.time/dist_units ;use stars within this limit (kpc/dist_units)
alignHI,stars,dark,gas,limit;,hi = hi

num = h 
num.nstar = n_elements(stars)
num.ngas = n_elements(gas)
num.ndark = n_elements(dark)
num.n = num.nstar + num.ngas + num.ndark
fileout = tfile + '.halo.1.std'

IF (1) THEN BEGIN
    if keyword_set(ascii) then $
      wtipsy, fileout, num, gas, dark, stars $
    else wtipsy,fileout,num,gas,dark,stars,/standard
ENDIF ELSE wtipsy,fileout,num,gas,dark,stars

indd = indd + h.ngas
inds = inds + h.ngas + h.ndark 
ind = [indg,indd,inds]

if keyword_set(outarray) then begin
    FOR k = 0, N_ELEMENTS(outarray) - 1 DO BEGIN
        IF (outarray[k] eq 'decomp' OR outarray[k] eq 'iord' ) THEN read_tipsy_arr,tfile+'.' +outarray[k],h,array,type = 'long' ELSE read_tipsy_arr,tfile+'.' +outarray[k],h,array,type = 'float'
        fileout = tfile + '.halo.1.' +outarray[k] ;+'.HI'
        openw,lunhi,fileout,/get_lun
        printf,lunhi,num.n
        for j=0L,n_elements(ind)-1 DO $
          IF (outarray[k] eq 'decomp' OR outarray[k] eq 'iord' ) THEN printf,lunhi,array[ind[j]],format = '(I)'  ELSE printf,lunhi,array[ind[j]],format = '(g)'
        close,lunhi
    ENDFOR
endif
IF FILE_TEST(tfile+'.H2') THEN molecularH = 1 ELSE molecularH = 0
IF KEYWORD_SET(writemacro) THEN writeHImacro,tfile + '.halo.' + strtrim(string(i),2),dist_units = dist_units,mass_unit = mass_unit,molecularH = molecularH
END
