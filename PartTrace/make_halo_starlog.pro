;CC 8/8/12
;If a starlog file is not avaliable, this program can make one sort of
;like it for a given halo

;It also tells you where a star was at the last output, as opposed to
;when formed (as the starlog file does)

;filebase = 'h603.cosmo50cmb.3072g14HBWK'
;make_halo_starlog,filebase,/molecularH

PRO make_halo_starlog,filebase,finalid = finalid,molecularH = molecularH
IF NOT keyword_set(finalid) THEN finalid = '1'
units = tipsyunits(filebase + '.param')
halodat = mrdfits('grp' + finalid + '.alignment.fits',1)
files = halodat.file ;going forward in time
haloid = strtrim(string(fix(halodat.haloid)),2)
halodat.xc = (halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0)*1000.0 - units.lengthunit/2.0 
halodat.yc = (halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0)*1000.0 - units.lengthunit/2.0 
halodat.zc = (halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0)*1000.0 - units.lengthunit/2.0
halodat.rvir = halodat.rvir;*units.lengthunit

ageUniverse = 13.7346*1e9

starttime = 0
iord_arr     = [0]
timeform_arr = [0]
massform_arr = [0]
x_arr = [0]
y_arr = [0]
z_arr = [0]
n = n_elements(files)

rtipsy,files[n-1],h,g,d,s;,/justhead
;readarr,files[n-1] + '.iord',h,iord_orig,part = 'star',/ascii, type = 'long'
;readarr,files[n-1] + '.massform',h,massform_orig,part = 'star',/ascii, type = 'float'
;readarr,files[n-1] + '.timeform',h,timeform_orig,part = 'star',/ascii, type = 'float'
read_tipsy_arr,files[n-1] + '.iord',h,iord_orig,part = 'star',type = 'long'
;IF file_test(filebase + 'starlog') THEN BEGIN
;   sldata = rstarlog(filebase + 'starlog',molecularH = molecularH)
;   massform_orig = sldata.massform
;   timeform_orig = sldata.timeform
;   stop
;ENDIF ELSE BEGIN
   read_tipsy_arr,files[n-1] + '.massform',h,massform_orig,part = 'star',type = 'float'
;   read_tipsy_arr,files[n-1] + '.timeform',h,timeform_orig,part =   'star',type = 'float'
   timeform_orig = s.tform
;ENDELSE 

FOR i = 0, n - 1 DO BEGIN
;    rtipsy,files[i],h_step,g_step,d_step,s_step,/justhead
;    scale = h_step.time
   scale = 1./(1. + halodat[i].z)
;    readarr,files[i]+ '.halo.' + haloid[i] + '.iord',h_halo,iord_halo,part = 'star',/ascii, type = 'long'
    data_orig = mrdfits(files[i] + '.grp' + finalid + '.star.history.fits',1)
;    data = data_orig[where(data_orig.haloid EQ haloid[i] AND data_orig.rho EQ 0)]
    data = data_orig[where((sqrt((data_orig.x - halodat[i].xc*scale)^2 + (data_orig.y - halodat[i].yc*scale)^2 + (data_orig.z - halodat[i].zc*scale)^2) LE halodat[i].rvir) $
                           AND data_orig.rho EQ 0)]
    starnow = data_orig[where(data_orig.rho EQ 0)]

    indrecent = where(timeform_orig gt starttime)
    iord     = iord_orig[indrecent]
    massform = massform_orig[indrecent]
    timeform = timeform_orig[indrecent]
   
    match,iord,data.iord,ind,ind_halo
    iord     = iord[ind]
    massform = massform[ind]
    timeform = timeform[ind]
    gpos = [[data[ind_halo].x - halodat[i].xc*scale],[data[ind_halo].y - halodat[i].yc*scale],[data[ind_halo].z - halodat[i].zc*scale]]
 
    az = [halodat[i].xa,halodat[i].ya,halodat[i].za]
    az0 = az[0]
    az1 = az[1]
    az2 = az[2]
    ax = [az2/sqrt(az0*az0 + az2*az2),0,-1.0*az0/sqrt(az0*az0 + az2*az2)]
    ay = crossp(az,ax)          ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
    basis = [[ax],[ay],[az]]
    gpos = transpose(transpose(basis)#transpose(gpos))

    iord_arr     = [iord_arr,iord]
    massform_arr = [massform_arr,massform]   
    timeform_arr = [timeform_arr,timeform]
    x_arr = [x_arr,reform(gpos[*,0])]
    y_arr = [y_arr,reform(gpos[*,1])]
    z_arr = [z_arr,reform(gpos[*,2])]

    match,iord_orig,starnow.iord,ind,ind_now
    starttime = max(timeform_orig[ind])
    print,starttime,n_elements(data),n_elements(ind),n_elements(iord_arr)
;    histogramp,x_arr,min = min(x_arr),max = max(x_arr)
;    histogramp,reform(gpos[*,0]),/overplot,color = 100
;    stop
ENDFOR
nstars = n_elements(iord_arr)
iord_arr     = iord_arr[1:nstars - 1]
massform_arr = massform_arr[1:nstars - 1]
timeform_arr = timeform_arr[1:nstars - 1]
x_arr = x_arr[1:nstars - 1]
y_arr = y_arr[1:nstars - 1]
z_arr = z_arr[1:nstars - 1]

;window,0
;histogramp,x_arr,min = -20,max = 20, nbins = 200
;histogramp,y_arr,/overplot,linestyle = 2
;histogramp,z_arr,/overplot,linestyle = 1

;window,2
;plot,x_arr,y_arr,psym = 3,yrange = [-5,5],xrange = [-5,5]
IF keyword_set(molecularH) THEN starlog = {iorderstar:0L,iordergas:0L,timeform:0d,x:0d,y:0d,z:0d,vx:0d,vy:0d,vz:0d,massform:0d,rhoform:0d,tempform:0d,H2form:0d} ELSE starlog = {iorderstar:0,iordergas:0,timeform:0d,x:0d,y:0d,z:0d,vx:0d,vy:0d,vz:0d,massform:0d,rhoform:0d,tempform:0d}
starlogcut = replicate(starlog, nstars - 1)
starlogcut.iorderstar = iord_arr
starlogcut.massform   = massform_arr
starlogcut.timeform   = timeform_arr 
starlogcut.x = x_arr
starlogcut.y = y_arr
starlogcut.z = z_arr

mwrfits,starlogcut ,'starlog.' + finalid + '.fits',/create
END
