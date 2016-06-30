;CC, 12/12/11

;This program determins the characteristics (density, distance from
;center etc.) of the gas at the step
;prior to outflow

;data = fb_gas_character('h516.cosmo25cmb.3072g1MBWK',25000)
FUNCTION fb_gas_character,base,kpcunit,totalm = totalm
;loadct,39
spawn,'ls *512/*512.coolontime',file_coolon
spawn,'ls *512/*512.iord',file_iord
spawn,'ls *512/*512',file
rtipsy,file,h,g,d,s,/justhead
readarr,file_coolon,h,coolon,part = 'gas',/ascii
readarr,file_iord,h,iord,part = 'gas',/ascii

outflow_z    = mrdfits('grp1.rvir.outflow_z.fits')
outflow_iord = mrdfits('grp1.rvir.outflow_iord.fits')
outflow_mass = mrdfits('grp1.rvir.mass_at_outflow.fits')
outflow_disk = mrdfits('grp1.rvir.outflow_timeindisk.fits')

;discard all gas particles that are never expelled (z = 99)
trash = where(outflow_z eq 99.0,complement = keep)
outflow_z = outflow_z[keep]
outflow_iord = outflow_iord[keep]
outflow_mass = outflow_mass[keep]
outflow_disk = outflow_disk[keep]
match,iord,outflow_iord,ind_iord,ind_iordout
coolonout = where(coolon[ind_iord] ne 0)

outflow_z = outflow_z[ind_iordout[coolonout]]
outflow_iord = outflow_iord[ind_iordout[coolonout]]
outflow_mass = outflow_mass[ind_iordout[coolonout]]
outflow_disk = outflow_disk[ind_iordout[coolonout]]
IF 0 THEN BEGIN
    readcol,base+'.haloid.dat',files,haloid,format= '(A,F)'
    files = REVERSE(files)      ;going forward in time
    haloid = REVERSE(haloid)
    center = read_center(files = files + '.amiga.stat',halos = haloid)
    center[*,2] = center[*,2]*1000.0 - kpcunit/2.0 ;"center" is in Mpc for a box that goes from [0,0,0] to [kpcunit/1000,kpcunit/1000,kpcunit/1000]
    center[*,3] = center[*,3]*1000.0 - kpcunit/2.0
    center[*,4] = center[*,4]*1000.0 - kpcunit/2.0
ENDIF ELSE BEGIN
    readcol,'alignment.txt',files,haloid,time,z,xc,yc,zc,xa,ya,za,format = '(A,I,F,F,F,F,F,F,F,F)'
    center = [[time],[z],[xc],[yc],[zc]]
    center[*,2] = center[*,2]*kpcunit ;"center" is in Mpc for a box that goes from [0,0,0] to [kpcunit/1000,kpcunit/1000,kpcunit/1000]
    center[*,3] = center[*,3]*kpcunit
    center[*,4] = center[*,4]*kpcunit
    a = [[xa],[ya],[za]]
ENDELSE

FOR i = 0, N_ELEMENTS(files) - 2 DO BEGIN
    z = center[i,1]
    z_fut = center[i + 1,1]
    scale = 1.0/(1.0 + z)
    fb_set = where(outflow_disk le z AND outflow_disk gt z_fut)
    IF fb_set[0] NE -1 THEN BEGIN
        fb_set_iord = outflow_iord[fb_set]
        fb_set_mass = outflow_mass[fb_set]
        fb_set_z = outflow_z[fb_set]
        fb_set_disk = outflow_disk[fb_set]

        dir = (strsplit(files[i],'/',/extract))[0]
        gashistory = mrdfits(dir + '/' + dir + '.allgas.history.fits',1)
        
        match,fb_set_iord,gashistory.iord,fb_set_ind,gashistory_ind
        temp2 = gashistory[gashistory_ind]
        temp2z = fb_set_z[fb_set_ind]
        temp2m = fb_set_mass[fb_set_ind]
        IF (where(temp2.temp le 2e4))[0] ne -1 THEN BEGIN
            temp = temp2[where(temp2.temp le 2e4)]
            tempz = temp2z[where(temp2.temp le 2e4)]
            tempm = temp2m[where(temp2.temp le 2e4)]
            temp2.x = (temp2.x - center[i,2])*scale
            temp2.y = (temp2.y - center[i,3])*scale
            temp2.z = (temp2.z - center[i,4])*scale
            gpos2 = [[temp2.x],[temp2.y],[temp2.z]]
            temp.x = (temp.x - center[i,2])*scale
            temp.y = (temp.y - center[i,3])*scale
            temp.z = (temp.z - center[i,4])*scale
            gpos = [[temp.x],[temp.y],[temp.z]]
            
            az = reform(a[i,*])
            az1 = a[i,0]
            az2 = a[i,1]
            az3 = a[i,2]
            ax = [az3/sqrt(az1*az1 + az3*az3),0,-1.0*az1/sqrt(az1*az1 + az3*az3)]
            ay = crossp(az,ax) ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
            basis = [[ax],[ay],[az]]
            
            gpos2 = transpose(basis#transpose(gpos2))
            temp2.x = gpos2[*,0]
            temp2.y = gpos2[*,1]
            temp2.z = gpos2[*,2]
            gpos = transpose(basis#transpose(gpos))
            temp.x = gpos[*,0]
            temp.y = gpos[*,1]
            temp.z = gpos[*,2]

            IF N_ELEMENTS(fbhistory) eq 0 THEN fbhistory = temp $
            ELSE fbhistory = [fbhistory,temp]
 
;            histogramp,temp2.z,title = files[i]
;            histogramp,temp.z,title = files[i],xrange = [-10,10],nbins = 100
;            histogramp,temp.x,/overplot,color = 240,nbins = 100
;            histogramp,temp.y,/overplot,color = 50, nbins = 100
;            histogramp,temp2.z,/overplot,linestyle = 2
;            histogramp,temp2.x,/overplot,color = 240,linestyle = 2
;            histogramp,temp2.y,/overplot,color = 50,linestyle = 2

             IF N_ELEMENTS(total) eq 0 THEN BEGIN
                total = temp
                totalz = tempz
                totalm = tempm
            ENDIF ELSE BEGIN
                total = [total,temp]
                totalz = [totalz,tempz]
                totalm = [totalm,tempm]
            ENDELSE
;            stop
        ENDIF
    ENDIF
ENDFOR
;histogramp,total.z,xrange = [-10,10],linestyle = 2
;histogramp,total.y,color = 240,/overplot,linestyle = 2
;histogramp,total.x,color = 50,/overplot,linestyle = 2
;histogramp,SQRT(total.x*total.x + total.y*total.y),/overplot
;stop

;histogramp,abs(total.z),xrange = [0,10],linestyle = 2
;histogramp,abs(total.y),color = 240,/overplot,linestyle = 2
;histogramp,abs(total.x),color = 50,/overplot,linestyle = 2
;histogramp,SQRT(total.x*total.x + total.y*total.y),/overplot
;stop

;histogramp,SQRT(total.x*total.x + total.y*total.y),xrange = [0,10]
return,total
END

PRO fb_gas_character_master,outplot = outplot
files = ['h516.cosmo25cmb.3072g1MBWK','h516.cosmo25cmb.3072g14HBWK']
dir = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK','/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK']
kpcunit = [25000,25000]
colors = [50,240]

files = ['h986.cosmo50cmb.3072g1bwK','h986.cosmo50cmb.3072gs1MbwK','h986.cosmo50cmb.3072g14HBWK']
dir = ['/astro/store/nbody3/fabio/h986/3072g1bwK',$
       '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK',$
       '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK']
kpcunit = [50000,50000,50000]
colors = [50,120,240]
thicks = [6,6,6]

files = ['h603.cosmo50cmb.3072g1bwK','h603.cosmo50cmb.3072gs1MbwK','h603.cosmo50cmb.3072g14HBWK']
dir = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.3072g1bwK',$
       '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK',$
       '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK']
kpcunit = [50000,50000,50000]
colors = [50,120,240]
thicks = [6,6,6]

formatplot,outplot=outplot,/thick
IF KEYWORD_SET(outplot) THEN device,filename = outplot,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch,/color

FOR i = 0, N_ELEMENTS(files)-1 DO BEGIN
    cd,dir[i]
    fb_gas = fb_gas_character(files[i],kpcunit[i])
    IF i EQ 0 THEN histogramp,SQRT(fb_gas.x*fb_gas.x + fb_gas.y*fb_gas.y),xrange = [0,30],/normalize,yrange = [0,0.1],xtitle = 'Radius [kpc]',ytitle = '1/N dN/dr'
    histogramp,SQRT(fb_gas.x*fb_gas.x + fb_gas.y*fb_gas.y),/overplot,/normalize,color = colors[i],thick = thicks[i]
ENDFOR

IF KEYWORD_SET(outplot) THEN device,/close

END
