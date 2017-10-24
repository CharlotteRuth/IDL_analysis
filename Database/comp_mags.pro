;Compares the magnitude of different mass galaxies found by different methods
;comp_mags,'X2X2'

pro comp_mags, outputfile
;msol = 2.362e5
msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 1.69875e16 ; Solar mass in system units
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 50000.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70
loadct,39
outputbase = '/astro/net/scratch2/christensen/Database/figures'

;*************************************************Reading in Files********************************************************

;list1 = 'halo_list'+outputfile+'.dat'
;readcol, list1, format = 'A60', files1
mag_file_leo = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/h1201-h1755-'+outputfile+'g1bwK.00512.amiga_vir.halos.leo_ab.Mv.fits'
mag_file_star = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/h1201-h1755-'+outputfile+'g1bwK.00512.amiga_vir.halos.star99_ab.Mv.fits'
mag_file_unscat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/512/X2.unscat.mag'
mag_file_scat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/512/X2.scat.mag'
;mag_file_unscat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/512/X3.unscat.mag'
;mag_file_scat ='../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/512/X3.scat.mag'
;mag_file_unscat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/512/X5.unscat.mag'
;mag_file_scat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/512/X5.scat.mag'
;mag_file_leo = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/h1201-h1755-'+outputfile+'g1bwK.00084.amiga_vir.halos.leo_ab.Mv.fits'
;mag_file_star = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/h1201-h1755-'+outputfile+'g1bwK.00084.amiga_vir.halos.star99_ab.Mv.fits'
;mag_file_unscat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/084/X5.unscat.mag'
;mag_file_scat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/084/X5.scat.mag'
;mag_file_leo = '../Fabio_Runs/h1201-h1755-'+outputfile+'bwK/h1201-h1755-'+outputfile+'bwK.00512.amiga.halos.leo_ab.Mv.fits'
;mag_file_star = '../Fabio_Runs/h1201-h1755-'+outputfile+'bwK/h1201-h1755-'+outputfile+'bwK.00512.amiga.halos.star99_ab.Mv.fits'
;mag_file_unscat = '../Fabio_Runs/h1201-h1755-'+outputfile+'bwK/Sunrise/512/X2g3.unscat.mag'
;mag_file_scat = '../Fabio_Runs/h1201-h1755-'+outputfile+'bwK/Sunrise/512/X2g3.scat.mag'

;outputfile = outputfile + '_z3'

;print,mag_file_star
halo_mags_leo = MRDFITS(mag_file_leo,1)
halo_mags_star = MRDFITS(mag_file_star,1)
;badv = [-23.0020,-23.3235,-22.8198,-24.0543,-24.4341,-24.5059,-24.5326,-24.4630,-22.7303,-23.9468]
;badind = where(halo_mags1.v eq badv
;halo_mags1[badind].id = -1
slice = where(halo_mags_leo.id NE -1)
id_idl = STRING(STRTRIM(FIX(halo_mags_leo[slice].id),2))
imag_leo = halo_mags_leo[slice].i
imag_star = halo_mags_star[slice].i
rmag_leo = halo_mags_leo[slice].r
rmag_star = halo_mags_star[slice].r
bmag_leo = halo_mags_leo[slice].b
bmag_star = halo_mags_star[slice].b
kmag_leo = halo_mags_leo[slice].k
kmag_star = halo_mags_star[slice].k
vmag_leo = halo_mags_leo[slice].v
vmag_star = halo_mags_star[slice].v
n_cor_halos = N_ELEMENTS(slice)

readcol,mag_file_unscat,id_unscat,u_unscat,b_unscat,v_unscat,r_unscat,i_unscat,us_unscat,gs_unscat,rs_unscat,is_unscat,zs_unscat,i1_unscat,i2_unscat,m2_unscat,m7_unscat,m1_unscat,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
readcol,mag_file_scat  ,id_scat  ,u_scat,b_scat,v_scat,r_scat,i_scat,us_scat,gs_scat,rs_scat,is_scat,zs_scat,i1_scat,i2_scat,m2_scat,m7_scat,m1_scat,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
;readcol,mag_file_unscat,id_unscat,u_unscat,b_unscat,v_unscat,r_unscat,i_unscat,format = 'A,F,F,F,F,F'
;readcol,mag_file_scat  ,id_scat  ,u_scat,b_scat,v_scat,r_scat,i_scat,format='A,F,F,F,F,F'
id_unscat = id_unscat(where (i_unscat NE 0))
id_scat = id_scat(where (i_scat NE 0))
i_unscat = i_unscat(where (i_unscat NE 0))
i_scat = i_scat(where (i_scat NE 0))
b_unscat = b_unscat(where (i_unscat NE 0))
b_scat = b_scat(where (i_scat NE 0))
r_unscat = r_unscat(where (i_unscat NE 0))
r_scat = r_scat(where (i_scat NE 0))
v_unscat = v_unscat(where (i_unscat NE 0))
v_scat = v_scat(where (i_scat NE 0))


filename = 'Data_files/'+outputfile+'randv.dat'
readcol,filename,id_rad,r2001, v2001, sm2001, bm2001, r11, v11, sm11, bm11, r21, v21, sm21, bm21,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent


;***********************************************Matching Halos******************************************************************
id_cor_SR = strarr(n_cor_halos)
id_cor_Rad = strarr(n_cor_halos)
imag_unscat = fltarr(n_cor_halos) ;unscat
imag_scat = fltarr(n_cor_halos) ;scat
bmag_unscat = fltarr(n_cor_halos)
bmag_scat = fltarr(n_cor_halos)
rmag_unscat = fltarr(n_cor_halos)
rmag_scat = fltarr(n_cor_halos)
vmag_unscat = fltarr(n_cor_halos)
vmag_scat = fltarr(n_cor_halos)
v11_match = fltarr(n_cor_halos)
v21_match = fltarr(n_cor_halos)
v2001_match = fltarr(n_cor_halos)
sm11_match = fltarr(n_cor_halos)
bm11_match = fltarr(n_cor_halos)
bm2001_match = fltarr(n_cor_halos)
matchSR_idl = intarr(n_cor_halos) - 1
matchRad_idl = intarr(n_cor_halos) - 1

FOR ct = 0, n_elements(id_rad)-1 DO BEGIN
    IF (strpos(id_rad[ct],'0') EQ 0) THEN id_rad[ct] = STRMID(id_rad[ct],1,STRLEN(id_rad[ct])-1)
    IF (strpos(id_rad[ct],'0') EQ 0) THEN id_rad[ct] = STRMID(id_rad[ct],1,STRLEN(id_rad[ct])-1)
    IF (strpos(id_rad[ct],'0') EQ 0) THEN id_rad[ct] = STRMID(id_rad[ct],1,STRLEN(id_rad[ct])-1)
ENDFOR
FOR ct = 0, n_elements(id_idl)-1 DO BEGIN
;Loop through the magnitude determined halos to mark which correspond
;to those that are Sunrised and those that have velocities
    index = where(strcmp(id_idl[ct],id_unscat) EQ 1)
    IF (index  NE -1) THEN BEGIN
;If there is a match between idl magnitudes and sunrised magnitudes
        id_cor_SR[ct]= id_unscat[index]
        imag_unscat[ct] = i_unscat[index] ;unscat
        rmag_unscat[ct] = r_unscat[index] ;unscat
        bmag_unscat[ct] = b_unscat[index] ;unscat
        vmag_unscat[ct] = v_unscat[index] ;unscat
        imag_scat[ct] = i_scat[index] ;scat
        rmag_scat[ct] = r_scat[index] ;scat
        bmag_scat[ct] = b_scat[index] ;scat
        vmag_scat[ct] = v_scat[index] ;scat
        matchSR_idl[ct] = index
;        stop
    ENDIF ELSE matchSR_idl[ct] = -1

    index = where(strcmp(id_idl[ct],id_rad) EQ 1)
    IF(index[0] NE -1) THEN matchRad_idl[ct] = index ELSE matchRad_idl[ct] = -1
;If there is a match between veolcities and magnitudes
ENDFOR

matchidl_rad = intarr(N_ELEMENTS(id_rad))
FOR ct = 0, n_elements(id_rad)-1 DO BEGIN
    index = where(strcmp(id_rad[ct],id_idl) EQ 1)
    IF(index[0] NE -1) THEN matchidl_rad[ct] = index ELSE matchidl_rad[ct] = -1
;If there is a match between veolcities and magnitudes
ENDFOR

matchall_idl = Where(matchRad_idl NE -1 AND matchSR_idl NE -1)
matchall_SR = matchSR_idl[matchAll_idl]
matchall_Rad = matchRad_idl[matchAll_idl]

matchSR_idls = where(matchSR_idl NE -1)
matchidl_rads = where(matchidl_rad NE -1)

matchidl_SRs = matchSR_idl[matchSR_idls]
matchrad_idls = matchidl_rad[matchidl_rads]

print,'id_scat[matchall_SR]:  ';,id_scat[matchall_SR]
print,'id_idl[matchall_idl]:  ';,id_idl[matchall_idl]
print,'id_rad[matchall_rad]:  ';,id_rad[matchall_rad]
print,''
print,'id_scat[matchidl_SRs]: ';,id_scat[matchidl_SRs]
print,'id_idl[matchSR_idl]:   ';,id_idl[matchSR_idls]
print,''
print,'id_rad[matchidl_rads]: ';,id_rad[matchidl_rads]
print,'id_idl[matchrad_idls]: ';,id_idl[matchrad_idls]


print,'Number of Magnitude Determinded Halos: ',STRTRIM(N_ELEMENTS(id_idl),2)
print,'Number of Velocity Determinded Halos:  ',STRTRIM(N_ELEMENTS(id_rad),2),  ', (matches: ',STRTRIM(N_ELEMENTS(matchRad_idl),2),')'
print,'Number of Sunrised Halos:              ',STRTRIM(N_ELEMENTS(id_unscat),2),', (matches: ',STRTRIM(N_ELEMENTS(matchSR_idl),2),')'
print,'Number of Halos with complete data:    ',STRTRIM(N_ELEMENTS(matchaLL_SR),2)


;*****************************Write Results****************************************************
close,1
openw,1,'Data_files/' + outputfile+'_data.dat'
printf,1,'ID ','v200 ','v20 ','v10 ','Ileo ','Kleo ','Istar ','Kstar ','Isun '
;printf,1,TRANSPOSE([[halo_mags_leo[where(matchRad eq 1)].id],[v2001[where(matchRad eq 1)]],[v11[where(matchRad eq 1)]],[v21[where(matchRad eq 1)]],[imag1[where(matchRad eq 1)]],[kmag1[where(matchRad eq 1)]],[imag2[where(matchRad eq 1)]],[kmag2[where(matchRad eq 1)]],[imag3[where(matchRad eq 1)]]])
printf,1,TRANSPOSE([[id_idl[matchall_idl]],[v2001[matchall_rad]],[v11[matchall_rad]],[v21[matchall_rad]],[imag_leo[matchall_idl]],[kmag_leo[matchall_idl]],[imag_star[matchall_idl]],[kmag_star[matchall_idl]],[i_unscat[matchall_SR]]])
close,1

;stop

!Y.style = 1
!X.style = 1
;****************************Magnitude Comparison************************************************
set_plot,'x'
plot,[0,-24],[0,-24],xrange = [-23,-14],yrange = [-23,-14],xtitle = 'Starburst99 V mag',ytitle = 'Unscat Sunrise V mag',linestyle = 2
;oplot,imag_star[matchSR_idl],imag_unscat[matchidl_SRs],psym = 4,color = 240
oplot,vmag_star[matchSR_idl]-0.044,vmag_unscat[matchSR_idl],psym = 4,color = 240
set_plot, 'ps'
device, filename=outputbase + '/comp_mags'+outputfile+'_rep_andy.eps',/color,bits_per_pixel=8
plot,[0,-24],[0,-24],xrange = [-23,-14],yrange = [-23,-14],xtitle = 'Starburst99 V mag',ytitle = 'Unscat Sunrise V mag',linestyle = 2
oplot,vmag_star[matchSR_idl]-0.044,vmag_unscat[matchSR_idl],psym = 4,color = 240
stop

set_plot,'x'
miny = MIN([MIN(imag_leo[matchrad_idls] - imag_star[matchrad_idls]),MIN(bmag_leo[matchrad_idls] - bmag_star[matchrad_idls]),MIN(rmag_leo[matchrad_idls] - rmag_star[matchrad_idls]),-0.5])
maxy = MAX([MAX(imag_leo[matchrad_idls] - imag_star[matchrad_idls]),MAX(bmag_leo[matchrad_idls] - bmag_star[matchrad_idls]),MAX(rmag_leo[matchrad_idls] - rmag_star[matchrad_idls])])
plot,[1e6,1e14],[0,0],title = 'Magnitude Differences -- '+outputfile,xtitle = 'M200 -- Baryonic',ytitle = 'Delta Magnitude (Leo - Star99)',yrange = [miny,maxy],/xlog,thick = 2,xrange = [1e8,MAX(bm2001[matchidl_rads])]
oplot,bm2001[matchidl_rads],bmag_leo[matchrad_idls] - bmag_star[matchrad_idls],color = 12,psym = 1
oplot,bm2001[matchidl_rads],imag_leo[matchrad_idls] - imag_star[matchrad_idls],color = 190,psym = 1
oplot,bm2001[matchidl_rads],bmag_leo[matchrad_idls] - bmag_star[matchrad_idls],color = 50 ,psym = 1
oplot,bm2001[matchidl_rads],rmag_leo[matchrad_idls] - rmag_star[matchrad_idls],color = 240,psym = 1
legend,['B','V','R','I'],color=[50,120,190,240],linestyle = [0,0,0];,/right
set_plot, 'ps'
device, filename=outputbase + '/comp_mags'+outputfile+'_grids_ls.eps',/color,bits_per_pixel=8
plot,[1e6,1e14],[0,0],title = 'Magnitude Differences -- '+outputfile,xtitle = 'M200 -- Baryonic',ytitle = 'Delta Magnitude (Leo - Star99)',yrange = [miny,maxy],/xlog,thick = 2,xrange = [1e8,MAX(bm2001[matchidl_rads])]
oplot,bm2001[matchidl_rads],bmag_leo[matchrad_idls] - bmag_star[matchrad_idls],color = 50,psym = 1
oplot,bm2001[matchidl_rads],vmag_leo[matchrad_idls] - vmag_star[matchrad_idls],color = 120,psym = 1
oplot,bm2001[matchidl_rads],rmag_leo[matchrad_idls] - rmag_star[matchrad_idls],color = 190,psym = 1
oplot,bm2001[matchidl_rads],imag_leo[matchrad_idls] - imag_star[matchrad_idls],color = 240,psym = 1
legend,['B','V','R','I'],color=[50,120,190,240],linestyle = [0,0,0];,/right
device,/close
stop

set_plot,'x'
miny = MIN([MIN(imag_leo[matchall_idl] - i_unscat[matchall_SR]),MIN(bmag_leo[matchall_idl] - b_unscat[matchall_SR]),MIN(rmag_leo[matchall_idl] - r_unscat[matchall_SR]),-0.5])
maxy = MAX([MAX(imag_leo[matchall_idl] - i_unscat[matchall_SR]),MAX(bmag_leo[matchall_idl] - b_unscat[matchall_SR]),MAX(rmag_leo[matchall_idl] - r_unscat[matchall_SR])])
plot,[1e6,1e14],[0,0],title = 'Magnitude Differences -- '+outputfile,xtitle = 'M200 -- Baryonic',ytitle = 'Delta Magnitude (Leo - Unscattered Sunrise)',/xlog,thick = 2,yrange=[miny,maxy],xrange = [1e8,MAX(bm2001[matchidl_rads])]
oplot,bm2001[matchall_rad],bmag_leo[matchall_idl] - b_unscat[matchall_SR],color = 50,psym = 1
oplot,bm2001[matchall_rad],vmag_leo[matchall_idl] - v_unscat[matchall_SR],color = 120,psym = 1
oplot,bm2001[matchall_rad],rmag_leo[matchall_idl] - r_unscat[matchall_SR],color = 190,psym = 1 
oplot,bm2001[matchall_rad],imag_leo[matchall_idl] - i_unscat[matchall_SR],color = 240,psym = 1
legend,['B','V','R','I'],color=[50,150,240],linestyle = [0,0,0];,/right
set_plot, 'ps'
device, filename=outputbase + '/comp_mags'+outputfile+'_sunriseUS_leo.eps',/color,bits_per_pixel=8
plot,[1e6,1e14],[0,0],title = 'Magnitude Differences -- '+outputfile,xtitle = 'M200 -- Baryonic',ytitle = 'Delta Magnitude (Leo - Unscattered Sunrise)',/xlog,thick = 2,yrange=[-1.5,maxy],xrange = [1e8,MAX(bm2001[matchidl_rads])]
oplot,bm2001[matchall_rad],bmag_leo[matchall_idl] - b_unscat[matchall_SR],color = 50,psym = 1
oplot,bm2001[matchall_rad],vmag_leo[matchall_idl] - v_unscat[matchall_SR],color = 120,psym = 1
oplot,bm2001[matchall_rad],rmag_leo[matchall_idl] - r_unscat[matchall_SR],color = 190,psym = 1 
oplot,bm2001[matchall_rad],imag_leo[matchall_idl] - i_unscat[matchall_SR],color = 240,psym = 1
legend,['B','V','R','I'],color=[50,150,240],linestyle = [0,0,0];,/right
device,/close
stop


;set_plot,'x'
;miny = MIN([MIN(imag_leo[matchall_idl] - i_scat[matchall_SR]),MIN(bmag_leo[matchall_idl] - b_scat[matchall_SR]),MIN(rmag_leo[matchall_idl] - r_scat[matchall_SR]),-0.5])
;maxy = MAX([MAX(imag_leo[matchall_idl] - i_scat[matchall_SR]),MAX(bmag_leo[matchall_idl] - b_scat[matchall_SR]),MAX(rmag_leo[matchall_idl] - r_scat[matchall_SR])])
;plot,[1e6,1e14],[0,0],title = 'Magnitude Differences -- '+outputfile,xtitle = 'M200 -- Baryonic',ytitle = 'Delta Magnitude (Leo - Scattered Sunrise)',/xlog,thick = 2,yrange=[miny,maxy],xrange = [1e8,MAX(bm2001[matchidl_rads])]
;oplot,bm2001[matchall_rad],imag_leo[matchall_idl] - i_scat[matchall_SR],color = 150,psym = 1
;oplot,bm2001[matchall_rad],bmag_leo[matchall_idl] - b_scat[matchall_SR],color = 50,psym = 1 
;oplot,bm2001[matchall_rad],rmag_leo[matchall_idl] - r_scat[matchall_SR],color = 240,psym = 1
;legend,['B','I','R'],color=[50,150,240],linestyle = [0,0,0];,/right
;set_plot, 'ps'
;device, filename=outputbase + '/comp_mags'+outputfile+'_sunriseS_leo.eps',/color,bits_per_pixel=8
;plot,[1e6,1e14],[0,0],title = 'Magnitude Differences -- '+outputfile,xtitle = 'M200 -- Baryonic',ytitle = 'Delta Magnitude (Leo - Scattered Sunrise)',/xlog,thick = 2,yrange=[-1.5,maxy],xrange = [1e8,MAX(bm2001[matchidl_rads])]
;oplot,bm2001[matchall_rad],imag_leo[matchall_idl] - i_scat[matchall_SR],color = 150,psym = 1
;oplot,bm2001[matchall_rad],bmag_leo[matchall_idl] - b_scat[matchall_SR],color = 50,psym = 1 
;oplot,bm2001[matchall_rad],rmag_leo[matchall_idl] - r_scat[matchall_SR],color = 240,psym = 1
;legend,['B','I','R'],color=[50,150,240],linestyle = [0,0,0];,/right
;device,/close
;stop

set_plot,'x'
miny = MIN([MIN(imag_star[matchall_idl] - i_unscat[matchall_SR]),MIN(bmag_star[matchall_idl] - b_unscat[matchall_SR]),MIN(rmag_star[matchall_idl] - r_unscat[matchall_SR]),-0.5])
maxy = MAX([MAX(imag_star[matchall_idl] - i_unscat[matchall_SR]),MAX(bmag_star[matchall_idl] - b_unscat[matchall_SR]),MAX(rmag_star[matchall_idl] - r_unscat[matchall_SR])])
plot,[1e6,1e14],[0,0],title = 'Magnitude Differences -- '+outputfile,xtitle = 'M200 -- Baryonic',ytitle = 'Delta Magnitude (Star99 - Unscattered Sunrise)',/xlog,thick = 2,yrange=[miny,maxy],xrange = [1e8,MAX(bm2001[matchall_rad])]
oplot,bm2001[matchall_rad],vmag_star[matchall_idl] - v_unscat[matchall_SR],color = 120,psym = 1
oplot,bm2001[matchall_rad],imag_star[matchall_idl] - i_unscat[matchall_SR],color = 190,psym = 1
oplot,bm2001[matchall_rad],bmag_star[matchall_idl] - b_unscat[matchall_SR],color = 50,psym = 1 
oplot,bm2001[matchall_rad],rmag_star[matchall_idl] - r_unscat[matchall_SR],color = 240,psym = 1
legend,['B','I','R'],color=[50,150,240],linestyle = [0,0,0];,/right
set_plot, 'ps'
device, filename=outputbase + '/comp_mags'+outputfile+'_sunriseUS_s99.eps',/color,bits_per_pixel=8
plot,[1e6,1e14],[0,0],title = 'Magnitude Differences -- '+outputfile,xtitle = 'M200 -- Baryonic',ytitle = 'Delta Magnitude (Star99 - Unscattered Sunrise)',/xlog,thick = 2,yrange=[-1.5,maxy],xrange = [1e8,MAX(bm2001[matchall_rad])]
oplot,bm2001[matchall_rad],vmag_star[matchall_idl] - v_unscat[matchall_SR],color = 120,psym = 1
oplot,bm2001[matchall_rad],imag_star[matchall_idl] - i_unscat[matchall_SR],color = 190,psym = 1
oplot,bm2001[matchall_rad],bmag_star[matchall_idl] - b_unscat[matchall_SR],color = 50,psym = 1 
oplot,bm2001[matchall_rad],rmag_star[matchall_idl] - r_unscat[matchall_SR],color = 240,psym = 1
legend,['B','I','R'],color=[50,150,240],linestyle = [0,0,0];,/right
device,/close
stop

;set_plot,'x'
;miny = MIN([MIN(imag_star[matchall_idl] - i_scat[matchall_SR]),MIN(bmag_star[matchall_idl] - b_scat[matchall_SR]),MIN(rmag_star[matchall_idl] - r_scat[matchall_SR]),-0.5])
;maxy = MAX([MAX(imag_star[matchall_idl] - i_scat[matchall_SR]),MAX(bmag_star[matchall_idl] - b_scat[matchall_SR]),MAX(rmag_star[matchall_idl] - r_scat[matchall_SR])])
;plot,[1e6,1e14],[0,0],title = 'Magnitude Differences -- '+outputfile,xtitle = 'M200 -- Baryonic',ytitle = 'Delta Magnitude (Star99 - Scattered Sunrise)',psym = 1,/xlog,thick = 2,yrange=[miny,maxy],xrange = [1e8,MAX(bm2001[matchall_rad])]
;oplot,bm2001[matchall_rad],imag_star[matchall_idl] - i_scat[matchall_SR],color = 150,psym = 1
;oplot,bm2001[matchall_rad],bmag_star[matchall_idl] - b_scat[matchall_SR],color = 50,psym = 1 
;oplot,bm2001[matchall_rad],rmag_star[matchall_idl] - r_scat[matchall_SR],color = 240,psym = 1
;oplot,[1e6,1e14],[0,0]
;legend,['B','I','R'],color=[50,150,240],linestyle = [0,0,0];,/right
;set_plot, 'ps'
;device, filename=outputbase + '/comp_mags'+outputfile+'_sunriseS_s99.eps',/color,bits_per_pixel=8
;plot,[1e6,1e14],[0,0],title = 'Magnitude Differences -- '+outputfile,xtitle = 'M200 -- Baryonic',ytitle = 'Delta Magnitude (Star99 - Scattered Sunrise)',psym = 1,/xlog,thick = 2,yrange=[-1.5,maxy],xrange = [1e8,MAX(bm2001[matchall_rad])]
;oplot,bm2001[matchall_rad],imag_star[matchall_idl] - i_scat[matchall_SR],color = 150,psym = 1
;oplot,bm2001[matchall_rad],bmag_star[matchall_idl] - b_scat[matchall_SR],color = 50,psym = 1 
;oplot,bm2001[matchall_rad],rmag_star[matchall_idl] - r_scat[matchall_SR],color = 240,psym = 1
;oplot,[1e6,1e14],[0,0]
;legend,['B','I','R'],color=[50,150,240],linestyle = [0,0,0];,/right
;device,/close
;stop


;*****************Tully Fisher Plots***************************************************************
mags = findgen((24.0-12.0)/0.1)*0.1 -24.0
slope = -0.128412
slope = -0.13176470588
velocities = slope*(mags + 15.5) + 1.5 ;1.12 /-8.5
set_plot,'x'
plot,velocities,mags,title = 'Tully Fisher Relation -- '+outputfile+' -- Sunrise -- Unscattered',xtitle = 'Log Velocity (km/s)',ytitle = 'MI  - 5 Log(h) ',xrange=[1.6,2.4],yrange = [-14,-23],thick = 2,pos = [0.2,0.1,0.8,0.9]
oplot,Alog10(v11[matchall_rad]),i_unscat[matchall_SR] - 5*ALOG10(0.7),color = 50,psym = 4
oplot,Alog10(v21[matchall_rad]),i_unscat[matchall_SR] - 5*ALOG10(0.7),color = 240,psym = 4
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle = [0,0,0]
set_plot, 'ps'
device, filename=outputbase + '/'+outputfile+'_TF_sunrise.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation -- '+outputfile+' -- Sunrise -- Unscattered',xtitle = 'Log Velocity (km/s)',ytitle = 'MI  - 5 Log(h) ',xrange=[1.6,2.4],yrange = [-14,-23],thick = 2,pos = [0.2,0.1,0.8,0.9]
oplot,Alog10(v11[matchall_rad]),i_unscat[matchall_SR] - 5*ALOG10(0.7),color = 50,psym = 4
oplot,Alog10(v21[matchall_rad]),i_unscat[matchall_SR] - 5*ALOG10(0.7),color = 240,psym = 4
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle = [0,0,0]
stop

set_plot,'x'
plot,velocities,mags,title = 'Tully Fisher Relation -- '+outputfile+' -- Sunrise -- Scattered',xtitle = 'Log Velocity (km/s)',ytitle = 'MI  - 5 Log(h) ',xrange=[1.6,2.4],yrange = [-14,-23],thick = 2,pos = [0.2,0.1,0.8,0.9]
oplot,Alog10(v11[matchall_rad]),i_scat[matchall_SR] - 5*ALOG10(0.7),color = 50,psym = 4
oplot,Alog10(v21[matchall_rad]),i_scat[matchall_SR] - 5*ALOG10(0.7),color = 240,psym = 4
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle = [0,0,0]
set_plot, 'ps'
device, filename=outputbase + '/'+outputfile+'_TF_sunrise_scatter.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation -- '+outputfile+' -- Sunrise -- Scattered',xtitle = 'Log Velocity (km/s)',ytitle = 'MI  - 5 Log(h) ',xrange=[1.6,2.4],yrange = [-14,-23],thick = 2,pos = [0.2,0.1,0.8,0.9]
oplot,Alog10(v11[matchall_rad]),i_scat[matchall_SR] - 5*ALOG10(0.7),color = 50,psym = 4
oplot,Alog10(v21[matchall_rad]),i_scat[matchall_SR] - 5*ALOG10(0.7),color = 240,psym = 4
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle = [0,0,0]
device,/close
stop

set_plot,'x'
plot,velocities,mags,title = 'Tully Fisher Relation -- '+outputfile+' -- Leo',xtitle = 'Log Velocity (km/s)',ytitle = 'MI  - 5 Log(h) ',xrange=[1.4,2.6],yrange = [-12,-26],thick = 2,pos = [0.2,0.1,0.8,0.9]
oplot,Alog10(v11[matchall_rad]),imag_leo[matchall_idl] - 5*ALOG10(0.7),color = 50,psym = 4
oplot,Alog10(v21[matchall_rad]),imag_leo[matchall_idl] - 5*ALOG10(0.7),color = 240,psym = 4
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle = [0,0,0]
set_plot, 'ps'
device, filename=outputbase + '/'+outputfile+'_TF_leo.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation -- '+outputfile+' -- Leo',xtitle = 'Log Velocity (km/s)',ytitle = 'MI  - 5 Log(h) ',xrange=[1.4,2.6],yrange = [-12,-26],thick = 2,pos = [0.2,0.1,0.8,0.9]
oplot,Alog10(v11[matchall_rad]),imag_leo[matchall_idl] - 5*ALOG10(0.7),color = 50,psym = 4
oplot,Alog10(v21[matchall_rad]),imag_leo[matchall_idl] - 5*ALOG10(0.7),color = 240,psym = 4
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle = [0,0,0]
stop

set_plot,'x'
plot,velocities,mags,title = 'Tully Fisher Relation -- '+outputfile+' -- Starburst99',xtitle = 'Log Velocity (km/s)',ytitle = 'MI  - 5 Log(h) ',xrange=[1.4,2.6],yrange = [-12,-26],thick = 2,pos = [0.2,0.1,0.8,0.9]
oplot,Alog10(v11[matchall_rad]),imag_leo[matchall_idl] - 5*ALOG10(0.7),color = 50,psym = 4
oplot,Alog10(v21[matchall_rad]),imag_leo[matchall_idl] - 5*ALOG10(0.7),color = 240,psym = 4
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle = [0,0,0]
set_plot, 'ps'
device, filename=outputbase + '/'+outputfile+'_TF_star99.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation -- '+outputfile+' -- Starburst99',xtitle = 'Log Velocity (km/s)',ytitle = 'MI  - 5 Log(h) ',xrange=[1.4,2.6],yrange = [-12,-26],thick = 2,pos = [0.2,0.1,0.8,0.9]
oplot,Alog10(v11[matchall_rad]),imag_leo[matchall_idl] - 5*ALOG10(0.7),color = 50,psym = 4
oplot,Alog10(v21[matchall_rad]),imag_leo[matchall_idl] - 5*ALOG10(0.7),color = 240,psym = 4
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle = [0,0,0]

stop
END
