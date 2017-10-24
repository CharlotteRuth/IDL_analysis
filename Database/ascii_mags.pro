pro ascii_mags, outputfile
loadct,39
set_plot,'x'
mag_file_star = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/h1201-h1755-'+outputfile+'g1bwK.00512.amiga_vir.halos.star99_ab.Mv.fits'
;mag_file_star = '../Fabio_Runs/h1201-h1755-'+outputfile+'bwK/h1201-h1755-'+outputfile+'bwK.00512.amiga_vir.halos.star99_ab.Mv.fits'
;mag_file_star = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/h1201-h1755-'+outputfile+'g1bwK.00084.amiga_vir.halos.star99_ab.Mv.fits'
halo_mags_star = MRDFITS(mag_file_star,1)
slice=where(halo_mags_star.id NE -1)
id_star = halo_mags_star[slice].id
ids = findgen(max(id_star)+1)
data = fltarr(N_ELEMENTS(ids),6)
vmag_star = halo_mags_star[slice].v
data[*,0] = ids                     ;ID
data[id_star,1] = vmag_star         ;Star99 V magnitude

mag_file_leo = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/h1201-h1755-'+outputfile+'g1bwK.00512.amiga_vir.halos.leo_ab.Mv.fits'
;mag_file_leo ='../Fabio_Runs/h1201-h1755-'+outputfile+'bwK/h1201-h1755-'+outputfile+'bwK.00512.amiga_vir.halos.leo_ab.Mv.fits'
;;mag_file_leo = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/h1201-h1755-'+outputfile+'g1bwK.00084.amiga_vir.halos.leo_ab.Mv.fits'
halo_mags_leo = MRDFITS(mag_file_leo,1)
vmag_leo = halo_mags_leo[slice].v
data[id_star,2] = vmag_leo          ;Leo V magnitude

mag_file_unscat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/512/X2.unscat.mag'
;mag_file_unscat = '../Fabio_Runs/h1201-h1755-'+outputfile+'bwK/Sunrise/512/X2g3.unscat.mag'
;mag_file_unscat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/512/X3.unscat.mag'
;mag_file_unscat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/512/X5.unscat.mag'
;mag_file_unscat = '../Fabio_Runs/h1201-h1755-'+outputfile+'g1bwK/Sunrise/084/X5.unscat.mag'
readcol,mag_file_unscat,id_unscat,u_unscat,b_unscat,v_unscat,r_unscat,i_unscat,us_unscat,gs_unscat,rs_unscat,is_unscat,zs_unscat,i1_unscat,i2_unscat,m2_unscat,m7_unscat,m1_unscat,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
data[id_unscat,3] = v_unscat
data[*,4] = data[*,1] - data[*,3]
data[*,5] = data[*,2] - data[*,3]
intersect = WHERE(data[*,1] ne 0 AND data[*,3] ne 0)
data_merge = Transpose(data[intersect,*])
plot,data_merge[0,*],data_merge[4,*];,yrange=[-0.75,1.5]
oplot,data_merge[0,*],data_merge[5,*],linestyle = 2
oplot,data_merge[0,*],data_merge[1,*] - data_merge[2,*],color=240,linestyle=1
print,data_merge 
stop
filename = 'Data_files/'+outputfile+'randv.dat'
;filename = 'Data_files/'+outputfile+'_z3_randv.dat'
readcol,filename,id_rad,r2001, v2001, sm2001, bm2001, r11, v11, sm11, bm11, r21, v21, sm21, bm21,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
FOR ct = 0, n_elements(id_rad)-1 DO BEGIN
    IF (strpos(id_rad[ct],'0') EQ 0) THEN id_rad[ct] = STRMID(id_rad[ct],1,STRLEN(id_rad[ct])-1)
    IF (strpos(id_rad[ct],'0') EQ 0) THEN id_rad[ct] = STRMID(id_rad[ct],1,STRLEN(id_rad[ct])-1)
    IF (strpos(id_rad[ct],'0') EQ 0) THEN id_rad[ct] = STRMID(id_rad[ct],1,STRLEN(id_rad[ct])-1)
ENDFOR
halo_data = replicate({id:0, b:0., v:0., i:0., r:0., k:0., u:0., r200:0., v200:0., r200_11:0, v200_11:0},N_ELEMENTS(ids))
halo_data.id = ids
halo_data[id_star].b = halo_mags_star[slice].b
halo_data[id_star].v = halo_mags_star[slice].v
halo_data[id_star].i = halo_mags_star[slice].i
halo_data[id_star].r = halo_mags_star[slice].r
halo_data[id_star].k = halo_mags_star[slice].k
halo_data[id_star].u = halo_mags_star[slice].u

halo_data[id_rad].r200 = r2001
halo_data[id_rad].v200 = v2001
halo_data[id_rad].r200_11 = r11
halo_data[id_rad].v200_11 = v11
intersect = WHERE(halo_data.r200 ne 0 OR halo_data.v ne 0)
data_merge = halo_data[intersect]
outfile=outputfile+'_Star99_ABmags'
;outfile=outputfile+'z3__Star99_ABmags'
writecol,outfile+'.dat',data_merge.id,data_merge.b,data_merge.v,data_merge.i,data_merge.r,data_merge.k,data_merge.u,data_merge.r200,data_merge.v200,data_merge.r200_11,data_merge.v200_11,format='(I,F,F,F,F,F,F,F,F,F,F)'
MWRFITS, data_merge, outfile+'.fits', /create
end
