PRO compShear
;This program compares the shears for two sepporate simulations
loadct,39
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23

msol_per_sysmass = 1.36e17
kpc_per_syslength = 1e5
msol_per_sysmass       = 2.310e15
kpc_per_syslength        = 25000.
dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength^3/cm_per_kpc^3

;oldfile = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e4_zsol/H2_split_oldshear/MW_disk.00010'
;newfile = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e4_zsol/H2_split_newshear/MW_disk.00010'

;dir = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_repl/z1_lw3e7_cp10/'
;oldfile = dir + 'MW_disk.00010.curl'
;newfile = dir + 'MW_disk.00010.turb'

;dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/'
;oldfile = dir + 'h516.cosmo25cmb.1536g3HBWK.00512.halo.1.curl'
;newfile = dir + 'h516.cosmo25cmb.1536g3HBWK.00512.halo.1.turb'

;dir = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e4_zsol/'
;oldfile = dir + 'H2_split_oldshear/MW_disk.00100.turb'
;newfile = dir + 'H2_split_newshear/MW_disk.00100.turb'

dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g7HBWK_00480/'
oldfile = 'steps.JeansSF/h516.cosmo25cmb.1536g8HBWK.JeansSF.00408.00011.dir/h516.cosmo25cmb.1536g8HBWK.JeansSF.00408.00011'
newfile = 'steps.JeansSF.newShear/h516.cosmo25cmb.1536g8HBWK.JeansSF.newShear.00408.00011.dir/h516.cosmo25cmb.1536g8HBWK.JeansSF.newShear.00408.00011'
file3 = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00420.dir/h516.cosmo25cmb.1536g3HBWK.00420'

rtipsy,oldfile,h1,g1,d1,s1
g1.dens = g1.dens*dens_convert
rtipsy,newfile,h2,g2,d2,s2
g2.dens = g2.dens*dens_convert
rtipsy,file3,h3,g3,d3,s3
g3.dens = g3.dens*dens_convert
;rtipsy,dir + 'H2_split_oldshear/MW_disk.00010',h1,g1,d1,s1
;rtipsy,dir + 'H2_split_newshear/MW_disk.00010',h2,g2,d2,s2
;rtipsy,dir + 'h516.cosmo25cmb.1536g3HBWK.00512.halo.1',h1,g1,d1,s1
;rtipsy,dir + 'h516.cosmo25cmb.1536g3HBWK.00512.halo.1',h2,g2,d2,s2
;readarr,oldfile+'.shear',h1,shear1,part = 'gas',/ascii
;readarr,newfile+'.shear',h2,shear2,part = 'gas',/ascii
readarr,oldfile+'.correL',h1,correL1,part = 'gas',/ascii
correL1 = correL1 * kpc_per_syslength
readarr,newfile+'.correL',h2,correL2,part = 'gas',/ascii
correL2 = correL2 * kpc_per_syslength
readarr,file3+'.shear',h3,correL3,part = 'gas',/ascii
correL3 = (1.0/correL3*kpc_per_syslength)
readarr,oldfile+'.smoothlength',h1,smooth1,part = 'gas',/ascii
smooth1 = smooth1 * kpc_per_syslength
readarr,newfile+'.smoothlength',h2,smooth2,part = 'gas',/ascii
smooth2 = smooth2 * kpc_per_syslength
readarr,file3+'.smoothlength',h3,smooth3,part = 'gas',/ascii
smooth3 = smooth3 * kpc_per_syslength
;readarr,dir +'h516.cosmo25cmb.1536g3HBWK.00512.halo.1.smoothlength',h1,smooth1,part = 'gas',/ascii
;shear2_base = shear2
;shear2 = shear2/smooth2/smooth2
readarr,newfile+'.c',h2,sound2,part = 'gas',/ascii
readarr,oldfile+'.c',h1,sound1,part = 'gas',/ascii
;shear1a = shear1*sound1
;shear2a = shear2*sound2


plot,g1.dens,correL1*cm_per_kpc*h1.time,psym = 3,/xlog,/ylog,xrange = [1e-6,1e2]
oplot,g2.dens,correL2*cm_per_kpc*h2.time,psym = 3,color = 240
oplot,g1.dens,smooth1*cm_per_kpc*h1.time,psym = 3,color = 50
oplot,g2.dens,smooth2*cm_per_kpc*h2.time,psym = 3,color = 100
oplot,g3.dens,correL3*cm_per_kpc*h3.time,psym = 3,color = 140
oplot,g3.dens,smooth3*cm_per_kpc*h3.time,psym = 3,color = 200
stop


plot,sqrt(g1.x*g1.x + g1.y*g1.y)*kpc_per_syslength,shear1,psym = 3,/ylog,xtitle = 'Radius [kpc]',ytitle = 'Shear',yrange = [1e3,1e9],xrange = [0,20]
oplot,sqrt(g2.x*g2.x + g2.y*g2.y)*kpc_per_syslength,shear2,psym = 3,color = 240
stop

plot,sqrt(g1.x*g1.x + g1.y*g1.y)*kpc_per_syslength,shear2/shear1,psym = 3,/ylog,xtitle = 'Radius [kpc]',ytitle = 'Ratio of Shear',xrange = [0,20]
ind = where(g1.dens*dens_convert gt 1)
oplot,sqrt(g2[ind].x*g2[ind].x + g2[ind].y*g2[ind].y)*kpc_per_syslength,shear2[ind]/shear1[ind],psym = 3,color = 190
ind = where(g1.dens*dens_convert gt 5)
oplot,sqrt(g2[ind].x*g2[ind].x + g2[ind].y*g2[ind].y)*kpc_per_syslength,shear2[ind]/shear1[ind],psym = 3,color = 240
stop
window,0
plot,g1.dens*dens_convert,kpc_per_syslength*cm_per_kpc/shear1,psym = 3,/xlog,/ylog,yrange = [1e18,1e23],xtitle = 'Density',ytitle = 'Length',xrange = [1e-6,1e2]
;oplot,g1.dens*dens_convert,kpc_per_syslength*cm_per_kpc*smooth1,psym = 3,color = 100
oplot,g2.dens*dens_convert,kpc_per_syslength*cm_per_kpc/shear2,psym = 3,color = 240
;oplot,g2.dens*dens_convert,kpc_per_syslength*cm_per_kpc*smooth2,psym = 3,color =190
stop
window,1
plot,g1.dens*dens_convert,g1.dens*dens_convert*kpc_per_syslength*cm_per_kpc/shear1,psym = 3,/xlog,/ylog,xtitle = 'Density',ytitle = 'Column Density',xrange = [1e-6,1e2]
;oplot,g1.dens*dens_convert,g1.dens*dens_convert*kpc_per_syslength*cm_per_kpc*smooth1,psym = 3,color = 100
oplot,g2.dens*dens_convert,g2.dens*dens_convert*kpc_per_syslength*cm_per_kpc/shear2,psym = 3,color = 240
;oplot,g2.dens*dens_convert,g2.dens*dens_convert*kpc_per_syslength*cm_per_kpc*smooth2,psym=
;3,color = 190

readarr,dir + 'H2_split_newshear/MW_disk.00010.H2',h1,H2,part = 'gas',/ascii
plot,g1.dens*dens_convert*kpc_per_syslength*cm_per_kpc/shear1,H2,psym = 3,/ylog,/xlog
stop
END
