;Like the first part of allgas_z0.pro, this program identifies gas
;that was once in two halos. Unlike it, this program does it more efficiently.

PRO allgas_z0_loop_master
dir = '/nobackupp8/crchrist/MolecH/storm.cosmo25cmb.4096g/storm.cosmo25cmb.4096g5HbwK1BH/' 
filebase='storm.cosmo25cmb.4096g5HbwK1BH' 
haloids= ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60']

;dir = '/nobackupp8/crchrist/MolecH/elektra.cosmo25cmb.4096g/elektra.cosmo25cmb.4096g5HbwK1BH/'
;filebase='elektra.cosmo25cmb.4096g5HbwK1BH'
;haloids=['1','2','3','4','5','8','9','10','11','12','17','18']

dir='/nobackupp8/crchrist/MolecH/rogue.cosmo25cmb.4096g/rogue.cosmo25cmb.4096g5HbwK1BH/'
filebase='rogue.cosmo25cmb.4096g5HbwK1BH'
haloids = ['1','3','7','8','10','11','12','16','17','18','30','31','32','34']

;dir='/nobackupp8/crchrist/MolecH/cptmarvel.cosmo25.4096g/cptmarvel.cosmo25cmb.4096g5HbwK1BH/'
;filebase='cptmarvel.cosmo25cmb.4096g5HbwK1BH'
;haloids = ['1','2','4','5','6','7','10','11','13','14']

dir = '/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK'
filebase = 'h986.cosmo50cmb.3072g14HBWK'
step = '00512'
haloids = ['1','2','3','8', '15','16'] ;11, 14

dir = '/nobackupp8/crchrist/MolecH/h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/'
filebase = 'h242.cosmo50PLK.3072gst5HbwK1BH'
step = '004096'
haloids = ['1','9','11','24','29','33','39','40','45','72','75','76','425','457']
haloids = ['1','33','45','72'] ;Not satellites

dir = '/nobackupp8/crchrist/MolecH/h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/' 
filebase='h229.cosmo50PLK.3072gst5HbwK1BH' 
haloids = ['1','2','4','7','17','21','22','27','51','52','70','104','203']
haloids = ['1','2','4','7','17'] ;Not satellites
step = '004096'


allgas_z0_loop,dir,filebase,step,haloids
END

PRO allgas_z0_loop,dir,filebase,step,haloids
cd,dir
;rtipsy,filebase+'.'+step+'/'+filebase+'.'+step,h,g,d,s,/justhead
;iord = read_lon_array(filebase+'.'+step+'/'+filebase+'.'+step+'.iord')
;amigastat = read_stat_struc_amiga(filebase+'.'+step+'/'+filebase+'.'+step+'.amiga.stat')
units = tipsyunits(filebase+'.param')
print,'1, grp'+haloids[0]
allgas1 = mrdfits('grp'+haloids[0]+'.allgas.iord.fits',/silent)
allgas = allgas1
halos = intarr(n_elements(allgas1)) + fix(haloids[0])
;match,iord,allgas1,ind1_tip,ind1
;readcol,filebase + '.grp'+haloids[0]+'.haloid.dat',files1,pasthalos1,format='(A,A)'

FOR i = 1, n_elements(haloids) - 1 DO BEGIN
    print,strtrim(i+1,2),', grp'+haloids[i]
    allgasnew = mrdfits('grp'+haloids[i]+'.allgas.iord.fits',/silent) 
    match,allgas,allgasnew,temp1,temp2
    IF temp1[0] NE -1 THEN BEGIN
;        print,strtrim((halos[temp1])[uniq(halos[temp1],sort(halos[temp1]))],2)+' - '+strtrim(haloids[i],2)
        primehalos = (halos[temp1])[uniq(halos[temp1],sort(halos[temp1]))]
        FOR iprimehalo = 0, n_elements(primehalos) - 1 DO BEGIN
            print,strtrim(primehalos[iprimehalo],2)+' - '+strtrim(haloids[i],2)
            print,'Intersection: ',strtrim(n_elements(where(halos[temp1] EQ primehalos[iprimehalo])),2),'; in halo 1: ',strtrim(n_elements(where(halos EQ primehalos[iprimehalo])),2),', in halo 2: ',strtrim(n_elements(allgasnew),2)
        ENDFOR
    ENDIF
    allgas = [allgas,allgasnew]
    halos = [halos,intarr(n_elements(allgasnew)) + fix(haloids[i])]
ENDFOR


END
