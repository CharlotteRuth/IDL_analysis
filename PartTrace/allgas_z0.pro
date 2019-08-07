;

PRO MASTER

dir = '/nobackupp8/crchrist/MolecH/cptmarvel.cosmo25.4096g/cptmarvel.cosmo25cmb.4096g5HbwK1BH'
filebase = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
step = '004096'
haloids = ['1','2','4','5','6','7','10','11','13','14']

dir = '/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK'
filebase = 'h986.cosmo50cmb.3072g14HBWK'
step = '00512'
haloids = ['1','2','3','8', '15','16'] ;11, 14

;dir = '/nobackupp8/crchrist/MolecH/elektra.cosmo25cmb.4096g/elektra.cosmo25cmb.4096g5HbwK1BH/'
;filebase='elektra.cosmo25cmb.4096g5HbwK1BH'
;haloids=['1','2','3','4','5','8','9','10','11','12','17','18']
;haloids_int = ['1','5','11']
;step = '004096'

;dir='/nobackupp8/crchrist/MolecH/rogue.cosmo25cmb.4096g/rogue.cosmo25cmb.4096g5HbwK1BH/'
;filebase='rogue.cosmo25cmb.4096g5HbwK1BH'
;haloids = ['1','3','7','8','10','11','12','16','17','18','30','31','32','34','36']
;haloids_int = ['1','3','7','8','10','32','34']
;step = '004096'

;dir = '/nobackupp8/crchrist/MolecH/storm.cosmo25cmb.4096g/storm.cosmo25cmb.4096g5HbwK1BH/' 
;filebase='storm.cosmo25cmb.4096g5HbwK1BH' 
;haloids= ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60']
;haloids_int = ['1','4','8','35','43','49']
;step = '004096'

dir = '/nobackupp8/crchrist/MolecH/h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/'
filebase='h242.cosmo50PLK.3072gst5HbwK1BH' 
haloids= ['1','9','11','24','29','33','39','40','45','72','75','76','425','457'] 
haloids = ['1','33','45','72'] ;Limited to non-satellites
step = '004096'

dir = '/nobackupp8/crchrist/MolecH/h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/' 
filebase='h229.cosmo50PLK.3072gst5HbwK1BH' 
haloids = ['1','2','4','7','17','21','22','27','51','52','70','104','203']
haloids = ['1','2','4','7','17'] ;Not satellites
step = '004096'

allgas_z0,dir,filebase,step,haloids
END

PRO allgas_z0,dir,filebase,step,haloids
formatplot
loadct,39

cd,dir
rtipsy,filebase+'.'+step+'/'+filebase+'.'+step,h,g,d,s,/justhead
iord = read_lon_array(filebase+'.'+step+'/'+filebase+'.'+step+'.iord')
grp = read_lon_array(filebase+'.'+step+'/'+filebase+'.'+step+'.amiga.grp')
amigastat = read_stat_struc_amiga(filebase+'.'+step+'/'+filebase+'.'+step+'.amiga.stat')
units = tipsyunits(filebase+'.param')
print,'1, grp'+haloids[0]
allgas1 = mrdfits('grp'+haloids[0]+'.allgas.iord.fits')
match,iord,allgas1,ind1_tip,ind1
readcol,filebase + '.grp'+haloids[0]+'.haloid.dat',files1,pasthalos1,format='(A,A)'
IF n_elements(haloids) GE 2 THEN BEGIN & $
  allgas2 = mrdfits('grp'+haloids[1]+'.allgas.iord.fits') &  $
  print,'2, grp'+haloids[1] & $
  match,iord,allgas2,ind2_tip,ind2 & $
  match,allgas2,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'2-1',n_elements(temp1),n_elements(allgas2),n_elements(allgas1) & $
  readcol,filebase + '.grp'+haloids[1]+'.haloid.dat',files2,pasthalos2,format='(A,A)' & $
  print,'Unique halos at final step: ',(grp[ind2_tip])[uniq(grp[ind2_tip],sort(grp[ind2_tip]))] & $
ENDIF
IF n_elements(haloids) GE 3 THEN BEGIN & $
  allgas3 = mrdfits('grp'+haloids[2]+'.allgas.iord.fits')& $ 
  print,'3, grp'+haloids[2] & $
  match,iord,allgas3,ind3_tip,ind3 & $
  match,allgas3,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'3-1',n_elements(temp1),n_elements(allgas3),n_elements(allgas1) & $
  match,allgas3,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'3-2',n_elements(temp1),n_elements(allgas3),n_elements(allgas2) & $
  readcol,filebase + '.grp'+haloids[2]+'.haloid.dat',files3,pasthalos3,format='(A,A)' & $
ENDIF
IF n_elements(haloids) GE 4 THEN BEGIN & $
  allgas4 = mrdfits('grp'+haloids[3]+'.allgas.iord.fits') & $
  print,'4, grp'+haloids[3] & $
  match,iord,allgas4,ind4_tip,ind4 & $
  match,allgas4,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'4-1',n_elements(temp1),n_elements(allgas4),n_elements(allgas1) & $
  match,allgas4,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'4-2',n_elements(temp1),n_elements(allgas4),n_elements(allgas2) & $
  match,allgas4,allgas3,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'4-3',n_elements(temp1),n_elements(allgas4),n_elements(allgas3) & $
  readcol,filebase + '.grp'+haloids[3]+'.haloid.dat',files4,pasthalos4,format='(A,A)' & $
ENDIF
  IF n_elements(haloids) GE 5 THEN BEGIN & $
  allgas5  = mrdfits('grp'+haloids[4]+'.allgas.iord.fits') & $
  print,'5, grp'+haloids[4] & $ 
  match,iord,allgas5 ,ind5_tip,ind5 & $
  match,allgas5,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'5-1',n_elements(temp1),n_elements(allgas5),n_elements(allgas1) & $
  match,allgas5,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'5-2',n_elements(temp1),n_elements(allgas5),n_elements(allgas2) & $
  match,allgas5,allgas3,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'5-3',n_elements(temp1),n_elements(allgas5),n_elements(allgas3) & $
  match,allgas5,allgas4,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'5-4',n_elements(temp1),n_elements(allgas5),n_elements(allgas4) & $
  readcol,filebase + '.grp'+haloids[4]+'.haloid.dat',files5,pasthalos5,format='(A,A)' & $
ENDIF
IF n_elements(haloids) GE 6 THEN BEGIN & $
  allgas6 = mrdfits('grp'+haloids[5]+'.allgas.iord.fits') & $
  print,'6, grp'+haloids[5] & $ 
  match,iord,allgas6,ind6_tip,ind6 & $
  match,allgas6,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'6-1',n_elements(temp1),n_elements(allgas6),n_elements(allgas1) & $
  match,allgas6,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'6-2',n_elements(temp1),n_elements(allgas6),n_elements(allgas2) & $
  match,allgas6,allgas3,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'6-3',n_elements(temp1),n_elements(allgas6),n_elements(allgas3) & $
  match,allgas6,allgas4,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'6-4',n_elements(temp1),n_elements(allgas6),n_elements(allgas4) & $
  match,allgas6,allgas5,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'6-5',n_elements(temp1),n_elements(allgas6),n_elements(allgas5) & $
  readcol,filebase + '.grp'+haloids[5]+'.haloid.dat',files6,pasthalos6,format='(A,A)' & $
ENDIF
IF n_elements(haloids) GE 7 THEN BEGIN & $
  allgas7 = mrdfits('grp'+haloids[6]+'.allgas.iord.fits') & $
  print,'7, grp'+haloids[6] & $ 
  match,iord,allgas7,ind7_tip,ind7 & $
  match,allgas7,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'7-1',n_elements(temp1),n_elements(allgas7),n_elements(allgas1) & $
  match,allgas7,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'7-2',n_elements(temp1),n_elements(allgas7),n_elements(allgas2) & $
  match,allgas7,allgas3,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'7-3',n_elements(temp1),n_elements(allgas7),n_elements(allgas3) & $
  match,allgas7,allgas4,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'7-4',n_elements(temp1),n_elements(allgas7),n_elements(allgas4) & $
  match,allgas7,allgas5,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'7-5',n_elements(temp1),n_elements(allgas7),n_elements(allgas5) & $
  match,allgas7,allgas6,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'7-6',n_elements(temp1),n_elements(allgas7),n_elements(allgas6) & $
  readcol,filebase + '.grp'+haloids[6]+'.haloid.dat',files7,pasthalos7,format='(A,A)' & $
ENDIF
IF n_elements(haloids) GE 8 THEN BEGIN & $
  allgas8 = mrdfits('grp'+haloids[7]+'.allgas.iord.fits') & $
  print,'8, grp'+haloids[7] & $ 
  match,iord,allgas8,ind8_tip,ind8 & $
  match,allgas8,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'8-1',n_elements(temp1),n_elements(allgas8),n_elements(allgas1) & $
  match,allgas8,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'8-2',n_elements(temp1),n_elements(allgas8),n_elements(allgas2) & $
  match,allgas8,allgas3,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'8-3',n_elements(temp1),n_elements(allgas8),n_elements(allgas3) & $
  match,allgas8,allgas4,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'8-4',n_elements(temp1),n_elements(allgas8),n_elements(allgas4) & $
  match,allgas8,allgas5,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'8-5',n_elements(temp1),n_elements(allgas8),n_elements(allgas5) & $
  match,allgas8,allgas6,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'8-6',n_elements(temp1),n_elements(allgas8),n_elements(allgas6) & $
  match,allgas8,allgas7,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'8-7',n_elements(temp1),n_elements(allgas8),n_elements(allgas7) & $
  readcol,filebase + '.grp'+haloids[7]+'.haloid.dat',files8,pasthalos8,format='(A,A)' & $
  ENDIF
IF n_elements(haloids) GE 9 THEN BEGIN & $
  allgas9 = mrdfits('grp'+haloids[8]+'.allgas.iord.fits') & $
  print,'9, grp'+haloids[8] & $ 
  match,iord,allgas9,ind9_tip,ind9 & $
  match,allgas9,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'9-1',n_elements(temp1),n_elements(allgas9),n_elements(allgas1) & $
  match,allgas9,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'9-2',n_elements(temp1),n_elements(allgas9),n_elements(allgas2) & $
  match,allgas9,allgas3,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'9-3',n_elements(temp1),n_elements(allgas9),n_elements(allgas3) & $
  match,allgas9,allgas4,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'9-4',n_elements(temp1),n_elements(allgas9),n_elements(allgas4) & $
  match,allgas9,allgas5,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'9-5',n_elements(temp1),n_elements(allgas9),n_elements(allgas5) & $
  match,allgas9,allgas6,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'9-6',n_elements(temp1),n_elements(allgas9),n_elements(allgas6) & $
  match,allgas9,allgas7,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'9-7',n_elements(temp1),n_elements(allgas9),n_elements(allgas7) & $
  match,allgas9,allgas8,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'9-8',n_elements(temp1),n_elements(allgas9),n_elements(allgas8) & $
  readcol,filebase + '.grp'+haloids[8]+'.haloid.dat',files9,pasthalos9,format='(A,A)' & $
  ENDIF
IF n_elements(haloids) GE 10 THEN BEGIN & $
  allgas10 = mrdfits('grp'+haloids[9]+'.allgas.iord.fits') & $
  print,'10, grp'+haloids[9] & $ 
  match,iord,allgas10,ind10_tip,ind10 & $
  match,allgas10,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'10-1',n_elements(temp1),n_elements(allgas10),n_elements(allgas1) & $
  match,allgas10,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'10-2',n_elements(temp1),n_elements(allgas10),n_elements(allgas2) & $
  match,allgas10,allgas3,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'10-3',n_elements(temp1),n_elements(allgas10),n_elements(allgas3) & $
  match,allgas10,allgas4,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'10-4',n_elements(temp1),n_elements(allgas10),n_elements(allgas4) & $
  match,allgas10,allgas5,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'10-5',n_elements(temp1),n_elements(allgas10),n_elements(allgas5) & $
  match,allgas10,allgas6,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'10-6',n_elements(temp1),n_elements(allgas10),n_elements(allgas6) & $
  match,allgas10,allgas7,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'10-7',n_elements(temp1),n_elements(allgas10),n_elements(allgas7) & $
  match,allgas10,allgas8,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'10-8',n_elements(temp1),n_elements(allgas10),n_elements(allgas8) & $
  match,allgas10,allgas9,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'10-9',n_elements(temp1),n_elements(allgas10),n_elements(allgas9) & $
  readcol,filebase + '.grp'+haloids[9]+'.haloid.dat',files10,pasthalos10,format='(A,A)' & $
  ENDIF
IF n_elements(haloids) GE 11 THEN BEGIN & $
  allgas11 = mrdfits('grp'+haloids[10]+'.allgas.iord.fits') & $
  print,'11, grp'+haloids[10] & $ 
  match,iord,allgas11,ind11_tip,ind11 & $
  match,allgas11,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'11-1',n_elements(temp1),n_elements(allgas11),n_elements(allgas1) & $
  match,allgas11,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'11-2',n_elements(temp1),n_elements(allgas11),n_elements(allgas2) & $
  match,allgas11,allgas3,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'11-3',n_elements(temp1),n_elements(allgas11),n_elements(allgas3) & $
  match,allgas11,allgas4,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'11-4',n_elements(temp1),n_elements(allgas11),n_elements(allgas4) & $
  match,allgas11,allgas5,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'11-5',n_elements(temp1),n_elements(allgas11),n_elements(allgas5) & $
  match,allgas11,allgas6,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'11-6',n_elements(temp1),n_elements(allgas11),n_elements(allgas6) & $
  match,allgas11,allgas7,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'11-7',n_elements(temp1),n_elements(allgas11),n_elements(allgas7) & $
  match,allgas11,allgas8,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'11-8',n_elements(temp1),n_elements(allgas11),n_elements(allgas8) & $
  match,allgas11,allgas9,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'11-9',n_elements(temp1),n_elements(allgas11),n_elements(allgas9) & $
  match,allgas11,allgas10,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'11-10',n_elements(temp1),n_elements(allgas11),n_elements(allgas10) & $
  readcol,filebase + '.grp'+haloids[10]+'.haloid.dat',files11,pasthalos11,format='(A,A)' & $
  ENDIF
IF n_elements(haloids) GE 12 THEN BEGIN & $
  allgas12 = mrdfits('grp'+haloids[11]+'.allgas.iord.fits') & $
  print,'12, grp'+haloids[11] & $ 
  match,iord,allgas12,ind12_tip,ind12 & $
  match,allgas12,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-1',n_elements(temp1),n_elements(allgas12),n_elements(allgas1) & $
  match,allgas12,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-2',n_elements(temp1),n_elements(allgas12),n_elements(allgas2) & $
  match,allgas12,allgas3,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-3',n_elements(temp1),n_elements(allgas12),n_elements(allgas3) & $
  match,allgas12,allgas4,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-4',n_elements(temp1),n_elements(allgas12),n_elements(allgas4) & $
  match,allgas12,allgas5,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-5',n_elements(temp1),n_elements(allgas12),n_elements(allgas5) & $
  match,allgas12,allgas6,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-6',n_elements(temp1),n_elements(allgas12),n_elements(allgas6) & $
  match,allgas12,allgas7,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-7',n_elements(temp1),n_elements(allgas12),n_elements(allgas7) & $
  match,allgas12,allgas8,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-8',n_elements(temp1),n_elements(allgas12),n_elements(allgas8) & $
  match,allgas12,allgas9,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-9',n_elements(temp1),n_elements(allgas12),n_elements(allgas9) & $
  match,allgas12,allgas10,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-10',n_elements(temp1),n_elements(allgas12),n_elements(allgas10) & $
  match,allgas12,allgas11,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'12-11',n_elements(temp1),n_elements(allgas12),n_elements(allgas11) & $
  readcol,filebase + '.grp'+haloids[11]+'.haloid.dat',files12,pasthalos12,format='(A,A)' & $
  ENDIF
IF n_elements(haloids) GE 13 THEN BEGIN & $
  allgas13 = mrdfits('grp'+haloids[12]+'.allgas.iord.fits') & $
  print,'13, grp'+haloids[12] & $ 
  match,iord,allgas13,ind13_tip,ind13 & $
  match,allgas13,allgas1,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-1',n_elements(temp1),n_elements(allgas13),n_elements(allgas1) & $
  match,allgas13,allgas2,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-2',n_elements(temp1),n_elements(allgas13),n_elements(allgas2) & $
  match,allgas13,allgas3,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-3',n_elements(temp1),n_elements(allgas13),n_elements(allgas3) & $
  match,allgas13,allgas4,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-4',n_elements(temp1),n_elements(allgas13),n_elements(allgas4) & $
  match,allgas13,allgas5,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-5',n_elements(temp1),n_elements(allgas13),n_elements(allgas5) & $
  match,allgas13,allgas6,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-6',n_elements(temp1),n_elements(allgas13),n_elements(allgas6) & $
  match,allgas13,allgas7,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-7',n_elements(temp1),n_elements(allgas13),n_elements(allgas7) & $
  match,allgas13,allgas8,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-8',n_elements(temp1),n_elements(allgas13),n_elements(allgas8) & $
  match,allgas13,allgas9,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-9',n_elements(temp1),n_elements(allgas13),n_elements(allgas9) & $
  match,allgas13,allgas10,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-10',n_elements(temp1),n_elements(allgas13),n_elements(allgas10) & $
  match,allgas13,allgas11,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-11',n_elements(temp1),n_elements(allgas13),n_elements(allgas11) & $
  match,allgas13,allgas12,temp1,temp2 & $
  IF temp1[0] NE -1 THEN print,'13-12',n_elements(temp1),n_elements(allgas13),n_elements(allgas12) & $
  readcol,filebase + '.grp'+haloids[12]+'.haloid.dat',files13,pasthalos13,format='(A,A)' & $
  ENDIF

stop
colors = (indgen(n_elements(haloids)) + 1)*(254/(n_elements(haloids)))
n_steps = n_elements(files1)
time = (findgen(n_steps)+1)/n_steps*13.7
rtipsy,filebase+'.'+step+'/'+filebase+'.'+step,h,g,d,s
g.x = g.x*h.time*units.lengthunit
g.y = g.y*h.time*units.lengthunit
g.z = g.z*h.time*units.lengthunit
s.x = s.x*h.time*units.lengthunit
s.y = s.y*h.time*units.lengthunit
s.z = s.z*h.time*units.lengthunit
loadct,39
xminmax = [-0.1,0.08]*h.time*units.lengthunit
yminmax = [-0.12,0.02]*h.time*units.lengthunit
zminmax = [-0.02,0.1]*h.time*units.lengthunit
print,xminmax
print,yminmax
print,zminmax
window,2,xsize = 1200,ysize = 600
multiplot,[2,1],/square
plot,s.x,s.y,psym = 3;,xrange = xminmax,yrange = yminmax,xstyle = 1,ystyle = 1,/nodata
oplot,s.x,s.y,psym = 3,color = 190
oplot,g[ind1_tip].x,g[ind1_tip].y,psym = 3,color = colors[0]
vir = circle(amigastat[where(amigastat.group EQ haloids[0])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[0])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[0])].rvir)
oplot,vir[0,*],vir[1,*]
IF n_elements(haloids) GE 2 THEN BEGIN & $
  oplot,g[ind2_tip].x,g[ind2_tip].y,psym = 3,color = colors[1] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[1])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[1])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[1])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 3 THEN BEGIN & $
  oplot,g[ind3_tip].x,g[ind3_tip].y,psym = 3,color = colors[2] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[2])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[2])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[2])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 4 THEN BEGIN & $
  oplot,g[ind4_tip].x,g[ind4_tip].y,psym = 3,color = colors[3]  & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[3])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[3])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[3])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 5 THEN BEGIN & $
  IF n_elements(ind5_tip) GT 1 THEN oplot,g[ind5_tip].x,g[ind5_tip].y,psym = 3,color = colors[4] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[4])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[4])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[4])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 6 THEN BEGIN & $
  IF n_elements(ind6_tip) GT 1 THEN oplot,g[ind6_tip].x,g[ind6_tip].y,psym = 3,color = colors[5] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[5])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[5])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[5])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 7 THEN BEGIN & $
  IF n_elements(ind7_tip) GT 1 THEN oplot,g[ind7_tip].x,g[ind7_tip].y,psym = 3,color = colors[6] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[6])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[6])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[6])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 8 THEN BEGIN & $
  IF n_elements(ind8_tip) GT 1 THEN oplot,g[ind8_tip].x,g[ind8_tip].y,psym = 3,color = colors[7] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[7])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[7])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[7])].rvir) & $
oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 9 THEN BEGIN & $
  IF n_elements(ind9_tip) GT 1 THEN oplot,g[ind9_tip].x,g[ind9_tip].y,psym = 3,color = colors[8] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[8])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[8])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[8])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 10 THEN BEGIN & $
  IF n_elements(ind10.skp_tip) GT 1 THEN oplot,g[ind10_tip].x,g[ind10_tip].y,psym = 3,color = colors[9] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[9])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[9])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[9])].rvir) & $
oplot,vir[0,*],vir[1,*] & $
ENDIF

oplot,s.x,s.y,psym = 3,color = 190
multiplot
plot,s.x,s.z,psym = 3;,xrange = xminmax,yrange = zminmax,xstyle = 1,ystyle = 1,/nodata
oplot,s.x,s.z,psym = 3,color = 190
oplot,g[ind1_tip].x,g[ind1_tip].z,psym = 3,color = colors[0]
vir = circle(amigastat[where(amigastat.group EQ haloids[0])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[0])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[0])].rvir)
oplot,vir[0,*],vir[1,*]
IF n_elements(haloids) GE 2 THEN BEGIN & $
  oplot,g[ind2_tip].x,g[ind2_tip].z,psym = 3,color = colors[1] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[1])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[1])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[1])].rvir) & $
  oplot,vir[0,*],vir[1,*],color = colors[1] & $
  ENDIF
IF n_elements(haloids) GE 3 THEN BEGIN & $
  oplot,g[ind3_tip].x,g[ind3_tip].z,psym = 3,color = colors[2] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[2])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[2])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[2])].rvir) & $
  oplot,vir[0,*],vir[1,*],color = colors[2] & $
  ENDIF
IF n_elements(haloids) GE 4 THEN BEGIN & $
  oplot,g[ind4_tip].x,g[ind4_tip].z,psym = 3,color = colors[3] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[3])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[3])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[3])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 5 THEN BEGIN & $
  oplot,g[ind5_tip].x,g[ind5_tip].z,psym = 3,color = colors[4] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[4])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[4])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[4])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 6 THEN BEGIN & $
  oplot,g[ind6_tip].x,g[ind6_tip].z,psym = 3,color = colors[5] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[5])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[5])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[5])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 7 THEN BEGIN & $
  oplot,g[ind7_tip].x,g[ind7_tip].z,psym = 3,color = colors[6] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[6])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[6])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[6])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 8 THEN BEGIN & $
  oplot,g[ind8_tip].x,g[ind8_tip].z,psym = 3,color = colors[7] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[7])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[7])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[7])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
ENDIF
IF n_elements(haloids) GE 9 THEN BEGIN & $
  oplot,g[ind9_tip].x,g[ind9_tip].z,psym = 3,color = colors[8] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[8])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[8])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[8])].rvir) & $
  oplot,vir[0,*],vir[1,*] & $
  ENDIF
IF n_elements(haloids) GE 10 THEN BEGIN & $
  oplot,g[ind10_tip].x,g[ind10_tip].z,psym = 3,color = colors[9] & $
  vir = circle(amigastat[where(amigastat.group EQ haloids[9])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[9])].zc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[9])].rvir) & $
oplot,vir[0,*],vir[1,*] & $
ENDIF

oplot,s.x,s.z,psym = 3,color = 190
multiplot,/reset

stop
data1 = mrdfits('grp'+haloids[0]+'.allgas.entropy.fits',1)
match,allgas2,allgas1,temp1,temp2
data12 = data1[temp2]
FOR is = n_steps - min([n_elements(pasthalos1),n_elements(pasthalos2)]), n_steps - 1 DO BEGIN &  $
    IF (where(data12[*].grp[is] EQ pasthalos1[n_elements(pasthalos1) - 1 - is]))[0] NE -1 THEN data12[where(data12[*].grp[is] EQ pasthalos1[n_elements(pasthalos1) - 1 - is])].grp[is] = long(haloids[0]) &  $
    IF (where(data12[*].grp[is] EQ pasthalos2[n_elements(pasthalos6) - 1 - is]))[0] NE -1 THEN data12[where(data12[*].grp[is] EQ pasthalos2[n_elements(pasthalos2) - 1 - is])].grp[is] = long(haloids[1]) &  $
    IF (where((data12[*].grp[is] NE long(haloids[0])) AND $
              (data12[*].grp[is] NE long(haloids[1])) AND $
              (data12[*].grp[is] NE 0)))[0] NE -1 THEN $
              data12[where( $
              data12[*].grp[is] NE long(haloids[0]) AND $
              data12[*].grp[is] NE long(haloids[1]) AND $
              data12[*].grp[is] NE 0)].grp[is] = -1 &  $
ENDFOR
match,allgas5,allgas1,temp1,temp2
data15 = data1[temp2]
FOR is = n_steps - min([n_elements(pasthalos1),n_elements(pasthalos5)]), n_steps - 1 DO BEGIN &  $
    IF (where(data15[*].grp[is] EQ pasthalos1[n_elements(pasthalos1) - 1 - is]))[0] NE -1 THEN data15[where(data15[*].grp[is] EQ pasthalos1[n_elements(pasthalos1) - 1 - is])].grp[is] = long(haloids[0]) &  $
    IF (where(data15[*].grp[is] EQ pasthalos5[n_elements(pasthalos5) - 1 - is]))[0] NE -1 THEN data15[where(data15[*].grp[is] EQ pasthalos5[n_elements(pasthalos5) - 1 - is])].grp[is] = long(haloids[4]) &  $
    IF (where((data15[*].grp[is] NE long(haloids[0])) AND $
              (data15[*].grp[is] NE long(haloids[4])) AND $
              (data15[*].grp[is] NE 0)))[0] NE -1 THEN $
              data15[where( $
              data15[*].grp[is] NE long(haloids[0]) AND $
              data15[*].grp[is] NE long(haloids[4]) AND $
              data15[*].grp[is] NE 0)].grp[is] = -1 &  $
ENDFOR
match,allgas6,allgas1,temp1,temp2
data16 = data1[temp2]
FOR is = n_steps - min([n_elements(pasthalos1),n_elements(pasthalos6)]), n_steps - 1 DO BEGIN &  $
    IF (where(data16[*].grp[is] EQ pasthalos1[n_elements(pasthalos1) - 1 - is]))[0] NE -1 THEN data16[where(data16[*].grp[is] EQ pasthalos1[n_elements(pasthalos1) - 1 - is])].grp[is] = long(haloids[0]) &  $
    IF (where(data16[*].grp[is] EQ pasthalos6[n_elements(pasthalos6) - 1 - is]))[0] NE -1 THEN data16[where(data16[*].grp[is] EQ pasthalos6[n_elements(pasthalos6) - 1 - is])].grp[is] = long(haloids[5]) &  $
    IF (where((data16[*].grp[is] NE long(haloids[0])) AND $
              (data16[*].grp[is] NE long(haloids[5])) AND $
              (data16[*].grp[is] NE 0)))[0] NE -1 THEN $
              data16[where( $
              data16[*].grp[is] NE long(haloids[0]) AND $
              data16[*].grp[is] NE long(haloids[5]) AND $
              data16[*].grp[is] NE 0)].grp[is] = -1 &  $
ENDFOR

data12colors = data12.grp
data12colors[where(data12.grp EQ -1)] = 60 ;Different halo
data12colors[where(data12.grp EQ 0)] = 20 ;Not in halo
data12colors[where(data12.grp EQ long(haloids[0]))] = 254 ;halo 1
data12colors[where(data12.grp EQ long(haloids[1]))] = 200 ;halo 2

data15colors = data15.grp
data15colors[where(data15.grp EQ -1)] = 60 ;Different halo
data15colors[where(data15.grp EQ 0)] = 20 ;Not in halo
data15colors[where(data15.grp EQ long(haloids[0]))] = 254 ;halo 1
data15colors[where(data15.grp EQ long(haloids[4]))] = 200 ;halo 2

data16colors = data16.grp
data16colors[where(data16.grp EQ -1)] = 60 ;Different halo
data16colors[where(data16.grp EQ 0)] = 20 ;Not in halo
data16colors[where(data16.grp EQ long(haloids[0]))] = 254 ;halo 1
data16colors[where(data16.grp EQ long(haloids[5]))] = 200 ;halo 2


window,2
plot,[.1,.1],[100,100],/xlog,/ylog,xrange = [1e-6,1e3],yrange = [10,1e7],xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]',/nodata
FOR ip = 0, n_elements(data12) - 1 DO oplot,data12[ip].rho,data12[ip].temp,color = ip
;FOR ip = 0, n_elements(data12) - 1 DO BEGIN & $
;  oplot,data12[ip].rho,data12[ip].temp,color = ip & $
;  wait,0.5 & $
;ENDFOR

plot,[0,0],[0,0],xrange = minmax(data12.x),yrange = minmax(data12.y)
FOR ip = 0, n_elements(data12) - 1 DO oplot,data12[ip].x,data12[ip].y,color = ip
oplot,(data12.x)[n_elements(data12[0].x) - 1,*],(data12.y)[n_elements(data12[0].x) - 1,*],psym = 4

plot,[.1,.1],[100,100],/xlog,/ylog,xrange = [1e-6,1e3],yrange = [10,1e7],xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]',/nodata
FOR ip = 0, n_elements(data16) - 1 DO oplot,data16[ip].rho,data16[ip].temp,color = ip
stop

plot,[0,0],[0,0],xrange = minmax(data12.x),yrange = minmax(data12.y)
FOR ip = 0, n_elements(data12) - 1 DO $ ;oplot,data12[ip].x,data12[ip].y,color = ip
  FOR is = 0, n_steps -2 DO $
    oplot,[data12[ip].x[is],data12[ip].x[is+1]],[data12[ip].y[is],data12[ip].y[is+1]],color = data12colors[is,ip] ;ip
oplot,(data12.x)[n_elements(data12[0].x) - 1,*],(data12.y)[n_elements(data12[0].x) - 1,*],psym = 4
vir = circle(amigastat[where(amigastat.group EQ haloids[1])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[1])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[1])].rvir)
oplot,vir[0,*],vir[1,*]
vir = circle(amigastat[where(amigastat.group EQ haloids[0])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[0])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[0])].rvir)
oplot,vir[0,*],vir[1,*]

plot,[0,0],[0,0],xrange = minmax(data15.x),yrange = minmax(data15.y)
FOR ip = 0, n_elements(data15) - 1 DO $ ;oplot,data15[ip].x,data15[ip].y,color = ip
  FOR is = 0, n_steps -2 DO $
    oplot,[data15[ip].x[is],data15[ip].x[is+1]],[data15[ip].y[is],data15[ip].y[is+1]],color = data15colors[is,ip] ;ip
oplot,(data15.x)[n_elements(data15[0].x) - 1,*],(data15.y)[n_elements(data15[0].x) - 1,*],psym = 4
vir = circle(amigastat[where(amigastat.group EQ haloids[4])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[4])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[4])].rvir)
oplot,vir[0,*],vir[1,*]
vir = circle(amigastat[where(amigastat.group EQ haloids[0])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[0])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[0])].rvir)
oplot,vir[0,*],vir[1,*]

plot,[0,0],[0,0],xrange = minmax(data16.x),yrange = minmax(data16.y)
FOR ip = n_elements(data16) - 1, 0, -1 DO $ ;oplot,data16[ip].x,data16[ip].y,color = ip
  FOR is = 0, n_steps -2 DO $
    oplot,[data16[ip].x[is],data16[ip].x[is+1]],[data16[ip].y[is],data16[ip].y[is+1]],color = data16colors[is,ip] ;ip
oplot,(data16.x)[n_elements(data16[0].x) - 1,*],(data16.y)[n_elements(data16[0].x) - 1,*],psym = 4
vir = circle(amigastat[where(amigastat.group EQ haloids[5])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[5])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[5])].rvir)
oplot,vir[0,*],vir[1,*]
vir = circle(amigastat[where(amigastat.group EQ haloids[0])].xc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[0])].yc*1000-units.lengthunit/2,amigastat[where(amigastat.group EQ haloids[0])].rvir)
oplot,vir[0,*],vir[1,*]
stop

plot,[0,0],[0,0],xrange = [0,max(time)],yrange = [1e-6,1e4],/ylog,xtitle = 'Time [Gyr]',ytitle = 'Density [amu/cc]'
FOR ip = 0, n_elements(data12) - 1 DO $
  FOR is = 0, n_steps -2 DO $
    oplot,[time[is],time[is+1]],[data12[ip].rho[is],data12[ip].rho[is+1]],color = data12colors[is,ip] ;ip


plot,[0,0],[0,0],xrange = [0,max(time)],yrange = [1e-6,1e4],/ylog,xtitle = 'Time [Gyr]',ytitle = 'Density [amu/cc]'
FOR ip = 0, n_elements(data15) - 1 DO $
  FOR is = 0, n_steps -2 DO $
    oplot,[time[is],time[is+1]],[data15[ip].rho[is],data15[ip].rho[is+1]],color = data15colors[is,ip] ;ip
stop

plot,[0,0],[0,0],xrange = [0,max(time)],yrange = [1e-6,1e4],/ylog,xtitle = 'Time [Gyr]',ytitle = 'Density [amu/cc]'
FOR ip = n_elements(data16) - 1, 0, -1 DO $
  FOR is = 0, n_steps -2 DO $
    oplot,[time[is],time[is+1]],[data16[ip].rho[is],data16[ip].rho[is+1]],color = data16colors[is,ip] ;ip
stop

END
