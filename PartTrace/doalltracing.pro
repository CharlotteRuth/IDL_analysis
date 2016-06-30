;file = 'h603.cosmo50cmb.3072g14HBWK'
;dir = '/nobackupp2/crchrist/MolecH/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK'
;haloids = ['1','2','3']

PRO doalltracing,dir,file,haloids
cd,dir
FOR i = 0, n_elements(haloids) - 1 DO BEGIN
   haloid = haloids[i]
   print,haloid
   print,'find_star_history'
   find_star_history_master,dir,file,haloid
   print,'find_gas_history'
   find_gas_history_master,file,haloid
   print,'disk_alignment'
   disk_alignment,file,haloid = haloid
   disk_alignment2,haloid
   print,'accrmode'
   accrmode,file,haloid
   print,'gas_gal_disk'
   gas_gal_disk_history,dir = dir,finalid = haloid
   print,'reeject_gas_character'
   reeject_gas_character,dir,ejecthistory,expellhistory,finalid = haloid
   reeject_gas_character,dir,ejecthistory,expellhistory,finalid = haloid,/postdisk
ENDFOR

END

