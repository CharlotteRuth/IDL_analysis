
PRO master_master
dir='/nobackupp8/crchrist/MolecH/cptmarvel.cosmo25cmb.4096g/cptmarvel.cosmo25cmb.4096g5HbwK1BH/'
filebase='cptmarvel.cosmo25cmb.4096g5HbwK1BH'
finalids = ['1','2','4','5','6','7','10','11','13','14']

dir='/nobackupp8/crchrist/MolecH/rogue.cosmo25cmb.4096g/rogue.cosmo25cmb.4096g5HbwK1BH/'
filebase='rogue.cosmo25cmb.4096g5HbwK1BH'
finalids = ['1','3','7','8','10','11','12','16','17','18','30','31','32','34']
finalids = ['31'] ;Still to do gas_gal_disk_history

dir = '/nobackupp8/crchrist/MolecH/elektra.cosmo25cmb.4096g/elektra.cosmo25cmb.4096g5HbwK1BH/'
filebase='elektra.cosmo25cmb.4096g5HbwK1BH'
finalids=['1','2','3','4','5','8','9','10','11','12','17','18'] ;still to do star accretion
finalids = ['18'] ;Still to do gas_gal_disk_history

dir = '/nobackupp8/crchrist/MolecH/storm.cosmo25cmb.4096g/storm.cosmo25cmb.4096g5HbwK1BH/' 
filebase='storm.cosmo25cmb.4096g5HbwK1BH' 
;finalids= ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60'] 
finalids = ['1','2','3','4','5','6','7','8','10','11','15','16','17','23','24','28','34','35','43','49','50','60'] ;still to do star accretion
finalids = ['15'] ;Still to do gas_gal_disk_history
finalids=['8'] ;Still to do with star accretion history

dir = '/nobackupp8/crchrist/MolecH/h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/'
filebase = 'h242.cosmo50PLK.3072gst5HbwK1BH'
finalids = ['1','9','11','24','29','33','39','40','45','72','75','76','425','457']

dir = '/nobackupp8/crchrist/MolecH/h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/' 
filebase='h229.cosmo50PLK.3072gst5HbwK1BH' 
finalids= ['1','2','4','7','17','21','22','27','51','52','70','104','203']

dir = '/nobackupp8/crchrist/MolecH/h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/'
filebase = 'h329.cosmo50PLK.3072gst5HbwK1BH'
finalids = ['1','9','32','126','129']

master_gaspass,dir,filebase,finalids

END

PRO master_gaspass, dir, filebase, finalids 
cd,dir
;start by having amiga run
;ahf_grp_stat, tipsyfile, boxsize=18.24873, munit = 2.310e15, vunit = 630.48, h0 = 0.73, multiplefile = 'ahf_grp_stat_input.txt'

print,'find_star_history_master'
find_star_history_master, dir, filebase, finalid = finalids,/rewrite;/u/crchrist/IDL/IDL_analysis/PartTrace/find_star_history_master.pro
;use haloid and mergertree to get correct merger tree (*halo_step.out
;and *haloid.dat)
IF 1 THEN BEGIN
   print,'haloid'
   FOR i=0, n_elements(finalids)-1 DO haloid, filebase, finalid = finalids[i]
   print,'mergertree'
   FOR i=0, n_elements(finalids)-1 DO mergertree, filebase, dir = dir, gal = finalids[i],/plot ;/u/crchrist/IDL/IDL_analysis/Twins/mergertree.pro
;/home1/crchrist/IDL/Twins/mergertree.pro
   FOR i=1, n_elements(finalids)-1 DO edit_haloid, filebase, finalids[i]
   print,'Stop and compare merger trees to haloid file'
   stop
   FOR i=0, n_elements(finalids)-1 DO mergertree_check,filebase, finalid = finalids[i]
   stop
ENDIF
print,'find_all_gas'
FOR i=0, n_elements(finalids)-1 DO find_all_gas, filebase, finalid = finalids[i]
print,'find_gas_history_master'
;FOR i=0, n_elements(finalids)-1 DO find_gas_history_master, filebase, finalid = finalids[i]
find_gas_history_master, filebase, finalid = finalids
print,'accrmode'
FOR i=0, n_elements(finalids)-1 DO accrmode, filebase, finalid = finalids[i]
print,'disk_alignment'
;.r /home6/crchrist/IDL/IDL_analysis/idl4tipsy/align.pro
FOR i=0, n_elements(finalids)-1 DO disk_alignment,filebase, finalid = finalids[i]
;print,'disk_alignment2' ;incorporated into disk_alignment
;FOR i=0, n_elements(finalids)-1 DO disk_alignment2,finalid = finalids[i];,/debug
;print,'make_halo_starlog'
;FOR i=0, n_elements(finalids)-1 DO make_halo_starlog,filebase,finalid = finalids[i],/molecularH
;print,'cutstarlog'
;FOR i=0, n_elements(finalids)-1 DO cutstarlog,filebase,finalid = finalids[i],/molecularH
FOR i=0, n_elements(finalids)-1 DO star_accretion,filebase,finalids[i],/molecularH,laststep = '4096',/big
;------------------------- Gas History ----------------------
print,'gas_gal_disk_history'
FOR i=0, n_elements(finalids)-1 DO gas_gal_disk_history,laststep = '4096',finalid = finalids[i];,/plots,/debug,/steps_debug;
FOR i=0, n_elements(finalids)-1 DO gas_gal_disk_history,finalid = finalids[i],/cooloncut
FOR i=0, n_elements(finalids)-1 DO gas_gal_disk_history,finalid = finalids[i],/rvircut_fifth
FOR i=0, n_elements(finalids)-1 DO gas_gal_disk_history,finalid = finalids[i],/rvircut_half



print,'metal_history'
FOR i=0, n_elements(finalids)-1 DO metal_history,finalid = finalids[i]
print,'track_mass_2'
FOR i=0, n_elements(finalids)-1 DO track_mass_2,finalid = finalids[i]
print,'reeject_gas_character, expell'
FOR i=0, n_elements(finalids)-1 DO reeject_gas_character,ejecthistory,finalid = finalids[i],/expell
FOR i=0, n_elements(finalids)-1 DO reeject_gas_character,ejecthistory,finalid = finalids[i],/rvircut_half
FOR i=0, n_elements(finalids)-1 DO reeject_gas_character,ejecthistory,finalid = finalids[i],/rvircut_fifth
FOR i=0, n_elements(finalids)-1 DO reeject_gas_character,ejecthistory,finalid = finalids[i],/heat
print,'reeject_gas_character, postdisk'
FOR i=0, n_elements(finalids)-1 DO reeject_gas_character,ejecthistory,finalid = finalids[i],/postdisk
print,'reaccr_gas_character'
FOR i=0, n_elements(finalids)-1 DO reaccr_gas_character,accrhistory,finalid = finalids[i];,/heat
FOR i=0, n_elements(finalids)-1 DO time_cycling_exp,[dir],[filebase],halo = [finalids[i]]
FOR i=0, n_elements(finalids)-1 DO time_cycling,[dir],[filebase],halo = [finalids[i]]
FOR i=0, n_elements(finalids)-1 DO track_vc,finalid = finalids[i]
;---------------------------------------------
IF file_test(dir + '/' + filebase + '.' + '00080' + '/' + filebase + '.' + '00080') $
THEN firststep = '00080' $
ELSE firststep = '00084'      
FOR i=0, n_elements(finalids)-1 DO earlygas,dir,filebase,firststep,finalids[i]
END
