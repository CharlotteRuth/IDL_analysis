
PRO master_master
dir='/nobackupp8/crchrist/MolecH/cptmarvel.cosmo25.4096g/cptmarvel.cosmo25cmb.4096g5HbwK1BH/'
filebase='cptmarvel.cosmo25cmb.4096g5HbwK1BH'
finalids = ['1','2','4','5','6','7','10','11','13','14']

master_gaspass,dir,filebase,finalids
END

PRO master_gaspass, dir, filebase, finalids 
cd,dir
;start by having amiga run
;ahf_grp_stat, tipsyfile, boxsize=18.24873, munit = 2.310e15, vunit = 630.48, h0 = 0.73, multiplefile = 'ahf_grp_stat_input.txt'

print,'find_star_history_master'
FOR i=0, n_elements(finalids)-1 DO find_star_history_master, dir, filebase, finalid = finalids[i];/u/crchrist/IDL/IDL_analysis/PartTrace/find_star_history_master.pro
;use haloid and mergertree to get correct merger tree
IF 1 THEN BEGIN
   print,'haloid'
   FOR i=0, n_elements(finalids)-1 DO haloid, filebase, finalid = finalids[i]
   print,'mergertree'
   FOR i=0, n_elements(finalids)-1 DO mergertree, filebase, dir = dir, gal = finalids[i],/plot ;/u/crchrist/IDL/IDL_analysis/Twins/mergertree.pro
;/home1/crchrist/IDL/Twins/mergertree.pro
   print,'Stop and compare merger trees to haloid file'
   stop
   FOR i=0, n_elements(finalids)-1 DO mergertree_check,filebase, finalid = finalids[i]
   stop
ENDIF
print,'find_all_gas'
FOR i=0, n_elements(finalids)-1 DO find_all_gas, filebase, finalid = finalids[i]
print,'find_gas_history_master'
FOR i=0, n_elements(finalids)-1 DO find_gas_history_master, filebase, finalid = finalids[i]
print,'accrmode'
FOR i=0, n_elements(finalids)-1 DO accrmode, filebase, finalid = finalids[i]
print,'disk_alignment'
FOR i=0, n_elements(finalids)-1 DO disk_alignment,filebase, finalid = finalids[i]
print,'disk_alignment2'
FOR i=0, n_elements(finalids)-1 DO disk_alignment2,finalid = finalids[i];,/debug
print,'metal_history'
FOR i=0, n_elements(finalids)-1 DO metal_history,finalid = finalids[i]
print,'track_mass_2'
FOR i=0, n_elements(finalids)-1 DO track_mass_2,finalid = finalids[i]
print,'make_halo_starlog'
FOR i=0, n_elements(finalids)-1 DO make_halo_starlog,filebase,finalid = finalids[i],/molecularH
print,'cutstarlog'
FOR i=0, n_elements(finalids)-1 DO cutstarlog,filebase,finalid = finalids[i],/molecularH
;------------------------- Gas History ----------------------
print,'gas_gal_disk_history'
FOR i=0, n_elements(finalids)-1 DO gas_gal_disk_history,finalid = finalids[i];,/plots,/debug,/steps_debug,/all_eject
FOR i=0, n_elements(finalids)-1 DO gas_gal_disk_history,finalid = finalids[i],/cooloncut
FOR i=0, n_elements(finalids)-1 DO gas_gal_disk_history,finalid = finalids[i],/rvircut_fifth
FOR i=0, n_elements(finalids)-1 DO gas_gal_disk_history,finalid = finalids[i],/rvircut_half
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
FOR i=0, n_elements(finalids)-1 DO star_accretion,filebase,finalids[i],/molecularH
END
