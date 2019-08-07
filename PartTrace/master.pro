
PRO master_master
dir='/nobackupp8/crchrist/MolecH/h2003.cosmo25cmb.4096g/h2003.cosmo25cmb.4096g14HBK/'
filebase='h2003.cosmo25cmb.4096g14HBK'
finalids = ['1']

dir='/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
filebase='h799.cosmo25cmb.3072g14HBWK'
;finalids = ['1']
finalids = ['1','4','6','10']
finalids = ['1','4','6']
;step = '00076'

dir='/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
filebase='h516.cosmo25cmb.3072g14HBWK'
finalids = ['1','2']
;finalids = ['3','6']
;finalids = ['1','2','3','6']

dir='/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
filebase='h986.cosmo50cmb.3072g14HBWK'
finalids = ['1','2','3','8']
finalids = ['15','16']
finalids = ['1','2','3','8','15','16']
finalids = ['11','27'] ;Done for Anna

;dir='/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK_short/'
;filebase='h986.cosmo50cmb.3072g14HBWK'
;finalids = ['1']

;dir='/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/'
;filebase='h986.cosmo50cmb.3072g1bwK'
;finalids = ['1']

;dir='/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/'
;filebase='h986.cosmo50cmb.3072gs1MbwK'
;finalids = ['1']

dir='/nobackupp8/crchrist/MolecH/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
filebase='h603.cosmo50cmb.3072g14HBWK'
finalids = ['1','2','3']
;finalids = ['7','13']

dir='/nobackupp8/crchrist/MolecH/h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
filebase='h285.cosmo50cmb.3072g14HMbwK'
;finalids = ['1','4','8']
finalids = ['1','4','9']

dir='/nobackupp8/crchrist/MolecH/h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
filebase='h258.cosmo50cmb.3072g14HMbwK'
finalids = ['1','4','8']

dir='/nobackupp8/crchrist/MolecH/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
filebase='h239.cosmo50cmb.3072g14HMbwK'
finalids=['1']

;dir='/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3BWK/'
;filebase='h799.cosmo25cmb.3072g3BWK'
;finalids = ['1']

;dir='/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g1MBWK/'
;filebase='h799.cosmo25cmb.3072g1MBWK'
;finalids = ['1']

;dir='/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK_short/'
;filebase='h799.cosmo25cmb.3072g14HBWK'
;finalids = ['1']
END

PRO master, dir, filebase, finalids 
cd,dir
;start by having amiga run
print,'find_star_history_master'
FOR i=0, n_elements(finalids)-1 DO find_star_history_master, dir, filebase, finalid = finalids[i]
;use haloid and mergertree get correct merger tree
IF 1 THEN BEGIN
   print,'haloid'
   FOR i=0, n_elements(finalids)-1 DO haloid, filebase, finalid = finalids[i]
   print,'mergertree'
   FOR i=0, n_elements(finalids)-1 DO mergertree, filebase, dir = dir, gal = finalids[i],/plot ;/home1/crchrist/IDL/Twins/mergertree.pro
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
print,'metal_history'x
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
