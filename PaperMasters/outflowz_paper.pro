;outplot = '/christensen/Plots/fbmass/hgal'
;outplot = '/home/christensen/Plots/fbmass/hgal_thick'
;outplot = '/home6/crchrist/Plots/hgal'
;outplot = '/home6/crchrist/Plots/hgal_thick'

PRO outflowz_paper,outplot = outplot,color = color,verbose = verbose,formatthick = formatthick,redo = redo
formatplot,outplot = outplot,thick = formatthick
IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(outplot) THEN bgcolor = 255 ELSE bgcolor = 0

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' $
ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) OR strcmp(hostname, 'ldan',4) THEN prefix = '/nobackupp8/crchrist/MolecH/' $
ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
prefix = '/nobackupp8/crchrist/MolecH/'
zsolar = 0.0130215

pre_get_halo_mags = 0           ;pre
pre_cuthalos = 0                ;pre
pre_make_halo_starlog = 0       ;pre
pre_reaccr_gas_character = 0    ;pre ;in Master
pre_event_gas_character = 0
pre_metal = 0                   ;run code to calculate metals in halo and disk at each step. Outputs *metals.txt
run_times_cycling = 0           ;pre ;in Master
plot_vcirc = 0                  ;pre
do_writeHImacro = 0             ;pre
do_line_profile = 0             ;pre
check_contam = 0
pre_earlygas = 0                 ;in Master
pre_stellarz = 0   ;Compute metals produced by stars and in stars
pre_z_at_ejection = 0
pre_mzr = 0

plot_multiHaloGuo = 0
plot_baryonicfrac = 0
plot_tully_fisher_obs = 0
plot_baryonic_tully_fisher_obs = 0
plot_galMS = 0
plot_mzr = 0
plot_mzr_stellar = 0

find_reaccr = 0
plot_track_mass =0
plot_fbcum = 0

plot_ejectz_v_mass = 0 ;1 Metal Mass loading and fraction of metals ejected & expelled
plot_metals_v_mass = 1 ;1 Fraction of metals produced that are retained
plot_times_ejected = 0 ;
plot_times_expelled = 0
plot_times_cycling = 0 ;
plot_reeject_r = 0 ;
plot_reaccr_r = 0
plot_angmom = 0 ;
plot_reeject_z = 0 ;1 Metallicity of Ejecta
plot_reeject_v = 0;
plot_sfh_reaccr = 0
plot_outflowr = 0

plot_inflow_outflow_historyz = 0 ;1 History of metal enrichment
plot_inflow_outflow_metallicity = 0 ;History of inflow and outflow metallicity. Relies on output from plot_inflow_outflow_historyz

check_coolontime = 0

mu_c50 = 1.84793e16
mu_c25 = 2.310e15
lu_c50 = 50000.
lu_c25 = 25000.

zstep = ['00084','00112','00260','00512'];116
;zstep = ['00084','00100','00112','00124','00136','00148','00160','00176','00184','00196','00208','00220','00232','00244','00260','00268','00280','00292','00304','00316','00328','00340','00352','00364','00376','00388','00404','00412','00424','00436','00448','00460','00472','00484','00496','00512']
;zstep = ['00512']

x = 7
CASE x OF
    1: BEGIN
        dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
        file799 = 'h799.cosmo25cmb.3072g14HBWK'
        key799 = 'h799'
        dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
        file516 = 'h516.cosmo25cmb.3072g14HBWK'
        key516 = 'h516'
        dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
        file986 = 'h986.cosmo50cmb.3072g14HBWK'
        key986 = 'h986'
        dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
        file603 = 'h603.cosmo50cmb.3072g14HBWK'
        key603 = 'h603'
        dir277 = prefix + 'h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/'
        file277 = 'h277.cosmo50cmb.3072g14HMbwK'
        key277 = 'h277'
        dirs    = [dir799,dir516,dir516,dir986,dir986,dir986,dir986,dir603,dir603,dir603,dir277,dir277]
        files    = [file799,file516,file516,file986,file986,file986,file986,file603,file603,file603,file277,file277]
        finalstep = '00512'
        outfiles = dirs + files + '.' + finalstep + '/' + files + '.' + finalstep
        haloid = ['1'   ,'1'   ,'2'   ,'1'   ,'2'   ,'3'   ,'9'   ,'1'   ,'2'    ,'3'  ,'1'   ,'2']
        outhalos = outfiles + '.halo.' + haloid 
        distunits = [fltarr(3) + lu_c25, fltarr(9) + lu_c50]
        massunits = [fltarr(3) + mu_c25, fltarr(9) + mu_c50]
        key  = [key799,key516,key516,key986,key986,key986,key986,key603,key603,key603,key277,key277] + ', ' + haloid
;        dirs    = [dir799,dir516,dir516,dir986,dir986,dir986,dir986,dir603,dir603,dir603,dir277]
;        files    = [file799,file516,file516,file986,file986,file986,file986,file603,file603,file603,file277]
;        haloid = ['1'   ,'1'   ,'2'   ,'1'   ,'2'   ,'3'   ,'9'   ,'1',   '2'    ,'3'  ,'2']

;        colors = fltarr(n_elements(files)) + fgcolor
;        loadct = fltarr(n_elements(files))
;        linestyles = fltarr(n_elements(files))
        masssort = [6,2,0,9,1,5,4,11,8,3,7,10]
;        psym = fltarr(n_elements(files)) + 16
        psym = [fltarr(3) + 14, fltarr(9) + 4]
        obscolor = 100
        obssym = 2
    END
    2: BEGIN
        dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
        file799 = 'h799.cosmo25cmb.3072g14HBWK'
        key799 = 'h799'
        dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
        file516 = 'h516.cosmo25cmb.3072g14HBWK'
        key516 = 'h516'
        dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
        file986 = 'h986.cosmo50cmb.3072g14HBWK'
        key986 = 'h986'
        dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
        file603 = 'h603.cosmo50cmb.3072g14HBWK'
        key603 = 'h603'
        dir277 = prefix + 'h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/'
        file277 = 'h277.cosmo50cmb.3072g14HMbwK'
        key277 = 'h277'
        dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
        file258 = 'h258.cosmo50cmb.3072g14HMbwK'
        key258 = 'h258'
        dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
        file285 = 'h285.cosmo50cmb.3072g14HMbwK'
        key285 = 'h285'
        dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
        file239 = 'h239.cosmo50cmb.3072g14HMbwK'
        key239 = 'h239'
        dirs    = [dir799,dir799,dir799,dir516,dir516,dir986,dir986,dir986,dir986,dir603,dir603,dir603,dir277,dir277,dir258,dir258,dir285,dir285,dir285,dir239]
        files    = [file799,file799,file799,file516,file516,file986,file986,file986,file986,file603,file603,file603,file277,file277,file258,file258,file285,file285,file285,file239]
        finalstep = '00512'
        outfiles = dirs + files + '.' + finalstep + '/' + files + '.' + finalstep
        haloid = ['1'   ,'4'      ,'6'     ,'1'   ,'2'   ,'1'   ,'2'   ,'3'   ,'9'   ,'1'   ,'2'    ,'3'  ,'1'   ,'2'  ,'1'   ,'4'   ,'1'   ,'4'    ,'9'   ,'1']
        haloid_split = [['1','4','6','0'],['1','2','0','0'],['1','2','3','9'],['1','2','3','0'],['1','2','0','0'],['1','4','0','0'],['1','4','9','0'],['1','0','0','0']]
        outhalos = outfiles + '.halo.' + haloid 
        distunits = [fltarr(5) + lu_c25, fltarr(15) + lu_c50]
        massunits = [fltarr(5) + mu_c25, fltarr(15) + mu_c50]
        key  = [key799,key799,key799,key516,key516,key986,key986,key986,key986,key603,key603,key603,key277,key277,key258,key258,key285,key285,key285,key239] + ', ' + haloid
        masssort = [8,1,0,15,18,4,2,11,17,3,7,6,13,10,5,9,12,14,16,19]
;        dirs    = [dir799,dir516,dir516,dir986,dir986,dir986,dir986,dir603,dir603,dir603,dir277]
;        files    = [file799,file516,file516,file986,file986,file986,file986,file603,file603,file603,file277]
;        haloid = ['1'   ,'1'   ,'2'   ,'1'   ,'2'   ,'3'   ,'9'   ,'1',   '2'    ,'3'  ,'2']

;        colors = fltarr(n_elements(files)) + fgcolor
;        loadct = fltarr(n_elements(files))
;        linestyles = fltarr(n_elements(files))
        psym = fltarr(n_elements(files)) + 14;16
        psym = [fltarr(5) + 14, fltarr(15) + 4]
        obscolor = 100
        obssym = 2
        vfinal = [54.4300, 32.5, 27.5, 65.0886, 39.9453, 100.913, 76.4687, 70.9076, 38.4851, 111.279, 99.5493, 53.7934, 213.539, 72.8713, 208.034, 40.2975, 203.057, 58.6933, 43.0594, 204.288] ;Made using vcirc.pro
        vfinal = [55.1222, 32.5, 27.5, 66.9988, 33.8313, 103.148, 76.7641, 76.2822, 35.0817, 115.102, 74.9589, 50.1720, 188.999, 76.5042, 181.909, 43.2878, 164.350, 64.0178, 52.6147, 165.470];h799 4 and 6 need to be updated
    END
   3: BEGIN
        dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
        file799 = 'h799.cosmo25cmb.3072g14HBWK'
        key799 = 'h799'
        dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
        file516 = 'h516.cosmo25cmb.3072g14HBWK'
        key516 = 'h516'
        dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
        file986 = 'h986.cosmo50cmb.3072g14HBWK'
        key986 = 'h986'
        dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
        file603 = 'h603.cosmo50cmb.3072g14HBWK'
        key603 = 'h603'
        dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
        file258 = 'h258.cosmo50cmb.3072g14HMbwK'
        key258 = 'h258'
        dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
        file285 = 'h285.cosmo50cmb.3072g14HMbwK'
        key285 = 'h285'
        dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
        file239 = 'h239.cosmo50cmb.3072g14HMbwK'
        key239 = 'h239'
        dirs    = [dir799,dir799,dir799,dir516,dir516,dir986,dir986,dir986,dir986,dir603,dir603,dir603,dir258,dir258,dir285,dir285,dir285,dir239]
;        dirs    = [dir799,dir516,dir516,dir986,dir986,dir986,dir986,dir603,dir603,       dir258,dir258,dir285,dir285,dir285,dir239]
        files    = [file799,file799,file799,file516,file516,file986,file986,file986,file986,file603,file603,file603,file258,file258,file285,file285,file285,file239]
;        files    = [file799,file516,file516,file986,file986,file986,file986,file603,file603,        file258,file258,file285,file285,file285,file239]
        finalstep = '00512'
        outfiles = dirs + files + '.' + finalstep + '/' + files + '.' + finalstep
        haloid = ['1'   ,'4'   ,'6'   ,'1'   ,'2'   ,'1'   ,'2'   ,'3'   ,'8'   ,'1'   ,'2'    ,'3'  ,'1'   ,'4'   ,'1'   ,'4'    ,'9'   ,'1']
;        haloid = ['1'   ,'1'   ,'2'   ,'1'   ,'2'   ,'3'   ,'8'   ,'1'   ,'2'          ,'1'   ,'4'   ,'1'   ,'4'    ,'9'   ,'1']
        haloid_split = [['1','4','6','0'],['1','2','0','0'],['1','2','3','8'],['1','2','3','0'],['1','4','0','0'],['1','4','9','0'],['1','0','0','0']]
;        haloid_split = [['1','0','0','0'],['1','2','0','0'],['1','2','3','8'],['1','2','0','0'],['1','4','0','0'],['1','4','9','0'],['1','0','0','0']]
        outhalos = outfiles + '.halo.' + haloid 
        nhighres = 5
        nlowres = 13
        distunits = [fltarr(nhighres) + lu_c25, fltarr(nlowres) + lu_c50]
        massunits = [fltarr(nhighres) + mu_c25, fltarr(nlowres) + mu_c50]
        key  = [key799,key799,key799,key516,key516,key986,key986,key986,key986,key603,key603,key603,key258,key258,key285,key285,key285,key239] + ', ' + haloid
;        key  = [key799,key516,key516,key986,key986,key986,key986,key603,key603,key258,key258,key285,key285,key285,key239] + ', ' + haloid
        masssort = [8,1,0,13,16, 4, 2, 11,15, 7, 3, 6, 10, 5, 9,11,14,17]
;        masssort = [6 ,9  ,12 ,2 ,0 ,11 ,1 ,5 ,4 , 8 ,3 , 7 ,    10 ,13]
;789

;        colors = fltarr(n_elements(files)) + fgcolor
;        loadct = fltarr(n_elements(files))
;        linestyles = fltarr(n_elements(files))
        psym = fltarr(n_elements(files)) + 14;16
        psym = [fltarr(nhighres) + 14, fltarr(nlowres) + 6]
        obscolor = 100
        obssym = 2
        gmass = [2.51987e+08, 5.50022e+08, 4.55808e+07, 3.45968e+09, 7.36412e+08, 5.44062e+08, 1.27984e+07, 4.23651e+09, 7.89399e+08, 1.77581e+08, 5.76264e+09, 1.46396e+09, 5.68696e+09, 6.20778e+07, 8.47653e+09, 1.49208e+08, 1.26554e+08, 6.19468e+09]
;        vfinal = [54.4300, 65.0886, 39.9453, 100.913, 76.4687, 70.9076, 38.4851, 111.279, 99.5493, 53.7934, 208.034, 40.2975, 203.057, 58.6933, 43.0594, 204.288] ;Made using vcirc.pro
;        vfinal = [55.1222, 66.9988, 33.8313, 103.148, 76.7641, 76.2822, 35.0817, 115.102, 74.9589, 50.1720, 181.909, 43.2878, 164.350, 64.0178, 52.6147, 165.470]

        vfinal = [54.4300, 32.5, 27.5, 65.0886, 39.9453, 100.913, 76.4687, 70.9076, 38.4851, 111.279, 99.5493,          208.034, 40.2975, 203.057, 58.6933, 43.0594, 204.288] ;Made using vcirc.pro
        vfinal = [55.1222, 32.5, 27.5, 66.9988, 33.8313, 103.148, 76.7641, 76.2822, 35.0817, 115.102, 74.9589,50.1720,          181.909, 43.2878, 164.350, 64.0178, 52.6147, 165.470] ;h799 4 and 6 need to be updated

    END
    4: BEGIN
;        dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
;        file799 = 'h799.cosmo25cmb.3072g14HBWK'
;        key799 = 'h799'
        dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
        file258 = 'h258.cosmo50cmb.3072g14HMbwK'
        key258 = 'h258'
        dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
        file285 = 'h285.cosmo50cmb.3072g14HMbwK'
        key285 = 'h285'
        dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
        file239 = 'h239.cosmo50cmb.3072g14HMbwK'
        key239 = 'h239'
        dirs    = [dir258,dir258,dir285,dir285,dir285,dir239];[dir799,dir258,dir258,dir285,dir285,dir285,dir239]
        files    = [file258,file258,file285,file285,file285,file239];[file799,file258,file258,file285,file285,file285,file239]
        finalstep = '00512'
        outfiles = dirs + files + '.' + finalstep + '/' + files + '.' + finalstep
        haloid = ['1'   ,'4'   ,'1'   ,'4'    ,'9'   ,'1'];['1'   ,'1'   ,'4'   ,'1'   ,'4'    ,'9'   ,'1']
        haloid_split = [['1','4','0','0'],['1','4','9','0'],['1','0','0','0']];[['1','0','0','0'],['1','4','0','0'],['1','4','9','0'],['1','0','0','0']]
        outhalos = outfiles + '.halo.' + haloid 
        distunits = [fltarr(6) + lu_c50];[fltarr(1) + lu_c25, fltarr(6) + lu_c50]
        massunits = [fltarr(6) + mu_c50];[fltarr(1) + mu_c25, fltarr(6) + mu_c50]
        key  = [key258,key258,key285,key285,key285,key239] + ', ' + haloid;[key799,key258,key258,key285,key285,key285,key239] + ', ' + haloid

;        colors = fltarr(n_elements(files)) + fgcolor
;        loadct = fltarr(n_elements(files))
;        linestyles = fltarr(n_elements(files))
        psym = fltarr(n_elements(files)) + 16
        psym = [fltarr(6) + 4]
        obscolor = 100
        obssym = 2
    END
    5: BEGIN
        dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
        file799 = 'h799.cosmo25cmb.3072g14HBWK'
        key799 = 'h799'
        dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
        file986 = 'h986.cosmo50cmb.3072g14HBWK'
        key986 = 'h986'
        dirs    = [dir799,dir986,dir986,dir986,dir986]
        files    = [file799,file986,file986,file986,file986]
        finalstep = '00512'
        outfiles = dirs + files + '.' + finalstep + '/' + files + '.' + finalstep
        haloid = ['1'   ,'1'   ,'2'   ,'3'   ,'8']
        haloid_split = [['1','0','0','0'],['1','2','3','8']]
        outhalos = outfiles + '.halo.' + haloid 
        nhighres = 1
        nlowres = 4
        distunits = [fltarr(nhighres) + lu_c25, fltarr(nlowres) + lu_c50]
        massunits = [fltarr(nhighres) + mu_c25, fltarr(nlowres) + mu_c50]
        key  = [key799,key986,key986,key986,key986] + ', ' + haloid
        masssort = [1,2,3,0,4]

        psym = fltarr(n_elements(files)) + 5;
        psym = [fltarr(nhighres) + 14, fltarr(nlowres) + 4]
        obscolor = 100
        obssym = 2

;        vfinal = [54.4300, 65.0886, 39.9453, 100.913, 76.4687, 70.9076, 38.4851, 111.279, 99.5493, 53.7934, 208.034, 40.2975, 203.057, 58.6933, 43.0594, 204.288] ;Made using vcirc.pro
;        vfinal = [55.1222, 66.9988, 33.8313, 103.148, 76.7641, 76.2822, 35.0817, 115.102, 74.9589, 50.1720, 181.909, 43.2878, 164.350, 64.0178, 52.6147, 165.470]

        vfinal = [54.4300, 100.913, 76.4687, 70.9076, 38.4851] ;Made using vcirc.pro
        vfinal = [55.1222, 103.148, 76.7641, 76.2822, 35.0817]
    END
    6: BEGIN
        dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
        file799 = 'h799.cosmo25cmb.3072g14HBWK'
        key799 = 'h799'
        dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
        file516 = 'h516.cosmo25cmb.3072g14HBWK'
        key516 = 'h516'
        dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
        file986 = 'h986.cosmo50cmb.3072g14HBWK'
        key986 = 'h986'
        dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
        file603 = 'h603.cosmo50cmb.3072g14HBWK'
        key603 = 'h603'
        dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
        file258 = 'h258.cosmo50cmb.3072g14HMbwK'
        key258 = 'h258'
        dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
        file285 = 'h285.cosmo50cmb.3072g14HMbwK'
        key285 = 'h285'
        dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
        file239 = 'h239.cosmo50cmb.3072g14HMbwK'
        key239 = 'h239'
        dirs    = [ dir799, dir799, dir799, dir516, dir516, dir986, dir986, dir986, dir986, dir986, dir986, dir603, dir603, dir603, dir603, dir603, dir258, dir258, dir258, dir285, dir285, dir285, dir239, dir239, dir239]
        files   = [file799,file799,file799,file516,file516,file986,file986,file986,file986,file986,file986,file603,file603,file603,file603,file603,file258,file258,file258,file285,file285,file285,file239,file239,file239]
        finalstep = '00512'
        haloid =  [ '1'   , '4'   , '6'   , '1'   , '2'   , '1'   , '2'   , '3'   , '8'   , '15'  , '16'  , '1'   , '2'    , '3'  , '7'   , '13'  , '1'   , '4'   , '8'   , '1'   , '4'   , '9'   , '1'   , '13'  , '30'  ]
        haloid_split = [['1','4','6','0','0','0'],['1','2','0','0','0','0'],['1','2','3','8','15','16'],['1','2','3','7','13','0'],['1','4','8','0','0','0'],['1','4','9','0','0','0'],['1','13','30','0','0','0']]
        gmass   = [2.615e8,1.860e7,1.323e7,5.315e8,4.801e7,3.471e9,7.646e8,5.694e8,1.632e7,2.253e7,6.264e6,4.903e9,9.429e8,2.424e8,1.915e7,8.253e6,7.570e9,6.809e7,5.857e6,1.01e10,1.717e8,1.305e8,9.480e9,1.812e7,2.103e6]
 ;       vfinal =  [54.4300,    32.5,  27.5,65.0886,39.9453,100.913,76.4687,70.9076,38.4851, 0     ,  0    ,111.279, 99.5493, 0     , 0    , 0     ,208.034,40.2975, 0     ,203.057,58.6933,43.0594, 204.288] ;Made using vcirc.pro
 ;       vfinal =  [55.1222,    32.5,  27.5,66.9988,33.8313,103.148,76.7641,76.2822,35.0817, 0     ,  0    ,115.102, 74.9589,50.1720, 0    , 0     ,181.909,43.2878, 0     ,164.350,64.0178,52.6147,165.470] ;h799 4 and 6 need to be updated
        vfinal =  [55.1222,33.3827,26.5413,66.6838,42.7685,102.305,78.8046,76.3437,36.4566,28.7297,26.7203,115.102, 74.9589,50.1720,28.4053,26.7649,181.909,43.2878,36.0320,163.970,52.1639,54.3400,165.470,34.7140,20.9011] ;calculated using vcirc program.  Note that I adjusted the max fit radius of h603.cosmo50cmb.3072g14HBWK.00512.halo.7 to be 1.4115108 kpc to get a good fit
;
        outfiles = dirs + files + '.' + finalstep + '/' + files + '.' + finalstep
        outhalos = outfiles + '.halo.' + haloid 
        nhighres = 5
        nlowres = 20
        distunits = [fltarr(nhighres) + lu_c25, fltarr(nlowres) + lu_c50]
        massunits = [fltarr(nhighres) + mu_c25, fltarr(nlowres) + mu_c50]
        key  =    [ key799, key799, key799, key516, key516, key986, key986, key986, key986, key986, key986, key603, key603, key603, key603, key603, key258, key258, key258, key285, key285, key285, key239, key239, key239] + ', ' + haloid
        mvir =    [2.4e10,  6.8e9,  4.4e9,  3.8e10, 1.5e10, 1.9e11, 5.9e10, 3.8e10, 1.1e10, 4.4e9,  3.1e9,  3.4e11, 1.0e11, 2.9e10, 4.1e9,  1.5e9,  7.7e11, 1.1e10, 3.2e9,  8.8e11, 3.4e10, 1.2e10, 9.1e11, 2.6e9,  1.0e9] 
        masssort =[24,      15,     23,     10,     18,     14,     2,      9,      1,      8,      17,     21,     4,      0,      13,     20,     3,      7,      6,      12,     5,      11,     16,     19,     22]

;        colors = fltarr(n_elements(files)) + fgcolor
;        loadct = fltarr(n_elements(files))
;        linestyles = fltarr(n_elements(files))
        psym = fltarr(n_elements(files)) + 14;16
        psym = [fltarr(nhighres) + 14, fltarr(nlowres) + 6]
        obscolor = 100
        obssym = 2
    END
    7: BEGIN
        dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
        file799 = 'h799.cosmo25cmb.3072g14HBWK'
        key799 = 'h799'
        dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
        file516 = 'h516.cosmo25cmb.3072g14HBWK'
        key516 = 'h516'
        dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
        file986 = 'h986.cosmo50cmb.3072g14HBWK'
        key986 = 'h986'
        dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
        file603 = 'h603.cosmo50cmb.3072g14HBWK'
        key603 = 'h603'
        dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
        file258 = 'h258.cosmo50cmb.3072g14HMbwK'
        key258 = 'h258'
        dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
        file285 = 'h285.cosmo50cmb.3072g14HMbwK'
        key285 = 'h285'
        dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
        file239 = 'h239.cosmo50cmb.3072g14HMbwK'
        key239 = 'h239'
        dirs    = [ dir799, dir799, dir799, dir516, dir516, dir986, dir986, dir986, dir986, dir986, dir986, dir603, dir603, dir603, dir258, dir258, dir285, dir285, dir285, dir239]
        files   = [file799,file799,file799,file516,file516,file986,file986,file986,file986,file986,file986,file603,file603,file603,file258,file258,file285,file285,file285,file239]
        finalstep = '00512'
        haloid =  [ '1'   , '4'   , '6'   , '1'   , '2'   , '1'   , '2'   , '3'   , '8'   , '15'  , '16'  , '1'   , '2'    , '3'  , '1'   , '4'   , '1'   , '4'   , '9'   , '1'   ]
        haloid_split = [['1','4','6','0','0','0'],['1','2','0','0','0','0'],['1','2','3','8','15','16'],['1','2','3','0','0','0'],['1','4','0','0','0','0'],['1','4','9','0','0','0'],['1','0','0','0','0','0']]
 ;       gmass = [2.615e8,1.860e7,1.323e7,5.315e8,4.801e7,3.471e9,7.646e8,5.694e8,1.632e7,2.253e7,6.264e6,4.903e9,9.429e8,2.424e8,7.570e9,6.809e7,1.01e10,1.717e8,1.305e8,9.480e9]
        gmass = [3.82108e+06, 6.97347e+06, 1.86596e+07, 1.09378e+07, 1.64170e+07, 6.24969e+07, 1.47324e+08, 2.69942e+07, 1.61342e+08, 1.50572e+08, 8.98666e+07, 3.63764e+08,  5.93802e+08, 7.79724e+08, 6.51192e+08, 2.59430e+09, 2.45403e+09, 5.56453e+09,  5.50414e+09, 3.96146e+09] ;Particles where at least half atomic/molecular H
        gmass = [7.54535e+06, 1.02993e+07, 2.59515e+07, 1.48630e+07, 2.13353e+07, 9.35337e+07, 1.86572e+08, 6.18723e+07, 4.26302e+08, 2.68833e+08, 1.10110e+08, 1.02407e+09, 8.61151e+08, 1.09502e+09, 8.31080e+08, 5.45538e+09, 5.56057e+09, 7.35684e+09, 7.94138e+09, 5.78996e+09] ;Particle where density is greater than 0.1
        vfinal =  [55.1222,33.3827,26.5413,66.6838,42.7685,102.305,78.8046,76.3437,36.4566,28.7297,26.7203,115.102, 74.9589,50.1720,181.909,43.2878,163.970,52.1639,54.3400,165.470]

        outfiles = dirs + files + '.' + finalstep + '/' + files + '.' + finalstep
        outhalos = outfiles + '.halo.' + haloid 
        nhighres = 5
        nlowres = 15
        distunits = [fltarr(nhighres) + lu_c25, fltarr(nlowres) + lu_c50]
        massunits = [fltarr(nhighres) + mu_c25, fltarr(nlowres) + mu_c50]
        key  =    [ key799, key799, key799, key516, key516, key986, key986, key986, key986, key986, key986, key603, key603, key603, key258, key258, key285, key285, key285, key239] + ', ' + haloid
        mvir =    [ 2.4e10,  6.8e9,  4.4e9, 3.8e10, 1.5e10, 1.9e11, 5.9e10, 3.8e10, 1.1e10, 4.4e9,  3.1e9,  3.4e11, 1.0e11, 2.9e10, 7.7e11, 1.1e10, 8.8e11, 3.4e10, 1.2e10, 9.1e11] 
        masssort =[10,      2,      9,      1,      8,      15,     18,     4,      0,      13,     17,     3,      7,      6,      12,     5,      11,     14,     16,     19]

;        colors = fltarr(n_elements(files)) + fgcolor
;        loadct = fltarr(n_elements(files))
;        linestyles = fltarr(n_elements(files))
        psym = fltarr(n_elements(files)) + 14;16
        psym = [fltarr(nhighres) + 14, fltarr(nlowres) + 4]
        psym = [fltarr(nhighres) + 15,fltarr(nlowres) + 14]
        obscolor = 100
        obssym = 2
        label = textoidl(['3.2\times10^9','4.4\times10^9','4.4\times10^9','6.8\times10^9','1.1\times10^{10}','1.1\times10^{10}','1.2\times10^{10}','1.4\times10^{10}','2.4\times10^{10}','2.9\times10^{10}','3.4\times10^{10}','3.8\times10^{10}','3.8\times10^{10}','5.9\times10^{10}','10^{11}','1.8\times10^{11}','3.4\times10^{11}','7.7\times10^{11}','8.8\times10^{11}','9.1\times10^{11}']) + 'M' + sunsymbol()
    END
ENDCASE

n = n_elements(files)
IF keyword_set(color) THEN BEGIN
;    loadct,39
    rainbow_colors
    IF NOT keyword_set(ctables) THEN ctables = fltarr(n) + 39
    IF NOT keyword_set(obscolor) THEN obscolor = fgcolor
    IF NOT keyword_set(colors) THEN  colors  = (findgen(n) + 1)*254/n ELSE colors = colors
    IF NOT keyword_set(psym) THEN psym = fltarr(n) + 4
    IF NOT keyword_set(obssym) THEN obssym = 2
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 4
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 1.5
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN ctables = fltarr(n)
    IF NOT keyword_set(obscolor) THEN obscolor = 100
    If NOT keyword_set(colors) THEN  colors = (fltarr(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(broadband)) + 5
    IF NOT keyword_set(psym) THEN  psym = (findgen(n)+2)*2
    IF NOT keyword_set(obssym) THEN obssym = 2
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 4 ;thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2) 
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 1.5
ENDELSE

IF keyword_set(outplot) THEN BEGIN
    xsize = 18                  ;10*n
    ysize = 12                  ;12
ENDIF ELSE BEGIN
    xsize = 400                 ;*n
    ysize = 350                 ;475
ENDELSE

; ------------------------------------------ Halo Mags -------------------------------------------------------------------------------
IF pre_get_halo_mags THEN BEGIN
    uniqdir = uniq(dirs + files)
    FOR i = 0, n_elements(uniqdir) - 1 DO BEGIN
;    For i = 0, n - 1 DO BEGIN
        FOR izstep = 0, n_elements(zstep) - 1 DO BEGIN
            x = file_test(dirs[uniqdir[i]] + files[uniqdir[i]] + '.' + zstep[izstep])
            zstep_temp = zstep[izstep]
            IF NOT x THEN BEGIN
                IF zstep[izstep] EQ '00084' THEN zstep_temp = '00088'
            ENDIF
            print,dirs[uniqdir[i]] + files[uniqdir[i]] + '.' + zstep_temp + '/'
            cd,dirs[uniqdir[i]] + files[uniqdir[i]] + '.' + zstep_temp + '/'
            IF NOT file_test(dirs[uniqdir[i]] + files[uniqdir[i]] + '.' + zstep[izstep]+'.amiga_r200.halos.star99_K_ab.Mv.fits') THEN get_halo_mags,massunits[uniqdir[i]],distunits[uniqdir[i]],mag_sys = 'ab',/multiple;,/virial_radius
        ENDFOR
    ENDFOR
ENDIF

;reaccr_gas_character
IF pre_reaccr_gas_character THEN BEGIN
    FOR i = 0,n-1 DO BEGIN
        print,files[masssort[i]] + '.halo.' + haloid[masssort[i]]
        cd,dirs[masssort[i]]
        ;IF NOT file_test('grp' + haloid[i] + '.reaccrdisk_history.fits') THEN 
        reaccr_gas_character,accrhistory,dir = dirs[masssort[i]], finalid = haloid[masssort[i]]
     ENDFOR
ENDIF

IF pre_stellarz THEN BEGIN
    FOR i = 13,n-1 DO BEGIN
        print,files[masssort[i]] + '.halo.' + haloid[masssort[i]]
        cd,dirs[masssort[i]]
        stellarmetal_history,dir = dirs[masssort[i]],finalid = haloid[masssort[i]]
     ENDFOR
ENDIF

IF pre_metal THEN BEGIN
   FOR i = 17,18 DO BEGIN
      print,files[masssort[i]] + '.halo.' + haloid[masssort[i]]
      cd,dirs[masssort[i]]
      metal_history,finalid = haloid[masssort[i]],/inhalo
   ENDFOR
ENDIF

;event_gas_character
IF pre_event_gas_character THEN BEGIN
    FOR i = 0, n-1 DO BEGIN
        print,files[masssort[i]] + '.halo.' + haloid[masssort[i]]
        cd,dirs[masssort[i]]
        IF NOT file_test('grp' + haloid[masssort[i]] + '.reaccrdiskall_history.fits') THEN event_gas_character,history,'reaccrdiskall',dir = dirs[masssort[i]], finalid = haloid[masssort[i]]
        IF NOT file_test('grp' + haloid[masssort[i]] + '.relost_history.fits') THEN event_gas_character,history,'relost',dir = dirs[masssort[i]], finalid = haloid[masssort[i]]
        IF NOT file_test('grp' + haloid[masssort[i]] + '.reaccr_history.fits') THEN event_gas_character,history,'reaccr',dir = dirs[masssort[i]], finalid = haloid[masssort[i]]
        IF NOT file_test('grp' + haloid[masssort[i]] + '.reoutflow_history.fits') THEN event_gas_character,history,'reoutflow',dir = dirs[masssort[i]], finalid = haloid[masssort[i]]
     ENDFOR
ENDIF

IF pre_earlygas THEN BEGIN
    FOR i = 0, n - 1 DO BEGIN
        print,files[masssort[i]] + '.halo.' + haloid[masssort[i]]
        cd,dirs[masssort[i]]
        IF file_test(dirs[masssort[i]] + '/' + files[masssort[i]] + '.' + '00080' + '/' + files[masssort[i]] + '.' + '00080') $
        THEN firststep = '00080' $
        ELSE firststep = '00084'    
        earlygas,dirs[masssort[i]],files[masssort[i]],firststep,haloid[masssort[i]]
     ENDFOR
ENDIF

; ------------------------------------------ Cut out halo ----------------------------
IF pre_cuthalos THEN BEGIN
    outarray =['amiga.grp','FeMassFrac','OxMassFrac','HI','H2','coolontime','iord','HeI','HeII']
    FOR i = 0, n - 1 DO BEGIN
        print,outhalos[i]
        stop
        close,/all
        IF NOT file_test(outhalos[i]) OR keyword_set(redo) THEN tipsysatshi, outfiles[i], haloid[i], distunits[i], massunits[i], /cutout_rad, outarray = outarray
    ENDFOR
 ENDIF

IF pre_z_at_ejection THEN BEGIN
   FOR i = 0, n - 1 DO BEGIN
        print,files[masssort[i]] + '.halo.' + haloid[masssort[i]]
        cd,dirs[masssort[i]]
        z_at_ejection,dirs[masssort[i]],files[masssort[i]],haloid[masssort[i]]
     ENDFOR
ENDIF

; ------------------------------------------ Check contamination -------------------
IF check_contam THEN BEGIN
   FOR i = 0, n - 1 DO BEGIN
      rtipsy,outfiles[i] + '.halo.' + haloid[i] + '.std',h,g,d,s
      print,files[i] + '.' + finalstep + '.halo.' + haloid[i],n_elements(where(d.mass GT MIN(d.mass))),total(d[where(d.mass GT MIN(d.mass))].mass)/total(d.mass)
   ENDFOR
ENDIF

; ------------------------------------------ Create starlog file for halo ----------------------------
IF pre_make_halo_starlog THEN BEGIN
    FOR i = 3, n - 1 DO BEGIN
        print,outhalos[i]
        cd,dirs[i]
        make_halo_starlog,files[i],halo = haloid[i],/molecularH
    ENDFOR
ENDIF

; ------------------------------------------- multiHaloGuo -----------------------------------------------------------------------------
;/astro/users/christensen/code/IDL/HIcubes/multiHaloGuo
IF plot_multiHaloGuo THEN BEGIN
    uniqdir = uniq(outfiles)
    multihaloguo,outfiles[uniqdir],key = key[uniqdir],psym = psym[uniqdir],color = colors[uniqdir], obscolor = obscolor, obssym = obssym, ctables = ctables[uniqdir],outplot = outplot,/things,/bw,/kcorrect,xrange = [9,12],yrange = [6,11],haloids = haloid_split,moster = 'moster.stars.z0'
ENDIF

; ------------------------------------------- baryonicfrac -----------------------------------------------------------------------------
;/astro/users/christensen/code/IDL/procedures/baryonicfrac.pro
IF plot_baryonicfrac THEN BEGIN
    IF NOT keyword_set(gmass) THEN gmass = tully_fisher_obs_gasmass(outfiles + '.halo.' + haloid, massunits,smass = smass)
    uniqdir = uniq(outfiles)
    baryonicfrac,outfiles[uniqdir],haloid_split,gmass,/kcorrect,xrange = [9,12],color = colors[uniqdir], ctables = ctables[uniqdir],psym = psym[uniqdir],symsizes = symsizes[uniqdir],thicks = thicks[uniqdir], outplot = outplot,yrange=[3e-3,1]
ENDIF

;--------------------------- Vcirc -------------
IF plot_vcirc THEN BEGIN
    yrange_vcirc = [0,300]
    maxdistance = 50
    verbose = 1
    nbins = 500
    linestyles_vcirc = fltarr(n_elements(outfiles))
    IF keyword_set(outplot) THEN BEGIN
        IF keyword_set(color) THEN vcirc,outfiles + '.halo.' + haloid,massunits,distunits,keys = key,thicks = thicks,linestyle = linestyles_vcirc,yrange = yrange_vcirc,maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,verbose = verbose,outfile = outplot + '_vcirc.eps',color = colors,nbins = nbins $
                              ELSE vcirc,outfiles + '.halo.' + haloid,massunits,distunits,keys = key,thicks = thicks,linestyle = linestyles_vcirc,yrange = yrange_vcirc,maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,verbose = verbose,outfile = outplot + '_vcirc.eps',nbins = nbins
    ENDIF ELSE BEGIN
        IF keyword_set(color) THEN vcirc,outfiles + '.halo.' + haloid,massunits,distunits,keys = key,thicks = thicks,linestyle = linestyles_vcirc,yrange = yrange_vcirc, maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,verbose = verbose,color = colors,nbins = nbins   $
                              ELSE vcirc,outfiles + '.halo.' + haloid,massunits,distunits,keys = key,thicks = thicks,linestyle = linestyles_vcirc,yrange = yrange_vcirc, maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,verbose = verbose,nbins = nbins
    ENDELSE
    print,vfinal
ENDIF

;---------------------------------HI macro-------------------------------
IF do_writeHImacro THEN BEGIN
    FOR i = 0, n - 1 DO BEGIN
        IF vfinal[i] LT 83 THEN dwarfv = 1 ELSE dwarfv = 0
        cd,dirs[i] + files[i] + '.' + finalstep
        writehimacro,files[i] + '.' + finalstep + '.halo.' + haloid[i],dist_unit = distunits[i],mass_unit = massunits[i],radius = 24.0,/molecularH,dwarfv = dwarfv
    ENDFOR
ENDIF

;----------------------------------
IF do_line_profile THEN BEGIN
    FOR i = 1, n - 1 DO BEGIN
        angle = 45
        angle_str = strtrim(angle,2) + '.'
        cfile = files[i] + '.' + finalstep + '.halo.' + haloid[i]
        print,cfile
        IF file_test(dirs[i] + files[i] + '.' + finalstep + '/' + cfile + '.cube.' + angle_str + 'fits') THEN BEGIN
            cube = read_cube_fits(dirs[i] + files[i] + '.' + finalstep + '/' + cfile + '.cube.' + angle_str + 'fits',header,expansion = expansion)
            shead = convertheader(header,max(cube),min(cube)) ;Converts header to strings in correct format to output to a fits file through MWRFITS
            total_line_profile,cube,header,vaxis,spectrum,w,outfile = cfile,doplot = doplot,xyrange = xyrange ;Calculate line width as if from single dish data (i.e., don't smooth and remove sensetivity as if from THINGS)
            width = get_width_fit(vaxis, spectrum, /doplot, gal = cfile,incl=angle, pa=0)                     ;,sparange = sparange)
        
            print,'Width [km/s]:',width
            print,'V_max [km/s]:',width/2/sin(angle*!pi/180)
            openw,1,'linewidths.txt'
            printf,1,'Width [km/s]             Width/2 [km/s]             V_max [km/s]'
            printf,1,width,width/2,width/2/sin(angle*!pi/180)
            close,1
        ENDIF ELSE print,'No file'
    ENDFOR
ENDIF

; ------------------------------------------ Baryonic Tully Fisher -----------------------------------------------------
IF plot_baryonic_tully_fisher_obs THEN BEGIN
    IF NOT keyword_set(vfinal) THEN BEGIN
        vfinal = fltarr(n)
        FOR i = 0, n - 1 DO BEGIN
            data = read_stat_struc_amiga(outfiles[i] + '.amiga.stat')
            ind = where(data.group eq haloid[i])
            vfinal[i] = data[ind].vc
        ENDFOR
    ENDIF
    IF NOT keyword_set(gmass) THEN gmass = tully_fisher_obs_gasmass(outfiles + '.halo.' + haloid, massunits)
    tully_fisher_obs_btf,outfiles,vfinal,gmass,halo = haloid,/kcorrect,symbols = psym,symsizes = symsizes,thicks = thicks,ctables = ctables,obscolor = obscolor,color = 60,outfile = outplot;,key = key
ENDIF

; ------------------------------------------ Tully Fisher -----------------------------------------------------
IF plot_tully_fisher_obs THEN BEGIN
    IF NOT keyword_set(vfinal) THEN BEGIN
        vfinal = fltarr(n)
        FOR i = 0, n - 1 DO BEGIN
            data = read_stat_struc_amiga(outfiles[i] + '.amiga.stat')
            ind = where(data.group eq haloid[i])
            vfinal[i] = data[ind].vc
        ENDFOR
    ENDIF
    IF NOT keyword_set(gmass) THEN gmass = tully_fisher_obs_gasmass(outfiles + '.halo.' + haloid, massunits)
    tully_fisher_obs,    outfiles,vfinal,gmass,halo = haloid,/kcorrect,symbols = psym,symsizes = symsizes,thicks = thicks,ctables = ctables,obscolor = obscolor,color = color,outfile = outplot
ENDIF

;------------------------------------------ MZR ------------------------------

IF pre_mzr THEN BEGIN
    uniqdir = uniq(dirs)
    stop
    FOR i = 0, n_elements(uniqdir) - 1 DO BEGIN
        FOR izstep = n_elements(zstep) - 1, n_elements(zstep) - 1 DO BEGIN ;3 DO BEGIN ;n_elements(zstep) - 1 DO BEGIN
            x = file_test(dirs[uniqdir[i]] + files[uniqdir[i]] + '.' + zstep[izstep])
            zstep_temp = zstep[izstep]
            IF NOT x THEN BEGIN
                IF zstep[izstep] EQ '00084' THEN zstep_temp = '00088'
            ENDIF
            filename = dirs[uniqdir[i]] + files[uniqdir[i]] + '.' + zstep_temp + '/' + files[uniqdir[i]] + '.' + zstep_temp
            print,filename
            IF file_test(filename + '.amiga.stat') THEN BEGIN
;                IF NOT file_test(filename + '.metals.fits') THEN BEGIN
;            IF file_test(filename + '.amiga_vir.halos.star99_K_ab.Mv.fits') THEN metals = mzr(filename,/obs) ELSE 
                    metals = mzr(filename,/obs) 
                    mwrfits,metals,filename + '.metals.fits',/create
;                ENDIF
            ENDIF
        ENDFOR
    ENDFOR
ENDIF

IF plot_galMS THEN BEGIN
;    uniqdir = uniq(outfiles)
    halos = intarr(1,n)
    halosz = intarr(1,n,n_elements(zstep))
    redshiftz = findgen(n_elements(zstep))
    nstop = n_elements(dirs) - 1 ;10
    FOR izstep = 0, n_elements(zstep) - 1 DO BEGIN
        filename = dirs + files + '.'
        FOR i = 0, nstop DO BEGIN ;n - 1 DO BEGIN
            x = file_test(dirs[i] + files[i] + '.' + zstep[izstep])
            zstep_temp = zstep[izstep]
            IF NOT x THEN BEGIN
                IF zstep[izstep] EQ '00084' THEN zstep_temp = '00088'
            ENDIF 
            filename[i] = filename[i] + zstep_temp + '/' + files[i] + '.' + zstep_temp
            spawn,'ls ' + dirs[i] + '*.grp' + haloid[i] + '.haloid.dat',files_haloid
            readcol,files_haloid[0],files_steps,halos_steps,format='a,l',/silent
            indstep = where(files_steps EQ files[i] + '.' + zstep_temp + '/' + files[i] + '.' + zstep_temp)
  ;          IF (indstep[0] EQ -1) OR (n_elements(indstep) GT 1) THEN stop
            halos[0,i] = halos_steps[indstep]
            halosz[0,i,izstep] = halos_steps[indstep]
        ENDFOR
        IF NOT keyword_set(filenamez) THEN filenamez = filename[0:nstop] ELSE filenamez = [filenamez, filename[0:nstop]]
;        stop
;        mzr_plot,filename[0:10],halos = halos[0,0:10],psym = psym[0:10],symsize = symsizes[0:10],thicks = thicks[0:10],ctables = ctables[0:10],outfile = outplot,obsct = obsct,/readfile,redshift = redshift,/color;,/simStellarMass;,key = key
        IF zstep[izstep] EQ '00084' THEN redshift = 3
        IF zstep[izstep] EQ '00116' THEN redshift = 2.2
        IF zstep[izstep] EQ '00112' THEN redshift = 2.3
        IF zstep[izstep] EQ '00260' THEN redshift = 0.8
        IF zstep[izstep] EQ '00512' THEN redshift = 0
        redshiftz[izstep] = redshift
    ENDFOR
    galMS,filenamez,halos = halosz[0,0:nstop,*],psym = psym[0:nstop],symsize = symsizes[0:nstop],thicks = thicks[0:nstop],ctables = ctables[0:nstop],outfile = outplot,/readfile,/color,xrange = [6.0,12],yrange = [-4,4],multiz = n_elements(zstep),redshift= redshiftz
    galMS,filenamez,halos = halosz[0,0:nstop,*],psym = psym[0:nstop],symsize = symsizes[0:nstop],thicks = thicks[0:nstop],ctables = ctables[0:nstop],outfile = outplot,/readfile,/color,xrange = [6.0,12],yrange = [-11,-7],multiz = n_elements(zstep),redshift= redshiftz,/specific
ENDIF


IF plot_mzr THEN BEGIN
;    uniqdir = uniq(outfiles)
    halos = intarr(1,n)
    halosz = intarr(1,n,n_elements(zstep))
    nstop = n_elements(dirs) - 1
    FOR izstep = 0, n_elements(zstep) - 1 DO BEGIN
        filename = dirs + files + '.'
        FOR i = 0, nstop DO BEGIN ;n - 1 DO BEGIN
            x = file_test(dirs[i] + files[i] + '.' + zstep[izstep])
            zstep_temp = zstep[izstep]
            IF NOT x THEN BEGIN
                IF zstep[izstep] EQ '00084' THEN zstep_temp = '00088'
            ENDIF 
            filename[i] = filename[i] + zstep_temp + '/' + files[i] + '.' + zstep_temp
            spawn,'ls ' + dirs[i] + '*.grp' + haloid[i] + '.haloid.dat',files_haloid
            readcol,files_haloid[0],files_steps,halos_steps,format='a,l',/silent
            indstep = where(files_steps EQ files[i] + '.' + zstep_temp + '/' + files[i] + '.' + zstep_temp)
  ;          IF (indstep[0] EQ -1) OR (n_elements(indstep) GT 1) THEN stop
            halos[0,i] = halos_steps[indstep]
            halosz[0,i,izstep] = halos_steps[indstep]
        ENDFOR
        redshift = 0
        IF zstep[izstep] EQ '00084' THEN redshift = 3
        IF zstep[izstep] EQ '00116' THEN redshift = 2.2
        IF zstep[izstep] EQ '00112' THEN redshift = 2.3
        IF zstep[izstep] EQ '00260' THEN redshift = 0.8
        IF NOT keyword_set(filenamez) THEN filenamez = filename[0:nstop] ELSE filenamez = [filenamez, filename[0:nstop]]
;        stop
        mzr_plot,filename[0:nstop],halos = halos[0,0:nstop],psym = psym[0:nstop],symsize = symsizes[0:nstop],thicks = thicks[0:nstop],ctables = ctables[0:nstop],outfile = outplot,obsct = obsct,/readfile,redshift = redshift,/color;,/simStellarMass;,key = key
    ENDFOR
;    mzr_plot,filenamez,halos = halosz[0,0:nstop,*],psym = psym[0:nstop],symsize = symsizes[0:nstop],thicks = thicks[0:nstop],ctables = ctables[0:nstop],outfile = outplot,obsct = obsct,/readfile,redshift = redshift,/color,multiz = n_elements(zstep)
ENDIF

IF plot_mzr_stellar THEN BEGIN
;    uniqdir = uniq(outfiles)
    halos = intarr(1,n)
    izstep = n_elements(zstep) - 1
    filename = dirs + files + '.'
    FOR i = 0, n - 1 DO BEGIN
        x = file_test(dirs[i] + files[i] + '.' + zstep[izstep])
        zstep_temp = zstep[izstep]
        IF NOT x THEN BEGIN
            IF zstep[izstep] EQ '00084' THEN zstep_temp = '00088'
        ENDIF 
        filename[i] = filename[i] + zstep_temp + '/' + files[i] + '.' + zstep_temp
        spawn,'ls ' + dirs[i] + '*.grp' + haloid[i] + '.haloid.dat',files_haloid
        readcol,files_haloid[0],files_steps,halos_steps,format='a,l',/silent
        indstep = where(files_steps EQ files[i] + '.' + zstep_temp + '/' + files[i] + '.' + zstep_temp)
        IF (indstep[0] EQ -1) OR (n_elements(indstep) GT 1) THEN stop
        halos[0,i] = halos_steps[indstep]
    ENDFOR
    print,filename
    mzr_plot,filename,halos = halos,psym = psym,symsize = symsizes,thicks = thicks,ctables = ctables,outfile = outplot,obsct = obsct,/readfile,/stellar,xrange = [3,12],/color;,/simStellarMass ;,key = key
 ENDIF

;------------------------------------------ Tracking the Mass ---------------
IF plot_track_mass THEN BEGIN
    track_mass_plot,dir,filenames,key = keys,linestyle = linestyles,thick = thicks
ENDIF

;------------------------------------------ Cumlative Feedback ------------
IF plot_fbcum THEN BEGIN
;    outflow_plot,dsirtop,outplot = outplot, keys = keys, color = colors, thicks = thicks, linestyles = linestyles,label = label,ctables = ctables,yrange_fbcum = yrange_fbcum
;    eject_plots,dirs,halo = haloid,outplot = outplot, keys = keys, color = colors, thicks = thicks,label = label,ctables = ctables,yrange_fbcum = [0,0.1],formatthick = formatthick,linestyles = linestyles ,/unscale

    eject_plots,dirs[masssort],halo = haloid[masssort],outplot = outplot, keys = key[masssort], color = colors, thicks = thicks,label = label,ctables = ctables,yrange_fbcum = [0,0.1],formatthick = formatthick,molecularH = molecularH,linestyles = linestyles,/massloading;,/unscale;/unscale;,/massloading
;    stop
;    eject_plots,dirs[masssort],halo = haloid[masssort],outplot = outplot, keys = key[masssort], color = colors[masssort], thicks = thicks[masssort],label = label[masssort],ctables = ctables[masssort],yrange_fbcum = [0,0.1],formatthick = formatthick,molecularH = molecularH,linestyles = linestyles
;    stop
;    eject_plots,dirs[masssort],halo = haloid[masssort],outplot = outplot, keys = key[masssort], color = colors, thicks = thicks,label = label,ctables = ctables,yrange_fbcum = [0,0.1],formatthick = formatthick,molecularH = molecularH,linestyles = linestyles,/unscale
ENDIF

IF find_reaccr THEN find_reaccr,dirs[masssort],files[masssort],halo = haloid[masssort]

;--------------------- Eject v Mass ------------------------------
IF plot_ejectz_v_mass THEN ejectz_v_mass,dirs[masssort[0:n_elements(masssort)-1]],files[masssort[0:n_elements(masssort)-1]],halo = haloid[masssort[0:n_elements(masssort)-1]],/colors,outplot = outplot,z_cut = [2,1,0.5,-1e-10],symbols = [psym[0],14,17,18,16,34],formatthick = formatthick,gmass = gmass;,/rewrite
IF plot_metals_v_mass THEN metals_v_mass,dirs,files,haloid,outplot = outplot,formatthick = formatthick,/color,/stellarmass
IF plot_metals_v_mass THEN $
  IF keyword_set(outplot) THEN metals_v_mass,dirs,files,haloid,outplot = outplot+'z2',formatthick = formatthick,/color,/stellarmass,finalstep = '00128' ELSE metals_v_mass,dirs,files,haloid,formatthick = formatthick,/color,/stellarmass,finalstep = '00128'

;noh603 = [0,1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19]
;IF plot_times_ejected THEN times_ejected,dirs,files,halo = haloid,colors = colors,outplot = outplot,/outflow
IF plot_times_ejected THEN times_ejected,dirs,files,halo = haloid,colors = colors,outplot = outplot,/outflow,/accrdisk
;IF plot_times_ejected THEN times_ejected,dirs,files,halo = haloid,colors = colors,outplot = outplot
IF plot_times_ejected THEN times_ejected,dirs,files,halo = haloid,colors = colors,outplot = outplot,/accrdisk
;IF plot_times_expelled THEN times_ejected,dirs,files,halo = haloid,colors = colors,outplot = outplot,/expelled
IF plot_times_expelled THEN times_ejected,dirs,files,halo = haloid,colors = colors,outplot = outplot,/accrdisk,/expelled
;IF run_times_cycling THEN time_cycling_exp,dirs,files,halo = haloid
IF run_times_cycling THEN time_cycling,dirs[masssort[0:n_elements(masssort) - 1]],files[masssort[0:n_elements(masssort) - 1]],halo = haloid[masssort[0:n_elements(masssort) - 1]],/nowrite
IF plot_times_cycling THEN plot_time_cycling,dirs,files,halo = haloid,colors = [50,254],outplot = outplot,symbols = [14,15];,/expelled 
;IF plot_times_cycling THEN plot_time_cycling,dirs,files,halo = haloid,colors = colors,outplot = outplot
;IF plot_reeject_r THEN plot_half_eject,dirs[masssort],files[masssort],halo = haloid[masssort],colors = colors,outplot = outplot,/normalize,symbols = [psym[0],psym[0] + 1],formatthick = formatthick,sfr = 0;,/stellarmass
IF plot_reeject_r THEN plot_half_eject,dirs[masssort],files[masssort],halo = haloid[masssort],colors = colors,outplot = outplot,/normalize,symbols = [psym[0],psym[0] + 1],formatthick = formatthick,sfr = 1;,/stellarmass
IF plot_reaccr_r THEN plot_half_accr,dirs[masssort],files[masssort],halo = haloid[masssort],colors = colors,outplot = outplot,formatthick = formatthick
IF plot_angmom THEN plot_angmom,dirs[masssort],files[masssort],halo = haloid[masssort],colors = colors,outplot = outplot,formatthick = formatthick
IF plot_reeject_z THEN BEGIN
    reeject_metallicity,reverse(dirs[masssort]),reverse(files[masssort]),finalid =reverse(haloid[masssort]),colors = reverse(colors),outplot = outplot,formatthick = formatthick,/normalize,avez = avez,red_avez = red_avez ;keys = reverse(key[masssort])
    reeject_metallicity,reverse(dirs[masssort]),reverse(files[masssort]),finalid =reverse(haloid[masssort]),colors = reverse(colors),outplot = outplot,formatthick = formatthick,/normalize,avez = avezabs,red_avez = red_avezabs,/absolute ;keys = reverse(key[masssort])

    formatplot,outplot = outplot,thick = formatthick
    IF keyword_set(outplot) THEN BEGIN
        xsize = 18
        ysize = 12
    ENDIF ELSE BEGIN
        xsize = 800
        ysize = 500
     ENDELSE
    loadct,0
    readcol,'~/etaflux.txt',mvir,vvir,eta_inner,etaz_inner,eta_outer,etaz_outer
    IF keyword_set(outplot) THEN device,filename = outplot + '_vm_zm.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*1.8,xoffset =  2,yoffset =  2 ELSE window, 0, xsize = xsize, ysize = ysize*1.8
    multiplot,[1,2]
    plot_z_m,avezabs,reverse(dirs[masssort]),reverse(files[masssort]),halo = reverse(haloid[masssort]),outplot = outplot,/normalize,formatthick = formatthick,colors = [colors[1]],symbols = [psym[0] + 1],red_avez = red_avezabs,/absolute,/legend;,/stellarmass
    stop
    multiplot
    plot_z_m,avez,   reverse(dirs[masssort]),reverse(files[masssort]),halo = reverse(haloid[masssort]),outplot = outplot,/normalize,formatthick = formatthick,colors = [colors[1]],symbols = [psym[0] + 1],red_ave = red_avez;,/stellarmass
    multiplot,/reset
    IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF
IF plot_reeject_v THEN BEGIN
;   FOR i = 0, n_elements(masssort) - 1 DO BEGIN
;      cd,dirs[masssort[i]]
;      track_vc,finalid = haloid[masssort[i]]
;   ENDFOR
    reeject_velocity,dirs[masssort],files[masssort],finalid =haloid[masssort],colors = colors,outplot = outplot,/normalize,formatthick = formatthick,medv_ej = medv_ej,     stdv_ej = stdv_ej,    medv_exp = medv_exp,      stdv_exp = stdv_exp,/expell
    reeject_velocity,dirs[masssort],files[masssort],finalid =haloid[masssort],colors = colors,outplot = outplot,/normalize,formatthick = formatthick,medv_ej = medvscale_ej,stdv_ej = stdvscale_ej,medv_exp = medvscale_exp,stdv_exp = stdvscale_exp,/vs_step
;   reeject_velocity,dirs[masssort],files[masssort],finalid =haloid[masssort],colors = colors,outplot = outplot,/normalize,formatthick = formatthick,medv_ej = medvscale_ej,stdv_ej = stdvscale_ej,medv_exp = medvscale_exp,stdv_exp = stdvscale_exp,vfinal = vfinal[masssort],/vscale
;   reeject_velocity,dirs[masssort],files[masssort],finalid =haloid[masssort],colors = colors,outplot = outplot,/normalize,formatthick = formatthick,medv_ej = medvscale_ej,stdv_ej = stdvscale_ej,medv_exp = medvscale_exp,stdv_exp = stdvscale_exp,/vscale
;   reeject_velocity,dirs[masssort],files[masssort],finalid =haloid[masssort],colors = colors,outplot = outplot,/normalize,formatthick = formatthick,medv_ej = medvscale_ej,stdv_ej = stdvscale_ej,medv_exp = medvscale_exp,stdv_exp = stdvscale_exp,/vs_step ,/expell
;   reeject_velocity,dirs[masssort],files[masssort],finalid =haloid[masssort],colors = colors,outplot = outplot,/normalize,formatthick = formatthick,medv_ej = medvscale_ej,stdv_ej = stdvscale_ej,medv_exp = medvscale_exp,stdv_exp = stdvscale_exp,/vscale,/expell
; reeject_velocity,reverse(dirs[masssort]),reverse(files[masssort]),finalid =reverse(haloid[masssort]),colors = reverse(colors),outplot = outplot,/normalize
    IF NOT keyword_set(outplot) THEN stop

    formatplot,outplot = outplot,thick = formatthick
    IF keyword_set(outplot) THEN BEGIN
        xsize = 18
        ysize = 12
    ENDIF ELSE BEGIN
        xsize = 800
        ysize = 500
     ENDELSE
    IF keyword_set(outplot) THEN  device,filename = outplot + '_vm_velm_multi.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*1.8,xoffset =  2,yoffset =  2 ELSE window, 0, xsize = xsize, ysize = ysize*1.8
    multiplot,[1,2]
    plot_v_m,medv_ej,     medv_exp,     stdv_ej,     stdv_exp,     dirs[masssort],files[masssort],halo = haloid[masssort],outplot = outplot,/normalize,formatthick = formatthick,colors = [colors[1]],/absolute;,vfinal = vfinal[masssort];,symbols = [psym[1]]
    multiplot
    plot_v_m,medvscale_ej,medvscale_exp,stdvscale_ej,stdvscale_exp,dirs[masssort],files[masssort],halo = haloid[masssort],outplot = outplot,/normalize,formatthick = formatthick,colors = [colors[1]];,vfinal = vfinal[masssort];,symbols = [psym[1]]
;    plot_v_m,medv_ej,medv_exp,stdv_ej,stdv_exp,dirs[masssort],files[masssort],halo = haloid[masssort],outplot = outplot,/normalize,formatthick = formatthick,colors = [colors[1]]
    multiplot,/reset
    IF keyword_set(outplot) THEN device, /close ELSE stop

    IF keyword_set(outplot) THEN  device,filename = outplot + '_vm_velm_abs.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 0, xsize = xsize, ysize = ysize
    plot_v_m,medv_ej,     medv_exp,     stdv_ej,     stdv_exp,     dirs[masssort],files[masssort],halo = haloid[masssort],outplot = outplot,/normalize,formatthick = formatthick,colors = [colors[1]],/absolute ;,vfinal = vfinal[masssort];,symbols = [psym[1]
    IF keyword_set(outplot) THEN device, /close ELSE stop

    IF keyword_set(outplot) THEN  device,filename = outplot + '_vm_velm.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 0, xsize = xsize, ysize = ysize
    plot_v_m,medvscale_ej,medvscale_exp,stdvscale_ej,stdvscale_exp,dirs[masssort],files[masssort],halo = haloid[masssort],outplot = outplot,/normalize,formatthick = formatthick,colors = [colors[1]];,vfinal = vfinal[masssort];,symbols = [psym[1]]
;    plot_v_m,medv_ej,medv_exp,stdv_ej,stdv_exp,dirs[masssort],files[masssort],halo = haloid[masssort],outplot = outplot,/normalize,formatthick = formatthick,colors = [colors[1]]
    IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF

IF check_coolontime THEN BEGIN
   coolontime_check,dirs[masssort],files[masssort],haloid[masssort]
ENDIF

IF plot_sfh_reaccr THEN BEGIN
   FOR i = 0, n - 1 DO BEGIN
      IF keyword_set(outplot) THEN sfh_reaccr,dirs[i],halo = haloid[i],formatthick = formatthick,color = color,outplot = files[i] ELSE sfh_reaccr,dirs[i],halo = haloid[i],formatthick = formatthick,color = color
   ENDFOR
ENDIF

IF plot_inflow_outflow_historyz THEN BEGIN
; i =4, 12, 13 and 15 have problems
    outind = [0,11,14,16]
    outind = [8];15
    outind = [1,6,10,17]
    inflow_outflow_historyz,dirs[masssort[outind]],haloid = haloid[masssort[outind]],outplot = outplot,colors = colors[outind],label = label[outind],/plot_halo,/debug
 ;  inflow_outflow_historyz,dirs[masssort],haloid = haloid[masssort],outplot = outplot,colors = colors,label = label,/plot_halo,/debug
   ;inflow_outflow_historyz,dirs[masssort[0:n_elements(masssort) - 1]],haloid = haloid[masssort[0:n_elements(masssort) - 1]],outplot = outplot,colors = colors,label = label,pmulti = [4,5];/plot_halo;,/debug;,pmulti = [4,5];,/debug
ENDIF

IF plot_inflow_outflow_metallicity THEN BEGIN
    outind = [0,11,14,16]
    outind = [0,1,2,3,5,6,7,8,9,10,12,13,15,17,18,19]
    outind = [1,6,10,17]
    inflow_outflow_metallicity_plot,dirs[masssort[outind]],haloid = haloid[masssort[outind]],outplot = outplot,colors = colors[outind],label = label[outind],pmulti = [0,2,2,0,0]
;    inflow_outflow_metallicity_plot,dirs[masssort],haloid = haloid[masssort],outplot = outplot,colors = colors,label = label
ENDIF

IF plot_outflowr THEN BEGIN
   plot_outflowr,dirs[masssort],files[masssort],finalid = haloid[masssort],outplot = outplot
ENDIF


IF NOT keyword_set(outplot) THEN stop
END
