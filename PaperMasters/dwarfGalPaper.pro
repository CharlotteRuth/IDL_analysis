;.r coolrate
;.r checkAbundMet.pro
;.r opticalRadii
;.r Halpha
;.r schmidtlaw_res_obs.pro
;.r schmidtlaw_global_obs.pro
;.r tully_fisher_obs
;.r sfh
;.r color_prof.pro
;.r densRadius.pro
;.r toomreRadius.pro
;.r dwarfGalPaper
;.r mzr.pro
;.r mass_evol.pro
;outplot = '~/plots/h516.cosmo25cmb.3072g.paper'
;outplot = '~/plots/h516.cosmo25cmb.3072g.papercolor'
;outplot = '~/plots/h516.cosmo25cmb.3072g.paperres'
;outplot = '~/plots/h516.cosmo25cmb.3072g.paper472'

pro dwarfGalPaper,outplot = outplot,color = color,verbose = verbose,formatthick = formatthick
idldir = '/astro/users/christensen/code/IDL/'
;resolve_routine,'multiHaloGuo'
;resolve_routine,'tully_fisher_obs'
;resolve_routine,'matchHalos',/is_function
;resolve_routine,'photometricProf'
;resolve_routine,'opticalRadii',/is_function
;resolve_routine,'Halpha',/is_function
;resolve_routine,'schmidtlaw_global_obs',/is_function
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23

plot_get_halo_mags = 0 ;pre
plot_coolrate = 0
plot_checkAbundMet = 0 ;dangerous
plot_checkAbundRes = 0 ;dangerous
plot_schmidtlaw_res_obs = 0 ;pre
plot_schmidtlaw_res_obs_master_out = 0
plot_schmidtlaw_global_obs = 0
plot_magcolor = 0
plot_tully_fisher_obs = 0
plot_sfh = 1
plot_color_prof = 0
plot_toomreRadius = 0
plot_fft_result = 0
plot_densRadius = 0;dangerous
plot_sfgas = 0;dangerous
plot_make_jpeg = 0
plot_mzr = 0

IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    obscolor = fgcolor
ENDIF ELSE BEGIN    loadct,0    
    obscolor = 150
ENDELSE
obssym = 2
obssymsize = 2;1.5

x=1
CASE x OF 
    1: BEGIN
       spawn,'hostname',hostname
       IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
;        steps = ['00492','00512']
       steps = ['00492','00492']
;       steps = ['00512','00512']
       dir = prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir',$
                        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir']
        files = prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[0],$
                          'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[1]]
        n = N_ELEMENTS(files)
        pfiles = prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
        sflog = prefix +  ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.starlog',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.starlog']
        filebases = ['h516.cosmo25cmb.3072g1MBWK.'  + steps[0] + '.halo.1',$
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.halo.1']
        filenames = ['h516.cosmo25cmb.3072g1MBWK.'  + steps[0], $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[1]]
        outext = ['met',$
                  'H2']
        useH2 = [0,1]
        halos = [1,1]
        halos_str = strtrim(halos,2)
        keys = ['DnoH2','DH2']
        keys = [textoidl('Metals, no H_2'),textoidl('Metals, H_2')]
        distunits = [25000.,25000.]
        massunits  = [2.310e15,2.310e15]
        angle = [45.0,45.0]
        extno = [15,15] ;45 degrees
        intq = [32,33]
        rotateAngle = [5,5]
        filternums = [13,13]    ;[13,14,13]
        cameras = extno       ;[14,15,14]

;        psym = [17,6]

        IF KEYWORD_SET(color) THEN BEGIN
            if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
            symbols = fltarr(n) + 4
            symsizes = fltarr(n) + 2 ;1.5
            thicks = fltarr(n) + 4 ;2
            linestyles = fltarr(n)

            colors = [80,245]
;            colors = [120,245]
            symbols = [16,16]
        ENDIF ELSE BEGIN
            colors = fltarr(n) + fgcolor
            symbols = (findgen(n)+2)*2
            symsizes = fltarr(n) + 2;1.5
            thicks = fltarr(n) + 2
            linestyles = [2,0]
        ENDELSE
        IF keyword_set(outplot) THEN thicks_psym = fltarr(n) + 6 ELSE thicks_psym = fltarr(n) + 2

        yrange_vcirc = [0,75]
        yrange_photo = [26,19]
        maxdistance_photo = 4
        maxdistance = 6
        yrange_SFH = [0,0.4]

        velocities_true = [53.2281,56.8234]
        velocities = velocities_true ;taken from flat part of rotation curve
        velocitiesHI = [53.406370,54.337259] ;from HI line widths        
        gh516_H = {B:-15.6996847198888, V:-16.0784155581285, alpha:0.81438442,alpha_gal:1.02964, alphaV:0.00000001,alpha_galV:0.000001} ;h516.cosmo25cmb.3072g14HBWK
        gh516_M = {B:-15.2652157189731, V:-15.7550592733375, alpha:0.79430061,alpha_gal:0.770795,alphaV:0.79689074,alpha_galV:0.777515}
        simdata = [gh516_M,gh516_H] 

    END
    2: BEGIN
        prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
        steps = ['00512','00492','00512','00512']
        dir = prefix + ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g2MBWK/steps/h516.cosmo25cmb.1536g2MBWK.'   + steps[0] + '.dir',$
                        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[1] + '.dir',$
                        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.' + steps[2] + '.dir',$
                        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[3] + '.dir']
        files = prefix + ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g2MBWK/steps/h516.cosmo25cmb.1536g2MBWK.'   + steps[0] + '.dir/h516.cosmo25cmb.1536g2MBWK.'  + steps[0],$
                          'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[1] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[1],$
                          'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.' + steps[2] + '.dir/h516.cosmo25cmb.1536g14HBWK.' + steps[2],$
                          'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[3] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[3]]
        pfiles = prefix + ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g2MBWK/h516.cosmo25cmb.1536g2MBWK.param',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.param',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
        sflog = prefix +  ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g2MBWK/h516.cosmo25cmb.1536g2MBWK.starlog',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.starlog',$
                           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.starlog',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.starlog']
        filebases = ['h516.cosmo25cmb.1536g2MBWK.'  + steps[0] + '.halo.1',$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[1] + '.halo.1',$
                     'h516.cosmo25cmb.1536g14HBWK.' + steps[2] + '.halo.1',$
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[3] + '.halo.1']
        filenames = ['h516.cosmo25cmb.1536g2MBWK.'  + steps[0],$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[1], $
                     'h516.cosmo25cmb.1536g14HBWK.' + steps[2],$
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[3]]
        outext = ['met, 1536',$
                  'met, 3072',$
                  'H2, 1536' ,$
                  'H2, 3072']
        n = N_ELEMENTS(files)
        useH2 = [0,0,1,1]
        halos = [1,1,1,1]
        halos_str = strtrim(halos,2)
        keys = ['DnoH2, lr','DnoH2','DH2 lr','DH2']
;        psym = [17,6]
        distunits = [25000.,25000.,25000.,25000.]
        massunits  = [2.310e15,2.310e15,2.310e15,2.310e15]
        linestyles = [2,2,0,0]
        thicks = [1,2,1,2]
        filternums = [13,13,13,13] ;[13,14,13]
        cameras = [14,14,14,14] ;[14,15,14]
        yrange_vcirc = [0,75]
        yrange_photo = [26,19]
        maxdistance_photo = 4
        yrange_SFH = [0,0.4]
        IF KEYWORD_SET(color) THEN BEGIN
            if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
            symbols = fltarr(n) + 4
            symsizes = fltarr(n) + 2
        ENDIF ELSE BEGIN
            colors = fltarr(n) + fgcolor
            symbols = [4,4,6,6] ;(findgen(n)+2)*2
            symsizes = [1,2,1,2]
        ENDELSE
        velocities_true = [56,61,53.2281,56.8234]
        velocities = velocities_true ;taken from flat part of rotation curve
        velocitiesHI = [54.552800,61.797691,53.406370,54.337259] ;from HI line widths
        
        gh516_H_lr = {B:-16.373721563452,  V:-16.6365793369207, alpha:2.1449371, alpha_gal:2.52505, alphaV:2.1781738, alpha_galV:2.56192} ;h516.cosmo25cmb.3072g14HBWK
        gh516_M_lr = {B:-15.7881139871463, V:-16.2747647221014, alpha:1.20930,   alpha_gal:1.20545, alphaV:1.0868989, alpha_galV:1.20930}
        gh516_H    = {B:-15.6996847198888, V:-16.0784155581285, alpha:0.81438442,alpha_gal:1.02964, alphaV:0.00000001,alpha_galV:0.000001} ;h516.cosmo25cmb.3072g14HBWK
        gh516_M    = {B:-15.2652157189731, V:-15.7550592733375, alpha:0.79430061,alpha_gal:0.770795,alphaV:0.79689074,alpha_galV:0.777515}
        simdata = [gh516_M_lr,gh516_H_lr,gh516_M,gh516_H]
        obssymsize = 2
        maxdistance = 8
    END
    3: BEGIN
      prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
        steps = ['00512','00512','00512','00512']
        dir = prefix + ['h516.cosmo25cmb.768g/h516.cosmo25cmb.768g14HBWK/steps/h516.cosmo25cmb.768g14HBWK.'   + steps[0] + '.dir',$
                        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.'   + steps[1] + '.dir',$
                        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.' + steps[2] + '.dir',$
                        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[3] + '.dir']
        files = prefix + ['h516.cosmo25cmb.768g/h516.cosmo25cmb.768g14HBWK/steps/h516.cosmo25cmb.768g14HBWK.'   + steps[0] + '.dir/h516.cosmo25cmb.768g14HBWK.'  + steps[0],$
                          'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.'   + steps[1] + '.dir/h516.cosmo25cmb.1536g14HBWK.'  + steps[1],$
                         'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.' + steps[2] + '.dir/h516.cosmo25cmb.2304g14HBWK.' + steps[2],$
                          'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[3] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[3]]
        pfiles = prefix + ['h516.cosmo25cmb.768g/h516.cosmo25cmb.768g14HBWK/h516.cosmo25cmb.768g14HBWK.param',$
                           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.param',$
                           'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/h516.cosmo25cmb.2304g14HBWK.param',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
        sflog = prefix +  ['h516.cosmo25cmb.768g/h516.cosmo25cmb.768g14HBWK/h516.cosmo25cmb.768g14HBWK.starlog',$
                           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.starlog',$
                          'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/h516.cosmo25cmb.2304g14HBWK.starlog',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.starlog']
        filebases = ['h516.cosmo25cmb.768g14HBWK.'  + steps[0] + '.halo.1',$
                     'h516.cosmo25cmb.1536g14HBWK.'  + steps[1] + '.halo.1',$
                     'h516.cosmo25cmb.2304g14HBWK.' + steps[2] + '.halo.1',$
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[3] + '.halo.1']
        filenames = ['h516.cosmo25cmb.768g14HBWK.'  + steps[0],$
                     'h516.cosmo25cmb.1536g14HBWK.'  + steps[1], $
                     'h516.cosmo25cmb.2304g14HBWK.' + steps[2],$
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[3]]
        outext = ['H2, 768',$
                  'H2, 1536' ,$
                  'H2, 2304' ,$
                  'H2, 3072']
        n = N_ELEMENTS(files)
        useH2 = [1,1,1,1]
        halos = [1,1,1,1]
        halos_str = strtrim(halos,2)
        keys = ['DH2, ulr','DH2 lr','DH2 mr','DH2, hr']
;        psym = [17,6]
        distunits = [25000.,25000.,25000.,25000.]
        massunits  = [2.310e15,2.310e15,2.310e15,2.310e15]
        linestyles = [0,0,0,0]
;        thicks = [1,2,4,6]
        filternums = [13,13,13,13] ;[13,14,13]
        cameras = [14,14,14,14] ;[14,15,14]
        yrange_vcirc = [0,75]
        yrange_photo = [26,19]
        maxdistance_photo = 4
        yrange_SFH = [0,0.4]
        IF KEYWORD_SET(color) THEN BEGIN
            if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
            symbols = fltarr(n) + 4
            symsizes = fltarr(n) + 2
        ENDIF ELSE BEGIN
            colors = fltarr(n) + fgcolor
            symbols = [4,4,4,4] ;(findgen(n)+2)*2
            symsizes = [1,2,3,4]
        ENDELSE
        obssymsize = 2
        maxdistance = 8
    END
    4: BEGIN
        prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
        steps = ['00492','00492','00512','00512','00512']
        dir = prefix + ['h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g2MBWK/steps/h516.cosmo25cmb.2304g2MBWK.'   + steps[0] + '.dir',$
                        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[1] + '.dir',$
                        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.' + steps[2] + '.dir',$
                        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g16HBWK/steps/h516.cosmo25cmb.2304g16HBWK.' + steps[3] + '.dir',$
                        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[4] + '.dir']
        files = prefix + ['h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g2MBWK/steps/h516.cosmo25cmb.2304g2MBWK.'   + steps[0] + '.dir/h516.cosmo25cmb.2304g2MBWK.'  + steps[0],$
                          'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[1] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[1],$
                          'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.' + steps[2] + '.dir/h516.cosmo25cmb.1536g14HBWK.' + steps[2],$
                          'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g16HBWK/steps/h516.cosmo25cmb.2304g16HBWK.' + steps[3] + '.dir/h516.cosmo25cmb.2304g16HBWK.' + steps[3],$
                          'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[4] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[4]]
        pfiles = prefix + ['h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g2MBWK/h516.cosmo25cmb.2304g2MBWK.param',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.param',$
                           'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g16HBWK/h516.cosmo25cmb.2304g16HBWK.param',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
        sflog = prefix +  ['h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g2MBWK/h516.cosmo25cmb.2304g2MBWK.starlog',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.starlog',$
                           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.starlog',$
                           'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g16HBWK/h516.cosmo25cmb.2304g16HBWK.starlog',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.starlog']
        filebases = ['h516.cosmo25cmb.2304g2MBWK.'  + steps[0] + '.halo.1',$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[1] + '.halo.1',$
                     'h516.cosmo25cmb.1536g14HBWK.' + steps[2] + '.halo.1',$
                     'h516.cosmo25cmb.2304g16HBWK.' + steps[3] + '.halo.1',$
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[4] + '.halo.1']
        filenames = ['h516.cosmo25cmb.2304g2MBWK.'  + steps[0],$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[1], $
                     'h516.cosmo25cmb.1536g14HBWK.' + steps[2],$
                     'h516.cosmo25cmb.2304g16HBWK.' + steps[3],$
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[4]]
        outext = ['met, 1536',$
                  'met, 3072',$
                  'H2, 1536' ,$
                  'H2, 2304, NO H2 SF' ,$                  
                  'H2, 3072']
        n = N_ELEMENTS(files)
        useH2 = [0,0,1,1,1]
        halos = [1,1,1,1,1]
        halos_str = strtrim(halos,2)
        keys = ['DnoH2, lr','DnoH2','DH2 lr','DH2, no SF','DH2']
;        psym = [17,6]
        distunits = [25000.,25000.,25000.,25000.,25000.]
        massunits  = [2.310e15,2.310e15,2.310e15,2.310e15,2.310e15]
        linestyles = [2,2,0,3,0]
        thicks = [1,2,1,1,2]
        filternums = [13,13,13,13,13] ;[13,14,13]
        cameras = [14,14,14,14,14] ;[14,15,14]
        yrange_vcirc = [0,75]
        yrange_photo = [26,19]
        maxdistance_photo = 4
        yrange_SFH = [0,0.4]
        IF KEYWORD_SET(color) THEN BEGIN
            if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
            symbols = fltarr(n) + 4
            symsizes = fltarr(n) + 2
        ENDIF ELSE BEGIN
            colors = fltarr(n) + fgcolor
            symbols = [4,4,6,8,6] ;(findgen(n)+2)*2
            symsizes = [1,2,1,1,2]
        ENDELSE
        velocities_true = [56,61,53.2281,56.8234]
        velocities = velocities_true ;taken from flat part of rotation curve
        velocitiesHI = [54.552800,61.797691,53.406370,54.337259] ;from HI line widths

        gh516_H_lr1 = {B:-15.4908047304537, V:-15.8834548099003, alpha:1.4043764, alpha_gal:0.615625, alphaV:2.1781738, alpha_galV:2.56192} ;h516.cosmo25cmb.3072g14HBWK        
        gh516_H_lr =  {B:-16.373721563452,  V:-16.6365793369207, alpha:2.0300389, alpha_gal:2.56192 , alphaV:2.1781738, alpha_galV:2.56192} ;h516.cosmo25cmb.3072g14HBWK
        gh516_M_lr =  {B:-15.7881139871463, V:-16.2747647221014, alpha:1.20930,   alpha_gal:1.20545, alphaV:1.0868989, alpha_galV:1.20930}
        gh516_H    =  {B:-15.6996847198888, V:-16.0784155581285, alpha:1.6039808, alpha_gal:1.62768, alphaV:0.00000001,alpha_galV:0.000001} ;h516.cosmo25cmb.3072g14HBWK
        gh516_M    =  {B:-15.2652157189731, V:-15.7550592733375, alpha:1.1288082, alpha_gal:0.770795,alphaV:0.79689074,alpha_galV:0.777515}
        simdata = [gh516_M_lr,gh516_M,gh516_H_lr,gh516_H_lr1,gh516_H]
        obssymsize = 2
        maxdistance = 8
    END
   5: BEGIN
        prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
        steps = ['00492','00512','00512','00512']
        dir = prefix +  ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir',$
                         'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir',$
                         'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.' + steps[2] + '.dir',$
                         'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.' + steps[3] + '.dir']

        files = prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir/h516.cosmo25cmb.3072g1MBWK.'   + steps[0],$
                         'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[1],$
                         'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.' + steps[2] + '.dir/h516.cosmo25cmb.2304g14HBWK.' + steps[2],$
                         'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.' + steps[3] + '.dir/h516.cosmo25cmb.1536g14HBWK.' + steps[3]]
        pfiles = prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param',$
                           'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/h516.cosmo25cmb.2304g14HBWK.param',$
                           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.param']
        sflog = prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.starlog',$
                           'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.starlog',$
                           'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/h516.cosmo25cmb.2304g14HBWK.starlog',$
                           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.starlog']

        filebases = prefix +  ['h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.halo.1',$
                               'h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.halo.1',$
                               'h516.cosmo25cmb.2304g14HBWK.' + steps[2] + '.halo.1',$
                               'h516.cosmo25cmb.1536g14HBWK.' + steps[3] + '.halo.1']
        filenames = prefix +  ['h516.cosmo25cmb.3072g1MBWK.'   + steps[0],$
                               'h516.cosmo25cmb.3072g14HBWK.' + steps[1],$
                               'h516.cosmo25cmb.2304g14HBWK.' + steps[2],$
                               'h516.cosmo25cmb.1536g14HBWK.' + steps[3]]
        outext = ['met, 3072',$
                  'H2, 3072',$
                  'H2, 2304',$                  
                  'H2, 1536']

        n = N_ELEMENTS(files)
        useH2 = [0,1,1,1]
        halos = [1,1,1,1]
        halos_str = strtrim(halos,2)
        keys = ['DnoH2','DH2','DH2_mr','DH2_lr']
;        psym = [17,6]
        distunits = [25000.,25000.,25000.,25000.]
        massunits  = [2.310e15,2.310e15,2.310e15,2.310e15]
        linestyles = [2,0,0,0]
        thicks = [5,5,3,1]
        filternums = [13,13,13,13] ;[13,14,13]
        cameras = [14,14,14,14] ;[14,15,14]
        yrange_vcirc = [0,75]
        yrange_photo = [26,19]
        maxdistance_photo = 4
        yrange_SFH = [0,0.4]
        IF KEYWORD_SET(color) THEN BEGIN
            if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
            symbols = fltarr(n) + 4
            symsizes = fltarr(n) + 2
        ENDIF ELSE BEGIN
            colors = fltarr(n) + fgcolor
            symbols = [4,6,6,6] ;(findgen(n)+2)*2
            symsizes = [3,3,2,1]
        ENDELSE
        velocities_true = [56,61,53.2281,56.8234]
        velocities = velocities_true ;taken from flat part of rotation curve
        velocitiesHI = [54.552800,61.797691,53.406370,54.337259] ;from HI line widths

        gh516_H_lr =  {B:-16.373721563452,  V:-16.6365793369207, alpha:2.0300389, alpha_gal:2.56192 , alphaV:2.1781738, alpha_galV:2.56192} ;h516.cosmo25cmb.3072g14HBWK
        gh516_H_mr =  {B:-15.7265738971342, V:-16.0888994327437, alpha:1.5484495, alpha_gal:1.79816, alphaV:0.0000001, alpha_galV:0.0000001}
        gh516_H    =  {B:-15.6996847198888, V:-16.0784155581285, alpha:0.96378927, alpha_gal:1.18978, alphaV:0.00000001,alpha_galV:0.000001} ;h516.cosmo25cmb.3072g14HBWK      ;1.6039808  1.62768
        gh516_M    =  {B:-15.2652157189731, V:-15.7550592733375, alpha:0.93222844, alpha_gal:1.14012,alphaV:0.79689074,alpha_galV:0.777515} ;4
;1.1288082  1.48637
        simdata = [gh516_M,gh516_H,gh516_H_mr,gh516_H_lr]
        obssymsize = 2
        maxdistance = 8
    END
ENDCASE
formatplot,outplot = outplot;,/thick

;zsol = 0.0130215
;ind = where(g.tempg lt 2e4)
;ind = where(g.dens*rhounit ge 0.1)
;histogramp,g[ind].zmetal/zsol,nbins = 100
;print,mean(g[ind].zmetal/zsol)

; Halo Mags -------------------------------------------------------------------------------------------------------
IF plot_get_halo_mags THEN BEGIN
    FOR i = 0, n - 1 DO BEGIN
        cd,dir[i]
        get_halo_mags,massunits[i],distunits[i],mag_sys = 'ab',/multiple
    ENDFOR
ENDIF


;-----------------------------------------------------------------------------------
;------------ Cooling Rate -----------------
; ~/code/IDL/MolecH/coolrate.pro
IF plot_coolrate THEN BEGIN
    prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/'
    cfiles = prefix +['h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.out/h516.cosmo25cmb.3072g14HBWK.halo.1','h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.out/h516.cosmo25cmb.3072g1MBWK.00444.halo.1','h516.cosmo25cmb.3072g15HBWK/steps/h516.cosmo25cmb.3072g15HBWK.00460.dir/h516.cosmo25cmb.3072g15HBWK.00460.halo.1']
    distunits1 = [25000.,25000.,25000.]
    massunits1  = [2.310e15,2.310e15,2.310e15]
    If KEYWORD_SET(outplot) THEN coolrate,cfiles,distunits1,massunits1,outplot = outplot ELSE coolrate,cfiles,distunits1,massunits1
ENDIF


;------------------------------------------------------------------------------------
;------------ Resolution/Surface Density Tests ----------
; ~/code/IDL/MolecH/checkAbundMet.pro
;Dangerous
IF plot_checkAbundMet THEN BEGIN
    if keyword_set(color) THEN BEGIN
        If KEYWORD_SET(outplot) THEN  checkAbundMet_master,outplot = outplot,/color ELSE checkAbundMet_master,/color 
    endif ELSE begin
        If KEYWORD_SET(outplot) THEN  checkAbundMet_master,outplot = outplot ELSE checkAbundMet_master
    endelse

endif

;------------------------------------------------------------------------------------
;------------ Resolution/Surface Density Tests ----------
; ~/code/IDL/MolecH/checkAbundRes.pro
;Dangerous
IF plot_checkAbundMet THEN BEGIN
    If KEYWORD_SET(outplot) THEN  checkAbundRes_master,/outfile ELSE checkAbundRes_master
endif

;-------------------------------------------------------------------------------------

;------------ Resolved K-S law data ----------------------
;~/code/HIcubes/schmidtlaw_res_obs.pro
IF plot_schmidtlaw_res_obs THEN BEGIN
    res = 0.75
    FOR i = 0, n - 1 DO BEGIN
        cd,dir[i]
        schmidtlaw_res_obs,filenames[i],pfiles[i],res = res,useH2 = useH2[i],outplot = outplot,angle = angle[i],extno = extno[i],verbose = verbose,rotateAngle = rotateAngle[i]
    ENDFOR
ENDIF

;------------- Resolved K-S law plot ----------------------
;~/code/HIcubes/schmidtlaw_res_obs.pro
skfilename = '/schmidtlaw_res_obs0.750000.dat'
;skfilename = '/schmidtlaw_res_obs_all_Ha0.750000.dat'
IF plot_schmidtlaw_res_obs_master_out THEN BEGIN
    IF KEYWORD_SET(color) THEN schmidtlaw_res_obs_master_out,dir + skfilename,outplot = outplot,color = colors,thick = thicks_psym,symbols = symbols,key = keys,symsize = symsize,/scaleHe,/lowz ELSE schmidtlaw_res_obs_master_out,dir + skfilename,outplot = outplot,thick = thicks_psym,symbols = symbols,key = keys,symsize = symsize,/scaleHe,/lowz
ENDIF

;------------- Global K-S Law -----------------------------
;/astro/users/christensen/code/IDL/HIcubes/schmidtlaw_global_obs.pro
IF plot_schmidtlaw_global_obs THEN BEGIN
    data=fltarr(2,n)
    FOR i = 0, n - 1 DO BEGIN
        data[*,i] = schmidtlaw_global_obs(dir[i],filenames[i],pfiles[i],useH2 = useH2[i],Halpha=1,tipsyfile = filebases[i],verbose = verbose,halo_str = halos_str[i],angle = angle[i],camera = extno[i],intq = intq[i])
    ENDFOR      

;   IF NOT KEYWORD_SET(color) THEN loadct,0
    IF (KEYWORD_SET(outplot)) THEN BEGIN
        set_plot,'ps'
        device,filename = outplot + '_globSchmidt.eps',/color,bits_per_pixel= 8,/times,ysize=18,xsize=18,xoffset =  2,yoffset =  2 
    ENDIF ELSE window,0
    
    xsigma = 10.0^(findgen(600)/100 - 1.) ;xsigmalow
    ysigma=2.5e-4*xsigma^1.4

    readcol,'~/code/Datafiles/HIcubes/ks98.dat',name,D,logHI,logH2,logH,logSFR,tdyn,co,HI,halpha ;,format='(A9D)'
    readcol,'~/code/Datafiles/HIcubes/uavpl.dat',logHdarf,logSFRdwarf
    plot,alog10(xsigma),alog10(ysigma),ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-0.5,2.5], yrange = [-4,-0.5], xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]"   ; -0.5
    oplot,logH,logSFR,psym = obssym,color = obscolor,symsize = 2
    oplot,logHdarf,logSFRdwarf,psym = 1,color = obscolor,symsize = 2 
    FOR i = 0, N_ELEMENTS(dir) - 1 DO oplot,[alog10(data[0,i]),alog10(data[0,i])],[alog10(data[1,i]),alog10(data[1,i])],psym = symbols[i],color = colors[i],symsize = symsizes[i],thick = thicks_psym[i]
    oplot,[-0.2,0.2],[-1,-1]
    oplot,[0,0],[-0.8,-1.2]
    legend,[keys,'Kennicutt 1998','Roychowdhury et al. 2009'],psym = [symbols,obssym,1],color = [colors,obscolor,obscolor],/bottom,/right,thick = [thicks_psym,1,1];,symsize = [symsizes,2,2]
    IF (KEYWORD_SET(outplot)) THEN  device,/close  ELSE stop
ENDIF

;------------- color, scale, and mag -----------------------------

; here
IF plot_magcolor THEN BEGIN
    datafile = '~/code/Datafiles/Sunrise_analysis/vanZee_short.txt'
    readcol,datafile,name1,name2,D25,d_25,I,PA,m_B,m_B_error,A_B,B_V, B_V_error,U_B,U_B_error,D,MB,Sigma_25,mu0B,mu0cB,alpha_kpc,alpha_arcsec,Halpha,Halpha_error,LF,format='(A,A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)' 


    IF (KEYWORD_SET(outplot)) THEN device,filename = outplot+'_Magcolor.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2  ELSE window,0,xsize = 712,ysize = 450
;    IF (KEYWORD_SET(outplot)) THEN device,filename = outplot+'_Magalphacolor.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 15,xoffset =  2,yoffset = 2  ELSE window,0,xsize = 712,ysize = 450
;    multiplot,[1,2],/doyaxis
    plot,MB,B_V,psym = obssym,ytitle = textoidl('B-V'),yrange=[0.201,0.8],ystyle = 1,symsize = obssymsize,xtitle = textoidl('M_B')
    oplot,MB,B_V,psym = obssym,color = obscolor,symsize = obssymsize
    FOR i=0,n -1 DO oplot,[simdata[i].B,simdata[i].B],[simdata[i].B - simdata[i].V, simdata[i].B - simdata[i].V],psym = symbols[i],color = colors[i],symsize = symsizes[i],thick = thicks_psym[i]
    IF KEYWORD_SET(keys) THEN legend,[keys,"van Zee 2000"],color = [colors,obscolor],psym = [symbols,obssym],/right,thick = [thicks_psym,1];symsize = [2,2,symsizes]
;    multiplot,/doyaxis
;    plot,MB,alog10(alpha_kpc),psym = obssym,xtitle = textoidl('M_B'),ytitle = textoidl('log \alpha [kpc]'),yrange = [-0.7,0.7],symsize = obssymsize
;    oplot,MB,alog10(alpha_kpc),psym = obssym,color = obscolor,symsize = obssymsize
;    FOR i=0,n -1 DO oplot,[simdata[i].B,simdata[i].B],alog10([simdata[i].alpha_gal,simdata[i].alpha_gal]),psym = symbols[i],color = colors[i],symsize = symsizes[i]
;    IF KEYWORD_SET(keys) THEN legend,['van Zee `00',keys],color = [obscolor,colors],psym = [obssym,symbols],symsize = [2,2,symsizes],/right
;    multiplot,/reset

    IF (KEYWORD_SET(outplot)) THEN BEGIN
        device,/close
        device,filename = outplot+'_Scalecolor.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
    ENDIF ELSE window,1,xsize = 712,ysize = 392

    plot,B_V,alog10(alpha_kpc),psym = obssym,xtitle = textoidl('B - V'),ytitle = textoidl('log \alpha [kcp]'),symsize = obssymsize ;,yrange=[0.2,0.8]
    oplot,B_V,alog10(alpha_kpc),psym = obssym,color = obscolor,symsize = obssymsize
    FOR i=0,n -1 DO oplot,[simdata[i].B - simdata[i].V, simdata[i].B - simdata[i].V],alog10([simdata[i].alpha_gal,simdata[i].alpha_gal]),psym = symbols[i],color = colors[i],symsize = symsizes[i]
    IF KEYWORD_SET(keys) THEN legend,['van Zee `00',keys],color = [obscolor,colors],psym = [obssym,symbols],symsize = [2,2,symsizes]

    if (KEYWORD_SET(outplot)) then device,/close else stop
ENDIF

;------------- Tully-Fisher -----------------------------
;/astro/users/christensen/code/IDL/HIcubes/tully_fisher_obs.pro
IF plot_tully_fisher_obs THEN BEGIN
    gmass = tully_fisher_obs_gasmass(files+'.halo.1', massunits, smass = smass, gmassall = gmassall)
    print,'True Gas Mass:     ',gmass
    print,'Total Gas Mass:    ',gmassall
    print,'True Stellar Mass: ',smass
    IF KEYWORD_SET(color) THEN tully_fisher_obs,files + '.1/broadband.fits', velocitiesHI, gmass, key = keys,outfile = outplot, symbols = symbols, symsizes = symsizes,color = colors,thicks = thicks_psym ELSE tully_fisher_obs,files + '.1/broadband.fits', velocitiesHI, gmass, key = keys,outfile = outplot, symbols = symbols, symsizes = symsizes,thicks = thicks_psym
ENDIF

;------------- SFH ----------------------
;/astro/users/christensen/code/IDL/MolecH/sfh.pro
IF plot_sfh THEN BEGIN
;    sfh,files+ '.halo.1',pfiles,outplot = outplot,key = keys,binsize = 1e7,yrange = [0,0.05],xrange = [0,1]
;stop
;    sfh,files+ '.halo.1',pfiles,outplot = outplot,key = keys,binsize = 5e7,yrange = [0,0.4],colors = colors,thick = thicks,linestyle = linestyles

;    sfh,files+ '.halo.1',pfiles,outplot = outplot,keys = keys,binsize = 5e7,colors = colors,thick = thicks,linestyle = linestyles,/cumlative,yrange = [1e5,2.7e8],/ylog,/redshift
   sfh_z,dir,files+'.halo.1',pfiles,keys = keys,colors = colors,linestyle = linestyles,outplot = outplot,thick = thicks,binsize = 5e7,/cumlative,yrange = [1e5,3e8],xrange = [0,13.7346],/ylog,/redshift,formatthick = formatthick
ENDIF

;------------- Color Prof ----------------------
;/astro/users/christensen/code/IDL/Sunrise_analysis/color_prof.pro
IF plot_color_prof THEN BEGIN
    IF KEYWORD_SET(color) $
      THEN color_prof, files + '.1/broadband.fits',key = keys,filt = filternums,dir = files + '.1',/r25,outplot = outplot,linestyle = linestyles,thick = thicks,maxdistance = maxdistance,color = colors $
      ELSE color_prof, files + '.1/broadband.fits',key = keys,filt = filternums,dir = files + '.1',/r25,outplot = outplot,linestyle = linestyles,thick = thicks,maxdistance = maxdistance
    FOR i = 0, n - 1 DO BEGIN
        cd,files[i] + '.1'
        filters= mrdfits('broadband.fits',filternums[i])
        print,'Mag_{B}: ',strtrim(filters[21].AB_MAG0,2),', B-V: ',strtrim(filters[21].AB_MAG0 - filters[22].AB_MAG0,2)
        rmax = opticalRadii(extno = cameras[i],B_num = 21,/verbose)
        print,'R_{25}: ',strtrim(rmax,2)    
    ENDFOR
ENDIF

;-------------- toomreRadius_Master--------------
;/astro/users/christensen/code/IDL/HIcubes/toomreRadius.pro
IF plot_toomreRadius THEN BEGIN
    if keyword_set(outplot) THEN BEGIN 
        IF KEYWORD_SET(color) THEN toomreRadius,files,distunits,massunits,outplot = outplot,maxdistance = maxdistance_photo,thicks = thicks,symbols = symbols,keys = keys,color = colors ELSE toomreRadius,files,distunits,massunits,outplot = outplot,maxdistance = maxdistance_photo,thicks = thicks,symbols = symbols,keys = keys
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(color) THEN toomreRadius,files,distunits,massunits,maxdistance = maxdistance_photo,thicks = thicks,symbols = symbols,keys = keys,color = colors ELSE toomreRadius,files,distunits,massunits,maxdistance = maxdistance_photo,thicks = thicks,symbols = symbols,keys = keys
    ENDELSE
ENDIF

;---------------- FFT -----------------
;/astro/users/christensen/code/IDL/HIcubes/fft_result.pro
IF plot_fft_result THEN BEGIN
    if keyword_set(outplot) THEN BEGIN
        IF KEYWORD_SET(color) THEN fft_result,dir,key = keys,thick = thicks,symbols = symbols,outplot = outplot,color = colors ELSE fft_result,dir,key = keys,thick = thicks,symbols = symbols,outplot = outplot
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(color) THEN fft_result,dir,key = keys,thick = thicks,symbols = symbols,color = colors ELSE fft_result,dir,key = keys,thick = thicks,symbols = symbols
    ENDELSE
ENDIF

;-------------- densRadius --------------
;/astro/users/christensen/code/IDL/HIcubes/densRadius.pro
;Dangerous
IF plot_densRadius THEN BEGIN
    if keyword_set(outplot) THEN BEGIN
        IF KEYWORD_SET(color) THEN densRadius,files,distunits,massunits,maxdistance = maxdistance_photo,outplot = outplot,color = colors ELSE densRadius,files,distunits,massunits,maxdistance = maxdistance_photo,outplot = outplot
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(color) THEN densRadius,files,distunits,massunits,maxdistance = maxdistance_photo,color = colors ELSE densRadius,files,distunits,massunits,maxdistance = maxdistance_photo
    ENDELSE
ENDIF

;------------- SF gas ----------------------
;/astro/users/christensen/code/IDL/MolecH/sfgas.pro
IF plot_sfgas THEN BEGIN
    densconvert = massunits[0] * gm_per_msol * amu_per_gm/distunits[0]^3/cm_per_kpc^3
    if keyword_set(outplot) THEN BEGIN
        IF KEYWORD_SET(color) THEN sfgas,sflog,dens_convert = densconvert,key = keys,thick = thicks,linestyle = linestyles,molecularH = useH2,outplot = outplot,color = colors,yrange = [0,0.035] ELSE sfgas,sflog,dens_convert = densconvert,key = keys,thick = thicks,linestyle = linestyles,molecularH = useH2,outplot = outplot,yrange = [0,0.035]
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(color) THEN sfgas,sflog,dens_convert = densconvert,key = keys,thick = thicks,linestyle = linestyles,molecularH = useH2,color = colors,yrange = [0,0.035] ELSE sfgas,sflog,dens_convert = densconvert,key = keys,thick = thicks,linestyle = linestyles,molecularH = useH2,yrange = [0,0.035]
    ENDELSE
ENDIF

;---------------------------- Picture -------------------
IF plot_make_jpeg THEN BEGIN
    cameraFO = cameras
    cameraEO = cameras + 2
;    IF NOT KEYWORD_SET(outplot) THEN outfile = files ELSE 
    IF KEYWORD_SET(outplot) THEN BEGIN
        outfile = outplot + outext
        FOR i = 0, n - 1 DO BEGIN
            print,files[i]
            cd,files[i] + '.' + halos_str[i]
;        center_temp = [-1.0*center[0,i],center[1,i]] 
            scales = [1.5,4.2,1.2]*0.5
            make_jpeg,bands = [4,3,2],scales = scales,cam = cameraFO[i],outfile = outfile[i] + '_FO.jpeg',range = maxdistance*2.0 ;range
            make_jpeg,bands = [4,3,2],scales = scales,cam = cameraEO[i],outfile = outfile[i] + '_EO.jpeg',range = maxdistance*2.0 ;,rotateAngle = rotateAngle,center = center_temp,range = 24;range
        ENDFOR
    ENDIF
ENDIF

;------------- Mass-Metallicity -----------------------------
;/astro/users/christensen/code/IDL/procedures/mzr.pro
IF plot_mzr THEN BEGIN
    IF KEYWORD_SET(color) THEN mzr_plot, files, halos = [1,1], key = keys, outfile = outplot, psym = symbols, symsize = symsizes, color = colors, thick = thicks, xrange = [5.5, 11], yrange = [7, 9],/readfile ELSE mzr_plot, files, halos = [1,1], key = keys, outfile = outplot, psym = symbols, symsize = symsizes, thick = thicks, xrange = [5.5, 11], yrange = [7, 9],/readfile
ENDIF

stop
END
