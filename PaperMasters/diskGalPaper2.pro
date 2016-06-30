;.r tully_fisher_obs
;.r matchHalos
;.r opticalRadii
;.r Halpha
;.r schmidtlaw_global_obs
;.r densRadius
;.r writeHImacro
;.r multiHaloGuo
;.r photometricProf
;.r diskGalPaper

;outplot = '~/plots/h603/h603.cosmo50cmb.00512'
;outplot = '~/plots/h603/h986.cosmo50cmb.00512'
;outplot = '~/plots/h603/h516.cosmo25cmb.00512'
;outplot = '~/plots/h606/hall.00512_color'

pro diskGalPaper,outplot = outplot,color = color,verbose = verbose
idldir = '/astro/users/christensen/code/IDL/'
;resolve_routine,idldir + 'HIcubes/tully_fisher_obs.pro'
;resolve_routine,idldir + 'idl4tipsy/matchHalos.pro',/is_function
;resolve_routine,idldir + 'Sunrise_analysis/photometricProf.pro'
;resolve_routine,idldir + 'HIcubes/schmidtlaw_global_obs.pro',/is_function

; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL

plot_get_halo_mags = 0 ;pre
plot_matchHalos = 0 ;pre
plot_twinsanalysis = 0
plot_multiHaloGuo = 0
plot_HImacro = 0
plot_analyzeHIcubes = 0
plot_vcirc = 1
plot_photometricProf = 0
plot_make_jpeg = 0
plot_tully_fisher_obs = 0
plot_densRadius = 0;dangerous
plot_sfh = 0
plot_schmidtlaw_global_obs = 0
plot_schmidtlaw_res_obs = 0 ;pre
plot_schmidtlaw_res_obs_master_out = 0
plot_fft_result = 0

;1: h603.72
;2: h603.328
;3: h603.512
;4: h986.72
;5: h986:512
;6: h516:72
;7: h516:512
;8: h603 (3), h986(3), h516(2)
;9: h603 (2), h986(3), h516(2)
x=9

CASE x OF
    1: BEGIN
        steps = ['00072','00072','00072' ]
        dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.' + steps[0] + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir']
        files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.' + steps[0] + '.dir/h603.cosmo50cmb.2304g2bwK.' + steps[0] ,$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[2]]
        pfiles =  ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param']
        filebases = ['h603.cosmo50cmb.2304g2bwK.' + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.halo.1']
        filenames = ['h603.cosmo50cmb.2304g2bwK.' + steps[0],$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2]]
        outext = ['',$
                  'met',$
                  'H2']
        useH2 = [0,0,1]
        halos = [1,1,1]
        halos_str = strtrim(halos,2)
;keys   = ['no Metals','Metals, no H2', 'H2']
        psym = [4,6,8]
        distunits = [50000.,50000.,50000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16]
        linestyles = [3,2,0]
        filternums = [13,13,13] ;[13,14,13]
        cameras = [14,14,14]    ;[14,15,14]
        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        yrange_SFH = [0,7]
        haloids = [[1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1],$
                   [1,2,3,4,5,6,7,8, 9,10,11,12,13,14,15,16,17],$
                   [1,2,3,5,6,4,7,8,10,13,11,14,12, 9,15,17,16]]
        center = [[0,0],$
                  [0,0],$
                  [0,0]]
    END
    2: BEGIN
        steps = ['00328','00324','00324']
;steps = ['00512','00512','00324']
        dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.' + steps[0]     + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir']
        files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.'     + steps[0] + '.dir/h603.cosmo50cmb.2304g2bwK.'   + steps[0],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[2]]
        pfiles =  ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param']
        filebases = ['h603.cosmo50cmb.2304g2bwK.'     + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.halo.1']
        filenames = ['h603.cosmo50cmb.2304g2bwK.'   + steps[0],$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2]]
        outext = ['orig',$
                  'met',$
                  'H2']
        useH2 = [0,0,1]
        halos = [1,1,1]
        halos_str = strtrim(halos,2)
;keys   = ['no Metals','Metals, no H2', 'H2']
        psym = [4,6,8]
        distunits = [50000.,50000.,50000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16]
        linestyles = [3,2,0]
        filternums = [13,13,13] ;[13,14,13]
        cameras = [14,14,14]    ;[14,15,14]
        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        yrange_SFH = [0,7]
;orig satellites:       17, 12                
;satellities (metal): 
;at least 20,000 DM particles
;    haloids = [[1,2,3,4,5,6,15,16,18,20],$ ;Check this first line
;               [1,2,3,4,5,6,15,16,18,20],$
;               [1,2,3,4,5,7,14,16,18,19]]
        haloids = [[1, 2, 3, 6, 17, 12, 15, 16, 25, 44,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1], $
                   [1, 2, 3, 4,  5,  6,  7,  8,  9, 15, 16, 18, 20, 21, 23, 26, 30, 31, 33, 34, 35, 36, 38, 39, 40, 41, 42, 43, 44, 51, 52],$
                   [1, 2, 3, 4,  5,  7,  6,  8,  9, 14, 16, 18, 19, 22, 24, 46, 31, 30, 33, 36, 34, 35, 39, 40, 41, 38, 43, 44, 49, 76, 53]]
        center = [[0,0],$
                  [0,0],$
                  [0,0]]
    END
    3: BEGIN
        steps = ['00512','00512','00512']
;steps = ['00512','00512','00324']
        dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.' + steps[0]     + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir']
        files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.'     + steps[0] + '.dir/h603.cosmo50cmb.2304g2bwK.'   + steps[0],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[2]]
        pfiles =  ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param']
        filebases = ['h603.cosmo50cmb.2304g2bwK.'     + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.halo.1']
        filenames = ['h603.cosmo50cmb.2304g2bwK.'   + steps[0],$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2]]
        outext = ['orig',$
                  'met',$
                  'H2']
        useH2 = [0,0,1]
        halos = [1,1,1]
        halos_str = strtrim(halos,2)
;keys   = ['no Metals','Metals, no H2', 'H2']
        psym = [4,6,8]
        distunits = [50000.,50000.,50000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16]
        linestyles = [3,2,0]
        filternums = [13,13,13] ;[13,14,13]
        cameras = [14,14,14]    ;[14,15,14]
        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        yrange_SFH = [0,7]
;orig satellites:       17, 12                
;satellities (metal): 
;at least 20,000 DM particles
;    haloids = [[1,2,3,4,5,6,15,16,18,20],$ ;Check this first line
;               [1,2,3,4,5,6,15,16,18,20],$
;               [1,2,3,4,5,7,14,16,18,19]]
;        haloids = [[1, 2, 3, 6, 17, 12, 15, 16, 25, 44,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1], $
;                   [1, 2, 3, 4,  5,  6,  7,  8,  9, 15, 16, 18, 20, 21, 23, 26, 30, 31, 33, 34, 35, 36, 38, 39, 40, 41, 42, 43, 44, 51, 52],$
;                   [1, 2, 3, 4,  5,  7,  6,  8,  9, 14, 16, 18, 19, 22, 24, 46, 31, 30, 33, 36, 34, 35, 39, 40, 41, 38, 43, 44, 49, 76, 53]]
        haloids = [[1,2,3,1,1,1],$
                   [1,2,3,7,8,16],$
                   [1,2,3,8,7,13]]
                   
        center = [[0,0],$
                  [0,0],$
                  [0,0]]
    END
    4: BEGIN
        steps = ['00072','00072']
        dir = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[0],$
;               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir',$
               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir']
       files = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '/h986.cosmo50cmb.3072gs1MbwK.' + steps[0],$
;               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir/h986.cosmo50cmb.3072gs1MbwK.' + steps[0],$
               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[1]]
        pfiles = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param']

        filebases = ['h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[1] + '.halo.1']
        filenames = ['h986.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[1]]
        outext = ['met',$
                  'H2']
        useH2 = [0,1]
        halos = [1,1]
        halos_str = strtrim(halos,2)
;keys   = ['no Metals','Metals, no H2', 'H2']
        psym = [6,8]
        distunits = [50000.,50000.]
        massunits  = [1.84793e16,1.84793e16]
        linestyles = [2,0]
        filternums = [13,13] ;[13,14,13]
        cameras = [16,14]    ;[14,15,14]
        intq = [38,33]
        recenter = [1,0]
        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        yrange_SFH = [0,7]
;orig satellites:       17, 12                
;satellities (metal): 
;at least 20,000 DM particles
        haloids = [[1,2,3,4,5,6,7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21],$
                   [1,2,3,4,5,7,6,15,11,10,20,13,17,16,19, 9,18,14,22,30,22]]

        center = [[0,0],$
                  [0,0]]
    END
    5: BEGIN
        steps = ['00512','00512','00512']
        dir = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[0],$
               '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[1],$
;               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir',$
               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir']
       files = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[0] + '/h986.cosmo50cmb.3072g1bwK.' + steps[0],$
               '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[1] + '/h986.cosmo50cmb.3072gs1MbwK.' + steps[1],$
;               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir/h986.cosmo50cmb.3072gs1MbwK.' + steps[0],$
               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[2]]
        pfiles = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.param',$
                  '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param']

        filebases = ['h986.cosmo50cmb.3072g1bwK.'     + steps[0] + '.halo.1',$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[2] + '.halo.1']
        filenames = ['h986.cosmo50cmb.3072g1bwK.'   + steps[0],$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[2]]
        outext = ['orig',$
                  'met',$
                  'H2']
        useH2 = [0,0,1]
        halos = [1,1,1]
        halos_str = strtrim(halos,2)
;keys   = ['no Metals','Metals, no H2', 'H2']
        psym = [4,6,8]
        distunits = [50000.,50000.,50000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16]
        linestyles = [3,2,0]
        filternums = [13,13,13] ;[13,14,13]
        cameras = [14,16,14]    ;[14,15,14]
        intq = [38,38,33]
        recenter = [0,1,0]
        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        yrange_SFH = [0,7]
;orig satellites:       17, 12                
;satellities (metal): 
;at least 20,000 DM particles
        haloids = [[1, 2, 4, 9, 11, 13, 14],$
                   [1, 2, 3, 4,  5,  7, 10],$
                   [1, 2, 3, 9, 11, 14, 17]]

        center = [[0,0],$
                  [-8,72],$
                  [0,0]]
    END
    6: BEGIN
        steps = ['00072','00080']
;steps = ['00512','00512','00324']
        dir =    ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir',$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir']
        files =  ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[0],$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[1]]
        pfiles = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
        filebases = ['h516.cosmo25cmb.3072g1MBWK.'  + steps[0] + '.halo.1', $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.halo.1']
        filenames = ['h516.cosmo25cmb.3072g1MBWK.'  + steps[0], $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[1]]
        outext = ['met',$
                  'H2']
        useH2 = [0,1]
        halos = [1,1]
        halos_str = strtrim(halos,2)
;keys   = ['Metals, no H2', 'H2']
        psym = [6,8]
        distunits = [25000.,25000.]
        massunits  = [2.310e15,2.310e15]
        linestyles = [2,0]
        filternums = [13,13]    ;[13,14,13]
        cameras = [14,14]       ;[14,15,14]
        yrange_vcirc = [0,75]
        yrange_photo = [26,19]
        maxdistance_photo = 4
        maxdistance = 6
        yrange_SFH = [0,0.4]       
        haloids = [[1,2,3,4,6],$
                   [2,1,4,3,6]]                   
        center = [[0,0],$
                  [0,0]]
    END
    7: BEGIN
        steps = ['00492','00512']
;steps = ['00512','00512','00324']
        dir =    ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir',$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir']
        files =  ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[0],$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[1]]
        pfiles = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
        filebases = ['h516.cosmo25cmb.3072g1MBWK.'  + steps[0] + '.halo.1', $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.halo.1']
        filenames = ['h516.cosmo25cmb.3072g1MBWK.'  + steps[0], $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[1]]
        outext = ['met',$
                  'H2']
        useH2 = [0,1]
        halos = [1,1]
        halos_str = strtrim(halos,2)
;keys   = ['Metals, no H2', 'H2']
        psym = [6,8]
        distunits = [25000.,25000.]
        massunits  = [2.310e15,2.310e15]
        linestyles = [2,0]
        filternums = [13,13]    ;[13,14,13]
        cameras = [14,14]       ;[14,15,14]
        yrange_vcirc = [0,75]
        yrange_photo = [26,19]
        maxdistance_photo = 4
        maxdistance = 6
        yrange_SFH = [0,0.4]       
        haloids = [[1,2,3,6,13],$
                   [1,2,3,6,13]]                   
        center = [[0,0],$
                  [0,0]]
    END
    8: BEGIN
       steps = ['00512','00512','00512','00512','00512','00512','00492','00512']
        dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.' + steps[0]     + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir',$
                 '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[3],$
                 '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[4],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[6] + '.dir',$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[7] + '.dir']

        files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.'     + steps[0] + '.dir/h603.cosmo50cmb.2304g2bwK.'   + steps[0],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[2],$
                 '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[3] + '/h986.cosmo50cmb.3072g1bwK.' + steps[3],$
               '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[4] + '/h986.cosmo50cmb.3072gs1MbwK.' + steps[4],$
               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[5],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[6] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[6],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[7] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[7]]

        pfiles =  ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
                   '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.param',$
                  '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']

        filebases = ['h603.cosmo50cmb.2304g2bwK.'     + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.halo.1',$
                     'h986.cosmo50cmb.3072g1bwK.'     + steps[3] + '.halo.1',$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[4] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.halo.1',$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[6] + '.halo.1', $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[7] + '.halo.1']

        filenames = ['h603.cosmo50cmb.2304g2bwK.'   + steps[0],$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2],$
                     'h986.cosmo50cmb.3072g1bwK.'   + steps[3],$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[4],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[5],$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[6], $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[7]]

        outext = ['orig',$
                  'met',$
                  'H2',$
                  'orig',$
                  'met',$
                  'H2',$
                  'met',$
                  'H2']

        useH2 = [0,0,1,0,0,1,0,1]
        halos = [1,1,1,1,1,1,1,1]
        halos_str = strtrim(halos,2)
;keys   = ['no Metals','Metals, no H2', 'H2']
        psym = [4,6,8,4,6,8,6,8]
        distunits = [50000.,50000.,50000.,50000.,50000.,50000.,25000.,25000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16,1.84793e16,1.84793e16,1.84793e16,2.310e15,2.310e15]
        linestyles = [3,2,0,3,2,0,2,0]
        filternums = [13,13,13,13,13,13,13,13]
        cameras =    [14,14,14,14,16,14,14,14]
        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        yrange_SFH = [0,7]
        haloids = [[1],[1],[1],[1],[1],[1],[1],[1]]
        intq = [38,38,33,33,38,33,32,32]
        center = [[0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [-8,72],$
                  [0,0],$
                  [0,0],$
                  [0,0]]
        vfinal =  [131.610, 120.751, 115.436, 96.2218, 109.535, 101.855, 52.5520, 56.3969]
    END
    9: BEGIN
       steps = ['00512','00512','00512','00512','00512','00512','00492','00512']
        dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir',$
                 '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[2],$
                 '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[3],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[4] + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/steps/h516.cosmo25cmb.3072g3BWK.'     + steps[5] + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[6] + '.dir',$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[7] + '.dir']

        files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[1],$
                 '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[2] + '/h986.cosmo50cmb.3072g1bwK.' + steps[2],$
                 '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[3] + '/h986.cosmo50cmb.3072gs1MbwK.' + steps[3],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[4] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[4],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/steps/h516.cosmo25cmb.3072g3BWK.'     + steps[5] + '.dir/h516.cosmo25cmb.3072g3BWK.'   + steps[5],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[6] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[6],$
                 '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[7] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[7]]

        pfiles =  ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
                   '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.param',$
                  '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/h516.cosmo25cmb.3072g3bwK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                   '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']

        filebases = ['h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.halo.1',$
                     'h986.cosmo50cmb.3072g1bwK.'   + steps[2] + '.halo.1',$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[3] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[4] + '.halo.1',$
                     'h516.cosmo25cmb.3072g3BWK.'   + steps[5] + '.halo.1',$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[6] + '.halo.1', $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[7] + '.halo.1']

        filenames = ['h603.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[1],$
                     'h986.cosmo50cmb.3072g1bwK.'   + steps[2],$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[3],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[4],$
                     'h516.cosmo25cmb.3072g3BWK.'   + steps[5],$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[6], $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[7]]

        outext = ['met',$
                  'H2',$
                  'CH',$
                  'met',$
                  'H2',$
                  'CH',$
                  'met',$
                  'H2']

        useH2 = [0,1,0,0,1,0,0,1]
        halos = [1,1,1,1,1,1,1,1]
        halos_str = strtrim(halos,2)
;keys   = ['no Metals','Metals, no H2', 'H2']
        psym = [6,8,4,6,8,4,6,8]
        distunits = [50000.,50000.,50000.,50000.,50000.,25000.,25000.,25000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16,1.84793e16,1.84793e16,2.310e15,2.310e15,2.310e15]
        linestyles = [2,0,3,2,0,3,2,0]
        filternums = [13,13,13,13,13,13,13,13]
        cameras =    [14,14,14,16,14,14,14,14]
        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        yrange_SFH = [0,7]
        haloids = [[1],[1],[1],[1],[1],[1],[1],[1]]
        intq = [38,33,38,38,33,32,32,32]
        center = [[0,0],$
                  [0,0],$
                  [0,0],$
                  [-8,72],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0]]
;        vfinal =  [120.751, 115.436, 96.2218, 109.535, 101.855, 52.5520, 56.3969]
        colors = [50,50,130,130,130,240,240,240]
        END
ENDCASE
n = N_ELEMENTS(files)
a = fltarr(n)
FOR i = 0, n - 1 DO BEGIN
    rtipsy,files[i],h,g,d,s,/justhead
    a[i] = h.time
ENDFOR

formatplot,outplot = outplot
IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    obscolor = 0 ;fgcolor
    if colors[0] eq 1 then  colors  = (findgen(n) + 1)*240/n else colors = colors
    IF NOT KEYWORD_SET(psym) THEN psym = fltarr(n) + 4
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
ENDIF ELSE BEGIN
    loadct,0    
    obscolor = 100
    colors = (findgen(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(broadband)) + 5
    IF NOT KEYWORD_SET(psym) THEN  psym = (findgen(n)+2)*2
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2 ;thicks = (findgen(n) + 1)*6/n - 1
    IF NOT KEYWORD_SET(linestyles) THEN linestyles = REVERSE(findgen(n)*2) 
ENDELSE

; Halo Mags -------------------------------------------------------------------------------------------------------
IF plot_get_halo_mags THEN BEGIN
    FOR i = 0, n - 1 DO BEGIN
        cd,dir[i]
        get_halo_mags,massunits[i],distunits[i],mag_sys = 'ab',/multiple
    ENDFOR
ENDIF

;matchHalos -----------------------------------------------------------------------------------------------------------------------------------
IF plot_matchHalos THEN BEGIN
;    FOR i = 1, n - 1 DO BEGIN
;        print,matchHalos([files[0] + '.amiga.stat',files[i] + '.amiga.stat'],minmass = 10^(8.5))
;    ENDFOR
    matches0 = matchHalos([files[1] + '.amiga.stat',files[0] + '.amiga.stat'],minmass = 10^(8.5))
    IF n gt 2 THEN matches1 = matchHalos([files[1] + '.amiga.stat',files[2] + '.amiga.stat'],minmass = 10^(8.5))
ENDIF

;--------------------------- Twins analysis -------------------------
IF plot_twinsanalysis THEN BEGIN
    twinsanalysis,outplot = outplot
ENDIF

;multiHaloGuo ----------------------------------------------------------------------------------------------------------------------------------
IF plot_multiHaloGuo THEN BEGIN
    IF KEYWORD_SEt(outplot) THEN multiHaloGuo,files,key = keys,yrange = [4,10.5],haloids = haloids,halomasses = halomasses,psym = psym,xrange = [8.5,12],outplot = outplot + '_MstarMhalo.eps' $
      ELSE multiHaloGuo,files,key = keys,yrange = [4,10.5],haloids = haloids,halomasses = halomasses,psym = psym,xrange = [8.5,12]

;   IF KEYWORD_SET(outplot) THEN device, filename= outplot + '_MstarMhaloRatio.eps',bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch ELSE stop
;   plot,alog10(halomasses[*,0]),alog10(halomasses[*,1]),psym = 5,xtitle = 'Log(Mass Halo)',ytitle = 'Log(Mstar_H2/Mstar_{no H2})',yrange = [-1,1],xrange = [9,11.5]
;   oplot,alog10(halomasses[*,0]),alog10(halomasses[*,1])
;   IF KEYWORD_SET(outplot) THEN device,/close
    stop
    IF KEYWORD_SEt(outplot) THEN multiHaloGuo,files,key = keys,yrange = [4,10.5],haloids = haloids,halomasses = halomasses,psym = psym,xrange = [8.5,12],/kcorrect,outplot = outplot + '_MstarMhalo_Obs.eps' $
      ELSE multiHaloGuo,files,key = keys,yrange = [4,10.5],haloids = haloids,halomasses = halomasses,psym = psym,xrange = [8.5,12],/kcorrec
ENDIF

;HImacro --------------------------------------------------------------------------------------------------------------------------------------
IF plot_HImacro THEN BEGIN
    cd,dir[0]
    writeHImacro,filebases[0],dist_units = distunits[0],mass_unit = massunits[0],radius = 24.0
    cd,dir[1]
    writeHImacro,filebases[1],dist_units = distunits[1],mass_unit = massunits[1],radius = 24.0
    cd,dir[2]
    writeHImacro,filebases[2],dist_units = distunits[2],mass_unit = massunits[2],radius = 24.0,/molecularH
ENDIF

;Smooth cube and moment maps --------------------------------------------------------------------------------------------------------------------------------------
IF plot_analyzeHIcubes THEN BEGIN
    analyzeHIcubes,dir[0],filebases[0],/physical_coord
    analyzeHIcubes,dir[1],filebases[1],/physical_coord
    analyzeHIcubes,dir[2],filebases[2],/physical_coord
ENDIF

;--------------------------- Vcirc -------------
IF 0 THEN BEGIN
IF plot_vcirc THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,outfile = outplot + '_vcirc.eps',maxdistance = 15,vfinal = vfinal,color = colors ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,outfile = outplot + '_vcirc.eps',maxdistance = 15,vfinal = vfinal
    ENDIF ELSE BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = 15,vfinal = vfinal,color = colors ELSE  vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = 15,vfinal = vfinal
    ENDELSE
ENDIF
ENDIF

IF plot_vcirc THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,outfile = outplot + '_vcircdark.eps',maxdistance = 15,vfinal = vfinal,type = 'dark',color = colors ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,outfile = outplot + '_vcircdark.eps',maxdistance = 15,vfinal = vfinal,type = 'dark'
    ENDIF ELSE BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = 15,vfinal = vfinal,type = 'dark',color = colors ELSE  vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = 15,vfinal = vfinal,type = 'dark'
    ENDELSE
    print,vfinal
ENDIF

IF plot_vcirc THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,outfile = outplot + '_vcircstar.eps',maxdistance = 15,vfinal = vfinal,type = 'star',color = colors ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,outfile = outplot + '_vcircstar.eps',maxdistance = 15,vfinal = vfinal,type = 'star'
    ENDIF ELSE BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = 15,vfinal = vfinal,type = 'star',color = colors ELSE  vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = 15,vfinal = vfinal,type = 'star'
    ENDELSE
ENDIF

IF plot_vcirc THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,outfile = outplot + '_vcircgas.eps',maxdistance = 15,vfinal = vfinal,type = 'gas',color = colors ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,outfile = outplot + '_vcircgas.eps',maxdistance = 15,vfinal = vfinal,type = 'gas'
    ENDIF ELSE BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = 15,vfinal = vfinal,type = 'gas',color = colors ELSE  vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = 15,vfinal = vfinal,type = 'gas'
    ENDELSE
ENDIF


;--------------------------- Photometric Profile ------------
IF plot_photometricProf THEN BEGIN
    IF KEYWORD_SET(outplot) THEN  photometricProf,files + '.' + halos_str + '/broadband.fits',keys = keys,thicks = thicks,linestyle = linestyles,filternums = filternums,cameras = cameras,yrange = yrange_photo,outplot = outplot + '_photoProf.eps',bands = [5,5,5],maxdistance = maxdistance_photo, recenter = recenter $
      ELSE photometricProf,files + '.' + halos_str + '/broadband.fits',keys = keys,thicks = thicks,linestyle = linestyles,filternums = filternums,cameras = cameras,yrange = yrange_photo,bands = [5,5,5],maxdistance = maxdistance_photo,recenter = recenter
ENDIF

;---------------------------- Picture -------------------
IF plot_make_jpeg THEN BEGIN
    cameraFO = cameras
    cameraEO = cameras + 2
    IF NOT KEYWORD_SET(outplot) THEN outfile = files ELSE outfile = outplot + outext
    FOR i = 0, n - 1 DO BEGIN
        cd,files[i] + '.' + halos_str[i]
        make_jpeg,bands = [4,3,2],scales = [5.0,5.0,5.0],cam = cameraFO[i],outfile = outfile[i] + '_FO.jpeg'
        make_jpeg,bands = [4,3,2],scales = [5.0,5.0,5.0],cam = cameraEO[i],outfile = outfile[i] + '_EO.jpeg'
    ENDFOR
ENDIF

;-------------------------------- TF ----------------------
IF plot_tully_fisher_obs THEN BEGIN
    vmax = fltarr(n)
    vmaxobs = fltarr(n)
    head = ' '    
    FOR i = 0, n - 1 DO BEGIN
        data = read_stat_struc_amiga(files[i] + '.amiga.stat')
        ind = where(data.group eq halos[i])
        vmax[i] = data[ind].vc
        openr,1,dir[i] + '/linewidths.txt'
        readf,1,head,format = '(A)'
        readf,1,width,width2,vmaxline
        close,1
        vmaxobs[i] = vmaxline ;width2
    ENDFOR
    gmass = tully_fisher_obs_gasmass(files + '.halo.' + halos_str, massunits,smass = smass)
 ;   IF KEYWORD_SET(outplot) THEN tully_fisher_obs,files + '.' + halos_str + '/broadband.fits',vmax,gmass,filternums = filternums,symbols = psym, key = keys, outfile = outplot + '_vtrue_' $
 ;   ELSE BEGIN
 ;       tully_fisher_obs,files + '.' + halos_str + '/broadband.fits',vmax,gmass,filternums = filternums,symbols = psym, key = keys
 ;       stop
 ;   ENDELSE
    IF KEYWORD_SET(outplot) THEN tully_fisher_obs,files + '.' + halos_str + '/broadband.fits',vmaxobs,gmass,filternums = filternums,symbols = psym, key = keys, outfile = outplot + '_vobs_' $
    ELSE BEGIN
        tully_fisher_obs,files + '.' + halos_str + '/broadband.fits',vmaxobs,gmass,filternums = filternums,symbols = psym, key = keys
        stop
    ENDELSE
;    IF KEYWORD_SET(outplot) THEN tully_fisher_obs_btf, files + '.' + halos_str + '/broadband.fits', vmaxobs, gmass, key = keys,smass_true = smass,velocities_true = vfinal,filternums = filternums,symbols = psym, outfile = outplot + '_vAll_' $
;    ELSE BEGIN
;        tully_fisher_obs_btf, files + '.' + halos_str + '/broadband.fits', vmaxobs, gmass, key = keys,smass_true = smass,velocities_true = vfinal,filternums = filternums,symbols = psym
;        stop
;    ENDELSE
    IF KEYWORD_SET(outplot) THEN tully_fisher_obs_btf, files + '.' + halos_str + '/broadband.fits', vfinal, gmass, key = keys,filternums = filternums,symbols = psym, outfile = outplot + '' $
    ELSE BEGIN
        tully_fisher_obs_btf, files + '.' + halos_str + '/broadband.fits', vfinal, gmass, key = keys,filternums = filternums,symbols = psym
        stop
    ENDELSE
ENDIF

;----------------------------------- Density Radius ----------
IF plot_densRadius THEN BEGIN
    densRadius,files  + '.halo.' + halos_str,distunits,massunits, outplot = outplot
ENDIF

;------------------------------------ SFH ----------------------
IF plot_sfh THEN BEGIN
    IF KEYWORD_SET(outplot) THEN sfh,files  + '.halo.' + halos_str,pfiles,key = keys,linestyle = linestyles,yrange = yrange_SFH, outplot = outplot + '_SFH.eps' ELSE sfh,files  + '.halo.' + halos_str,pfiles,key = keys,linestyle = linestyles,yrange = yrange_SFH
ENDIF

;-------------------------- Global SK law ----------------------------
IF  plot_schmidtlaw_global_obs THEN BEGIN
    data=fltarr(2,N_ELEMENTS(dir))
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
        print,dir[i]
        data[*,i] = schmidtlaw_global_obs(dir[i],filenames[i],pfiles[i],useH2 = useH2[i],tipsyfile = filebases[i],/Halpha,intq = intq[i], camera = cameras[i],center = center[*,i],verbose = verbose)
    ENDFOR 

    if (KEYWORD_SET(outplot)) then begin
        device,filename = outplot + '_SK.eps',/color,bits_per_pixel= 8,/times,ysize=18,xsize=18
    endif else begin
        set_plot,'x'
        window,1
    endelse
    xsigma = 10.0^(findgen(600)/100 - 1.) ;xsigmalow
    ysigma=2.5e-4*xsigma^1.4
    
    readcol,'~/code/HIcubes/ks98.dat',name,D,logHI,logH2,logH,logSFR,tdyn,co,HI,halpha ;,format='(A9D)'
    plot,alog10(xsigma),alog10(ysigma),ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-0.5,2.5], yrange = [-4,-0.5], xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]"   
    oplot,logH,logSFR,psym = 2,color = obscolor
    FOR i = 0, N_ELEMENTS(dir) - 1 DO oplot,[alog10(data[0,i]),alog10(data[0,i])],[alog10(data[1,i]),alog10(data[1,i])],psym = psym[i],color = colors[i] 
    oplot,[-0.2,0.2],[-1,-1]
    oplot,[0,0],[-0.8,-1.2]
    IF KEYWORD_SET(keys) THEN legend,[keys,'Kennicutt 98'],psym = [psym,2],color = [colors,obscolor],/bottom,/right
    if (KEYWORD_SET(outplot)) then device,/close else stop
ENDIF

;-------------------------- Res SK law ----------------------------
IF plot_schmidtlaw_res_obs THEN BEGIN
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
        cd,dir[i]
        schmidtlaw_res_obs,filenames[i],pfiles[i],useH2 = useH2[i]
        stop
    ENDFOR
ENDIF


;------------- Resolved K-S law plot ----------------------
;~/code/HIcubes/schmidtlaw_res_obs.pro
IF plot_schmidtlaw_res_obs_master_out THEN BEGIN
    IF KEYWORD_SET(color) THEN schmidtlaw_res_obs_master_out,outplot = outplot,color = colors,thick = thicks ELSE schmidtlaw_res_obs_master_out,outplot = outplot,thick = thicks
ENDIF

;---------------- FFT -----------------
;/astro/users/christensen/code/IDL/HIcubes/fft_result.pro
IF plot_fft_result THEN BEGIN
    if keyword_set(outplot) THEN BEGIN
        IF KEYWORD_SET(color) THEN fft_result,dir,key = keys,thick = thicks,symbols = psym,outplot = outplot,color = colors ELSE fft_result,dir,key = keys,thick = thicks,symbols = psym,outplot = outplot
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(color) THEN fft_result,dir,key = keys,thick = thicks,symbols = psym,color = colors ELSE fft_result,dir,key = keys,thick = thicks,symbols = psym
    ENDELSE
ENDIF

stop
end
