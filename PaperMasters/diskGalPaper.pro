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
;outplot = '~/plots/h603/h799.cosmo25cmb.00512'
;outplot = '~/plots/h603/hall.00512'

;outplot = '~/plots/h603/h516.cosmo25cmb.00512.fabio'
;outplot = '~/plots/h603/hall.00512.fabio'

;outplot = '~/plots/h603/h603.cosmo50cmb.00512.dashed'
;outplot = '~/plots/h603/h986.cosmo50cmb.00512.dashed'
;outplot = '~/plots/h603/h516.cosmo25cmb.00512.dashed'
;outplot = '~/plots/h603/h799.cosmo25cmb.00512.dashed'
;outplot = '~/plots/h603/hall.00512.dashed'

;diskGalPaper,/color,/formatthick,outplot = '~/plots/h603/h986.cosmo50cmb.00512.talk'

;3651884408
pro diskGalPaper,outplot = outplot,color = color,verbose = verbose,formatthick = formatthick
formatplot,outplot = outplot,thick = formatthick
IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(outplot) THEN bgcolor = 255 ELSE bgcolor = 0

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
IF hostname EQ 'ozma' THEN datafile_prefix = '/home/christensen/Code/Datafiles/' ELSE datafile_prefix = '/astro/users/christec/code/Datafiles/'

;resolve_routine,'tully_fisher_obs'
;resolve_routine,'matchHalos'
;resolve_routine,'photometricProf'
;resolve_routine,'schmidtlaw_global_obs',/is_function
;resolve_routine,'Halpha',/is_funtion
;resolve_routine,'opticalRadii',/is_funtion
idldir = '/astro/users/christensen/code/IDL/'
file_compile,idldir + 'HIcubes/tully_fisher_obs.pro'
file_compile,idldir + 'idl4tipsy/matchHalos.pro';,/is_function
file_compile,idldir + 'Sunrise_analysis/photometricProf.pro'
file_compile,idldir + 'HIcubes/opticalRadii.pro'
file_compile,idldir + 'HIcubes/Halpha.pro'
file_compile,idldir + 'HIcubes/schmidtlaw_global_obs.pro';,/is_function

; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL

plot_get_halo_mags = 0 ;pre
plot_matchHalos = 0 ;pre
plot_HImacro = 0 ;pre
plot_analyzeHIcubes = 0 ;pre
plot_twinsanalysis = 0
plot_multiHaloGuo = 0

plot_make_jpeg = 0
plot_gas_surface_den = 0;1

plot_schmidtlaw_res_obs = 0 ;pre
plot_vcirc = 0 ;pre
plot_DMprof = 0
plot_photometricProf = 0
plot_fft_result = 0
plot_j_test = 0
plot_densRadius = 0;dangerous
;-----------------------------
plot_sfh = 0
plot_track_mass = 0
plot_sigmar = 0
plot_fbrad = 0
plot_fbz = 0
plot_fbcum = 0
plot_sfgas = 0

plot_schmidtlaw_global_obs = 1
plot_schmidtlaw_res_obs_master_out = 0
plot_tully_fisher_obs = 0

;1: h603.72
;2: h603.328
;3: h603.512 (3)X
;4: h603.512 (2) 
;5: h986.72
;6: h986:512 (3) X
;7: h516:72
;8: h516:512 (2)
;9: h516:492 (3) X
;10: h603 (3), h986(3), h516(2)
;11: h603 (2), h986(3), h516(2)
;12: h603 (3), h986(3), h516(3) 
;13: h603 (2), h986(3), h516(2)
;14: h986:512 (2)
;15: h799 (3) X
;16: h603 (3), h986(3), h516(3), h799(3) X
;17: h603(1), h986(2) X
x = 16
CASE x OF
    1: BEGIN
        steps = ['00072','00072','00072' ]
        dir   =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.' + steps[0] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir']
        files =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.' + steps[0] + '.dir/h603.cosmo50cmb.2304g2bwK.' + steps[0] ,$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[2]]
        pfiles =  [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param']
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
        maxdistance = 8
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
        dir   =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.' + steps[0]     + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir']
        files =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.'     + steps[0] + '.dir/h603.cosmo50cmb.2304g2bwK.'   + steps[0],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[2]]
        pfiles =  [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param']
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
        maxdistance = 8
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
        dir   =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g1bwK/steps/h603.cosmo50cmb.3072g1bwK.' + steps[0]     + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir']
        files =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g1bwK/steps/h603.cosmo50cmb.3072g1bwK.'     + steps[0] + '.dir/h603.cosmo50cmb.3072g1bwK.'   + steps[0],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[2]]
        pfiles =  [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g1bwK/h603.cosmo50cmb.3072g1bwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param']
        filebases = ['h603.cosmo50cmb.3072g1bwK.'     + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.halo.1']
        filenames = ['h603.cosmo50cmb.3072g1bwK.'   + steps[0],$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2]]
        outext = ['warm',$
                  'met',$
                  'H2']
        useH2 = [0,0,1]
;        useH2 = [0,0,0]
        halos = [1,1,1]
        halos_str = strtrim(halos,2)
        distunits = [50000.,50000.,50000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16]

        filternums = [13,13,13] ;[13,14,13]
;FO        cameras = [14,14,14]    ;[14,15,14]
;        rotateAngle = [3,3,3]       
        cameras = [15,15,15] ;45
        rotateAngle = [5,5,5]  
        recenter = [0,0,0]
        center = [[0,0],$
                  [0,0],$
                  [0,0]]
        angle = [45.0,45.0,45.0]
        intq = [33,33,33]     
        resizefactor = [2,2,2]
        bands = [5,5,5]

        yrange_vcirc = [0,250]
        yrange_vcirc = [0,350]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        maxdistance = 8 ;10 
        maxdistance = 10
        yrange_SFH = [0,7]
        yrange_fb = [0,1.2]
        xrange_fb = [0,8.0]
;        normalize = 1
;        yrange_fb = [0,0.03]
        yrange_fbcum = [0,0.2];12]
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
        vfinal = [ 102.159  ,    122.812  ,    111.321]
             
        colors = [120,fgcolor,254]      
        linestyles = [2,0,0]
        IF KEYWORD_SET(outplot) THEN thicks = [5,5,8] ELSE thicks = [1,1,3]
;        linestyles = [3,2,0]
;        psym = [4,6,8]
;keys   = ['no Metals','Metals, no H2', 'H2']
        psym = [4,5,16]
        ctables = [0,0,39]
        keys   = ['no Metals, no ' + textoidl('H_2'),'Metals, no ' + textoidl('H_2'),  textoidl('H_2')]
        label = 'h603'
    END
    4: BEGIN
        steps = ['00512','00512']
;steps = ['00512','00512','00324']
        dir   =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir']
        files =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[1]]
        pfiles =  [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param']
        filebases = ['h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.halo.1']
        filenames = ['h603.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[1]]
        outext = ['met',$
                  'H2']
        useH2 = [0,1]
        halos = [1,1]
        halos_str = strtrim(halos,2)
;        keys   = ['no Metals, no ' + textoidl('H_2'),'Metals, no ' + textoidl('H_2'),  textoidl('H_2')]
        psym = [6,8]
        distunits = [50000.,50000.]
        massunits  = [1.84793e16,1.84793e16]
;        linestyles = [2,0]
        colors = [80,254]
        filternums = [13,13] ;[13,14,13]
        cameras = [14,14]    ;[14,15,14]
        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        maxdistance = 8
        yrange_SFH = [0,7]
        xrange_fb = [0,10]
        yrange_fb = [0,4.0]
        yrange_fbcum = [0,0.15]
       
;orig satellites:       17, 12                
;satellities (metal): 
;at least 20,000 DM particles
;    haloids = [[1,2,3,4,5,6,15,16,18,20],$ ;Check this first line
;               [1,2,3,4,5,6,15,16,18,20],$
;               [1,2,3,4,5,7,14,16,18,19]]
;        haloids = [[1, 2, 3, 6, 17, 12, 15, 16, 25, 44,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1], $
;                   [1, 2, 3, 4,  5,  6,  7,  8,  9, 15, 16, 18, 20, 21, 23, 26, 30, 31, 33, 34, 35, 36, 38, 39, 40, 41, 42, 43, 44, 51, 52],$
;                   [1, 2, 3, 4,  5,  7,  6,  8,  9, 14, 16, 18, 19, 22, 24, 46, 31, 30, 33, 36, 34, 35, 39, 40, 41, 38, 43, 44, 49, 76, 53]]
        haloids = [[1,2,3,7,8,16],$
                   [1,2,3,8,7,13]]
                   
        center = [[0,0],$
                  [0,0],$
                  [0,0]]
        rotateAngle = [3,3,3]
        label = 'h603'
    END
    5: BEGIN
        steps = ['00072','00072']
        dir =   ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[0],$
;                 prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir']
       files =   ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '/h986.cosmo50cmb.3072gs1MbwK.' + steps[0],$
;                 prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir/h986.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[1]]
        pfiles = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param']

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
        maxdistance = 8        
        yrange_SFH = [0,7]
;orig satellites:       17, 12                
;satellities (metal): 
;at least 20,000 DM particles
        haloids = [[1,2,3,4,5,6,7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21],$
                   [1,2,3,4,5,7,6,15,11,10,20,13,17,16,19, 9,18,14,22,30,22]]

        center = [[0,0],$
                  [0,0]]
    END
    6: BEGIN
        steps = ['00512','00512','00512']
        dir =    ['/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[0],$
;                 '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir']
       files =   ['/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[0] + '/h986.cosmo50cmb.3072g1bwK.' + steps[0],$
;                 '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[1] + '/h986.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h986.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[2]]
        pfiles = ['/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.param',$
                 ; '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param']
        filebases = ['h986.cosmo50cmb.3072g1bwK.'   + steps[0] + '.halo.1',$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[2] + '.halo.1']
        filenames = ['h986.cosmo50cmb.3072g1bwK.'   + steps[0],$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[2]]
        outext = ['warm',$
                  'met',$
                  'H2']
        useH2 = [0,0,1]
        halos = [1,1,1]
        halos_str = strtrim(halos,2)
        distunits = [50000.,50000.,50000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16]

        filternums = [13,13,13] ;[13,14,13]
;        cameras = [14,16,14]    ;[14,15,14]
;        rotateAngle = [3,5,3]
;        angle = [90.0,90.0,90.0]
;        intq = [38,38,33]
;FO        cameras = [14,16,14]    ;[14,15,14]
;        rotateAngle = [3,3,3]       
        cameras = [15,15,15] ;45
        rotateAngle = [5,5,5] 
        recenter = [0,0,0]
        center = [[0,0],$
                  [0,0],$
                  [0,0]]
        angle = [45.0,45.0,45.0]
        intq = [38,33,33]
        resizefactor = [2,2,2]
        bands        = [7,7,7]

        yrange_vcirc = [0,250]
;        yrange_vcirc = [0,350]
        yrange_photo = [26,15]
        maxdistance_photo = 8 ;16
        maxdistance = 8;10
        maxdistance = 10
        maxdistance_dens = 16
;        yrange_SFH = [0,7]
        yrange_SFH = [0,3.5]
        xrange_fb = [0,8]
        yrange_fb = [0,1.2]
;        normalize = 1
;        yrange_fb = [0,0.03]
        yrange_fbcum = [0,0.3]

;orig satellites:       17, 12                
;satellities (metal): 
;at least 20,000 DM particles
        haloids = [[1, 2, 4, 9, 11, 13, 14],$
                   [1, 2, 3, 4,  5,  7, 10],$
                   [1, 2, 3, 9, 11, 14, 17]]
        vfinal = [96.2218, 109.535, 101.855]

        colors = [120,fgcolor,254]      
;        colors = [50,120,245]
        IF KEYWORD_SET(outplot) THEN thicks = [5,5,8] ELSE thicks = [1,1,3]
        linestyles = [2,0,0]
;        linestyles = [0,0,0]
        ctables = [0,0,39]
;        ctables = [39,39,39]
        psym = [4,6,16]
        keys   = ['Primordial','Metals',textoidl('H_2')]
;        keys   = ['no Metals, no ' + textoidl('H_2'),'Metals, no ' + textoidl('H_2'),  textoidl('H_2')] 
       label = 'h986'
       normalize = 1
    END
    7: BEGIN
        steps = ['00072','00080']
;steps = ['00512','00512','00324']
        dir =    [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir']
        files =  [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[0],$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[1]]
        pfiles = [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
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
        maxdistance = 8
        yrange_SFH = [0,0.4]       
        haloids = [[1,2,3,4,6],$
                   [2,1,4,3,6]]                   
        center = [[0,0],$
                  [0,0]]
    END
    8: BEGIN
        steps = ['00492','00512']
;steps = ['00512','00512','00324']
        dir =    [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir']
        files =  [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[0] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[0],$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[1]]
        pfiles = [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
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
        maxdistance = 8
        yrange_SFH = [0,0.4]       
        haloids = [[1,2,3,6,13],$
                   [1,2,3,6,13]]                   
        center = [[0,0],$
                  [0,0]]
    END
    9: BEGIN
 ;       steps = ['00512','00492','00512']
        steps = ['00492','00492','00492']
;steps = ['00512','00512','00324']
        dir =    [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/steps/h516.cosmo25cmb.3072g3BWK.'     + steps[0] + '.dir',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[1] + '.dir',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[2] + '.dir']
        files =  [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/steps/h516.cosmo25cmb.3072g3BWK.'     + steps[0] + '.dir/h516.cosmo25cmb.3072g3BWK.'   + steps[0],$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[1] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[1],$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[2] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[2]]
        pfiles = [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/h516.cosmo25cmb.3072g3bwK.param',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
        filebases = ['h516.cosmo25cmb.3072g3BWK.'   + steps[0] + '.halo.1',$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[1] + '.halo.1', $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[2] + '.halo.1']
        filenames = ['h516.cosmo25cmb.3072g3BWK.'   + steps[0],$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[1], $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[2]]
        outext = ['warm',$
                  'met',$
                  'H2']
        useH2 = [0,0,1]
        halos = [1,1,1]
        halos_str = strtrim(halos,2)
        distunits = [25000.,25000.,25000.]
        massunits  = [2.310e15,2.310e15,2.310e15]

        filternums = [13,13,13]    ;[13,14,13]
;FO        cameras = [14,14,14]       ;[14,15,14]
;        rotateAngle = [3,3,3]
        cameras = [15,15,15]       ;[14,15,14]
        rotateAngle = [5,5,5]
        recenter =   [0,0,0]
        center = [[0,0],$
                  [0,0],$
                  [0,0]]        
        angle = [45.0,45.0,45.0]
;        intq = [33,32,32]
        intq= [33,32,33]
        resizefactor = [2,2,2]
        bands        = [5,5,5]

        yrange_vcirc = [0,75]
        yrange_photo = [26,19]
        maxdistance_photo = 4
        maxdistance = 8 ;10
        yrange_SFH = [0,0.4]
        xrange_fb = [0,4]
;        yrange_fb = [0,5e7]     
;        normalize = 1  
        yrange_fb = [0,0.2]
        yrange_fbcum = [0,0.15]

        haloids = [[1,2,4,6,16],$
                   [1,2,3,6,13],$
                   [1,2,3,6,13]]                   
        vfinal = [57.2425,52.5520,56.3969]

        colors = [120,fgcolor,254]
        IF KEYWORD_SET(outplot) THEN thicks = [5,5,8] ELSE thicks = [1,1,3]
;        linestyles = [3,2,0]
        linestyles = [2,0,0]
        psym = [4,6,16]
        ctables = [0,0,39]
;keys   = ['Metals, no H2', 'H2']
 ;       keys   = ['no Metals, no ' + textoidl('H_2'),'Metals, no ' + textoidl('H_2'),  textoidl('H_2')]
        label = 'h516'
    END
   10: BEGIN
       steps = ['00512','00512','00512','00512','00512','00512','00492','00512']
        dir   =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.' + steps[0]     + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir',$
                   '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[3],$
                   '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[4],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[6] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[7] + '.dir']

        files =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.'     + steps[0] + '.dir/h603.cosmo50cmb.2304g2bwK.'   + steps[0],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[2],$
                   '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[3] + '/h986.cosmo50cmb.3072g1bwK.' + steps[3],$
                   '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[4] + '/h986.cosmo50cmb.3072gs1MbwK.' + steps[4],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[5],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[6] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[6],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[7] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[7]]

        pfiles =  [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
                   '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.param',$
                  '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']

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

        outext = ['h603_' + steps[0] + '_orig',$
                  'h603_' + steps[1] + '_met',$
                  'h603_' + steps[2] + '_H2',$
                  'h986_' + steps[3] + '_orig',$
                  'h986_' + steps[4] + '_met',$
                  'h986_' + steps[5] + '_H2',$
                  'h516_' + steps[6] + '_met',$
                  'h516_' + steps[7] + '_H2']

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
        maxdistance = 10
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
    11: BEGIN
       steps = ['00512','00512','00512','00512','00512','00492','00512']
        dir   =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir',$
                   '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[2],$
                   '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[3],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[4] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[5] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[6] + '.dir']

        files =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[1],$
                   '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[2] + '/h986.cosmo50cmb.3072g1bwK.' + steps[2],$
                   '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[3] + '/h986.cosmo50cmb.3072gs1MbwK.' + steps[3],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[4] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[4],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[5] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[5],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[6] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[6]]

        pfiles =  [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
                   '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.param',$
                  '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']

        filebases = ['h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.halo.1',$
                     'h986.cosmo50cmb.3072g1bwK.'   + steps[2] + '.halo.1',$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[3] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[4] + '.halo.1',$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[5] + '.halo.1', $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[6] + '.halo.1']

        filenames = ['h603.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[1],$
                     'h986.cosmo50cmb.3072g1bwK.'   + steps[2],$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[3],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[4],$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[5], $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[6]]

        outext = ['h603_' + steps[0] + '_met',$
                  'h603_'+steps[1] + '_H2',$
                  'h986_'+steps[2] + '_orig',$
                  'h986_'+steps[3] + '_met',$
                  'h986_'+steps[4] + '_H2',$
                  'h516_'+steps[5] + '_met',$
                  'h516_'+steps[6] + '_H2']

        useH2 = [0,1,0,0,1,0,1]
        halos = [1,1,1,1,1,1,1]
        halos_str = strtrim(halos,2)
;keys   = ['no Metals','Metals, no H2', 'H2']
        psym = [6,8,4,6,8,6,8]
        distunits = [50000.,50000.,50000.,50000.,50000.,25000.,25000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16,1.84793e16,1.84793e16,2.310e15,2.310e15]
        linestyles = [2,0,3,2,0,2,0]
        filternums = [13,13,13,13,13,13,13]
        cameras =    [14,14,14,16,14,14,14]
        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        maxdistance = 8
        yrange_SFH = [0,7]
        haloids = [[1],[1],[1],[1],[1],[1],[1]]
        intq = [38,33,38,38,33,32,32]
        center = [[0,0],$
                  [0,0],$
                  [0,0],$
                  [-8,72],$
                  [0,0],$
                  [0,0],$
                  [0,0]]
        vfinal =  [120.751, 115.436, 96.2218, 109.535, 101.855, 52.5520, 56.3969]
    END
   12: BEGIN
       steps = ['00512','00512','00512','00512','00512','00512','00512','00492','00512']
        dir   =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g1bwK/steps/h603.cosmo50cmb.3072g1bwK.'     + steps[0] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir',$
                   '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[3],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[4] + '.dir',$
;                  '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[4],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/steps/h516.cosmo25cmb.3072g3BWK.'     + steps[6] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[7] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[8] + '.dir']
        files =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g1bwK/steps/h603.cosmo50cmb.3072g1bwK.'     + steps[0] + '.dir/h603.cosmo50cmb.3072g1bwK.'   + steps[0],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[2],$
                   '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[3] + '/h986.cosmo50cmb.3072g1bwK.' + steps[3],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[4] + '.dir/h986.cosmo50cmb.3072gs1MbwK.' + steps[4],$
;                  '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.' + steps[4] + '/h986.cosmo50cmb.3072gs1MbwK.' + steps[4],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[5],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/steps/h516.cosmo25cmb.3072g3BWK.'     + steps[6] + '.dir/h516.cosmo25cmb.3072g3BWK.'   + steps[6],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[7] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[7],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[8] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[8]]
        pfiles =  [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g1bwK/h603.cosmo50cmb.3072g1bwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
                   '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.param',$
                   '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/h516.cosmo25cmb.3072g3bwK.param',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
        filebases = ['h603.cosmo50cmb.3072g1bwK.'   + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.halo.1',$
                     'h986.cosmo50cmb.3072g1bwK.'   + steps[3] + '.halo.1',$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[4] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.halo.1',$
                     'h516.cosmo25cmb.3072g3BWK.'   + steps[6] + '.halo.1',$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[7] + '.halo.1', $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[8] + '.halo.1']
        filenames = ['h603.cosmo50cmb.3072g1bwK.'   + steps[0],$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2],$
                     'h986.cosmo50cmb.3072g1bwK.'   + steps[3],$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[4],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[5],$
                     'h516.cosmo25cmb.3072g3BWK.'   + steps[6],$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[7], $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[8]]
        outext = ['_h603_'+steps[0] + '_warm',$
                  '_h603_'+steps[1] + '_met',$
                  '_h603_'+steps[2] + '_H2',$
                  '_h986_'+steps[3] + '_warm',$
                  '_h986_'+steps[4] + '_met',$
                  '_h986_'+steps[5] + '_H2',$
                  '_h516_'+steps[6] + '_warm',$
                  '_h516_'+steps[7] + '_met',$
                  '_h516_'+steps[8] + '_H2']
        useH2 = [0,0,1,0,0,1,0,0,1]
        halos = [1,1,1,1,1,1,1,1,1]
        halos_str = strtrim(halos,2)
        distunits = [50000.,50000.,50000.,50000.,50000.,50000.,25000.,25000.,25000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16,1.84793e16,1.84793e16,1.84793e16,2.310e15,2.310e15,2.310e15]

        filternums = [13,13,13,13,13,13,13,13,13] ;number of the filters extension in broadband.fits
;       cameras =    [14,14,14,14,16,14,14,14,14]
;       rotateAngle = [3,3,3,3,3,3,3,3,3,3]
        cameras =    [15,15,15,15,15,15,15,15,15] ;number of the CAMERA-1-BROADBAND extension in broadband.fits 
        rotateAngle = [5,5,5,5,5,5,5,5,5] ;key for IDL to say how the Sunrise image must be rotated by to match the HI cube
        recenter =   [0,0,0,0,0,0,0,0,0]
        center = [[0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0]]
        angle = [45.0,45.0,45.0,45.0,45.0,45.0,45.0,45.0,45.0] ;angle the galaxy is at
        intq = [33,33,33,38,33,33,33,32,32] ; number of the spectrum extension is mcrx.fits
        resizefactor = [2,2,2,2,2,2,1,1,1] ; resize factor for the jpeg image
        bands        = [5,5,5,5,5,5,5,5,5] ; number for bands (2MASS H)

        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        maxdistance = 8;10        
        yrange_SFH = [0,7]

        haloids = [[1],[1],[1],[1],[1],[1],[1],[1],[1]]
;        vfinal =  [120.751, 115.436, 96.2218, 109.535, 101.855, 52.5520, 56.3969]
        vfinal =  [102.159, 122.812, 111.321, 96.2218, 109.535, 101.855,57.2425,52.5520,56.3969]

;       colors = [50,50,130,130,130,254,254,254]
        colors = [120,fgcolor,254,120,fgcolor,254,120,fgcolor]
        obscolor = 50
;       linestyles = [3,2,0,3,2,0,3,2,0]
        linestyles = [0,0,0,0,0,0,0,0,0]
;       psym = [4,6,8,4,6,8,4,6,8]
        psym = [4,6,16,4,6,16,4,6,16]
        ctables = [0,0,39,0,0,39,0,0,39]
        IF KEYWORD_SET(outplot) THEN thicks = [5,5,8,5,5,8,5,5,8] ELSE thicks = [1,1,3,1,1,3,1,1,3]
;        keys   = ['no Metals, no ' + textoidl('H_2'),'Metals, no ' + textoidl('H_2'),  textoidl('H_2')]
        keys   = ['Primordial','Metals','Metals + ' + textoidl('H_2')]
    END
    13: BEGIN
       steps = ['00512','00512','00512','00512','00492','00512']
        dir   =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[2] + '.dir',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[3] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[4] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[5] + '.dir']

        files =   [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[1],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[2] + '.dir/h986.cosmo50cmb.3072gs1MbwK.' + steps[2],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[3] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[3],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[4] + '.dir/h516.cosmo25cmb.3072g1MBWK.'  + steps[4],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[5] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[5]]

        pfiles =  [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
                  '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']

        filebases = ['h603.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[1] + '.halo.1',$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[2] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[3] + '.halo.1',$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[4] + '.halo.1', $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[5] + '.halo.1']

        filenames = ['h603.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[1],$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[2],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[3],$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[4], $
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[5]]

        outext = ['h603_' + steps[0] + '_met',$
                  'h603_'+steps[1] + '_H2',$
                  'h986_'+steps[2] + '_met',$
                  'h986_'+steps[3] + '_H2',$
                  'h516_'+steps[4] + '_met',$
                  'h516_'+steps[5] + '_H2']

        useH2 = [0,1,0,1,0,1]
        halos = [1,1,1,1,1,1]
        halos_str = strtrim(halos,2)
;keys   = ['no Metals','Metals, no H2', 'H2']
        psym = [6,8,6,8,6,8]
        distunits = [50000.,50000.,50000.,50000.,25000.,25000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16,1.84793e16,2.310e15,2.310e15]
        linestyles = [2,0,2,0,2,0]
        filternums = [13,13,13,13,13,13]
        cameras =    [14,14,16,14,14,14]
        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        maxdistance = 8
        yrange_SFH = [0,7]
        haloids = [[1],[1],[1],[1],[1],[1]]
        intq = [38,33,38,33,32,32]
        center = [[0,0],$
                  [0,0],$
                  [-8,72],$
                  [0,0],$
                  [0,0],$
                  [0,0]]
        vfinal =  [120.751, 115.436, 109.535, 101.855, 52.5520, 56.3969]
        rmax = [7.38417,8.01984,8,8,2,2]
    END
    14: BEGIN
        steps = ['00512','00512']
        dir =    [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir']
       files =   [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.dir/h986.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[1] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[1]]
        pfiles = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param']
        filebases = [ 'h986.cosmo50cmb.3072gs1MbwK.' + steps[0] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[1] + '.halo.1']
        filenames = [ 'h986.cosmo50cmb.3072gs1MbwK.' + steps[0],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[1]]
        outext = ['met',$
                  'H2']
        useH2 = [0,1]
        halos = [1,1]
        halos_str = strtrim(halos,2)
        distunits = [50000.,50000.]
        massunits  = [1.84793e16,1.84793e16]

        filternums = [13,13] ;[13,14,13]
;        cameras = [14,16,14]    ;[14,15,14]
;        rotateAngle = [3,5,3]
;        angle = [90.0,90.0,90.0]
;        intq = [38,38,33]
;FO        cameras = [14,16,14]    ;[14,15,14]
;        rotateAngle = [3,3,3]       
        cameras = [15,15] ;45
        rotateAngle = [5,5] 
        recenter = [0,0]
        center = [[0,0],$
                  [0,0]]
        angle = [45.0,45.0]
        intq = [33,33]
        resizefactor = [2,2]
        bands        = [5,5]

        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8 ;16
        maxdistance = 8;10
        yrange_SFH = [0,7]
        xrange_fb = [0,10]
        yrange_fb = [0,4.0]
;        normalize = 1
;        yrange_fb = [0,0.03]
        yrange_fbcum = [0,0.15]

;orig satellites:       17, 12                
;satellities (metal): 
;at least 20,000 DM particles
        haloids = [[1, 2, 3, 4,  5,  7, 10],$
                   [1, 2, 3, 9, 11, 14, 17]]
        vfinal = [109.535, 101.855]

        colors = [80,245]
        thicks = [6,6]
        linestyles = [0,0]
        psym = [6,16]
;        keys   = ['Metals, no ' + textoidl('H_2'),  textoidl('H_2')] 
;        label = 'h986'
    END
    15: BEGIN
        steps = ['00512','00512','00512']
        dir =    [prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3BWK/steps/h799.cosmo25cmb.3072g3BWK.'     + steps[0] + '.dir', $
;prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3bwdK/steps/h799.cosmo25cmb.3072g3bwdK.' + steps[0] + '.dir',$ 
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g1MBWK/steps/h799.cosmo25cmb.3072g1MBWK.'   + steps[1] + '.dir', $
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.' + steps[2] + '.dir']
        files =  [prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3BWK/steps/h799.cosmo25cmb.3072g3BWK.'     + steps[0] + '.dir/h799.cosmo25cmb.3072g3BWK.' + steps[0],$
;prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3bwdK/steps/h799.cosmo25cmb.3072g3bwdK.' + steps[0] + '.dir/h799.cosmo25cmb.3072g3bwdK.' + steps[0],$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g1MBWK/steps/h799.cosmo25cmb.3072g1MBWK.'   + steps[1] + '.dir/h799.cosmo25cmb.3072g1MBWK.' + steps[1], $
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.' + steps[2] + '.dir/h799.cosmo25cmb.3072g14HBWK.' + steps[2]]
        pfiles = [prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3BWK/h799.cosmo25cmb.3072g3BWK.param',$
;prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3bwdK/h799.cosmo25cmb.3072g3bwdK.param',$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g1MBWK/h799.cosmo25cmb.3072g1MBWK.param', $
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/h799.cosmo25cmb.3072g14HBWK.param']
        filebases = ['h799.cosmo25cmb.3072g3BWK.' + steps[0] + '.halo.1',$
;'h799.cosmo25cmb.3072g3bwdK.' + steps[0] + '.halo.1',$
                     'h799.cosmo25cmb.3072g1MBWK.' + steps[1] + '.halo.1',$
                     'h799.cosmo25cmb.3072g14HBWK.' + steps[2] + '.halo.1']
        filenames = ['h799.cosmo25cmb.3072g3BWK.' + steps[0],$
;'h799.cosmo25cmb.3072g3bwdK.' + steps[0],$
                     'h799.cosmo25cmb.3072g1MBWK.' + steps[1],$
                     'h799.cosmo25cmb.3072g14HBWK.' + steps[2]]
        outext = ['warm',$
                  'met',$
                  'H2']
        useH2 = [0,0,1]
        halos = [1,1,1]
        halos_str = strtrim(halos,2)
        distunits = [25000.,25000.,25000.]
        massunits  = [2.310e15,2.310e15,2.310e15]

        filternums = [13,13,13] ;[13,14,13]
;        cameras = [14,16,14]    ;[14,15,14]
;        rotateAngle = [3,5,3]
;        angle = [90.0,90.0,90.0]
;        intq = [38,38,33]
;FO        cameras = [14,16,14]    ;[14,15,14]
;        rotateAngle = [3,3,3]       
        cameras = [15,15,15] ;45
        rotateAngle = [5,5,5] 
        recenter = [0,0,0]
        center = [[0,0],$
                  [0,0],$
                  [0,0]]
        angle = [45.0,45.0,45.0]
        intq = [33,33,33]
        resizefactor = [2,2,2]
        bands        = [5,5,5]

        yrange_vcirc = [0,75]
        yrange_photo = [26,19]
        maxdistance_photo = 4 ;16
        maxdistance = 8;10
        yrange_SFH = [0,0.4]
        yrange_fbcum = [0,0.15];12]
        xrange_fb = [0,4]
        yrange_fb = [0,0.12]
;        normalize = 1;1  
;        yrange_fb = [0,0.07]

;orig satellites:       17, 12                
;satellities (metal): 
;at least 20,000 DM particles

        haloids = [[1,6,11,15,16],$
                   [1,6,11,15,16],$
                   [1,6,10,15,16]]
;        vfinal = [52.7254,53.9138,50.9242]
        vfinal = [57.6719,60.0327,54.4270]
        ctables = [0,0,39]
        IF KEYWORD_SET(outplot) THEN thicks = [5,5,8] ELSE thicks = [1,1,3]
        colors = [120,fgcolor,254]
        obscolor = 50
        linestyles = [2,0,0]
        psym = [4,6,16]
;        keys   = ['no Metals, no ' + textoidl('H_2'),'Metals, no ' + textoidl('H_2'),  textoidl('H_2')] 
        keys   = ['Primordial','Metals',textoidl('H_2')]
        label = 'h799'
    END
   16: BEGIN
       steps = ['00512','00512','00512','00512','00512','00512','00492','00492','00492','00512','00512','00512']
        dir   =   [$
;prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g1bwK/steps/h603.cosmo50cmb.3072g1bwK.'     + steps[0] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$                   
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/steps/h986.cosmo50cmb.3072g1bwK.' + steps[3] + '.dir',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[4] + '.dir',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/steps/h516.cosmo25cmb.3072g3BWK.'     + steps[6] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[7] + '.dir',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[8] + '.dir',$
                   prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3BWK/steps/h799.cosmo25cmb.3072g3BWK.'     + steps[9] + '.dir',$
;                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3bwdK/steps/h799.cosmo25cmb.3072g3bwdK.'   + steps[9] + '.dir',$ 
                   prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g1MBWK/steps/h799.cosmo25cmb.3072g1MBWK.'   + steps[10] + '.dir',$
                   prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.' + steps[11] + '.dir']
        files =   [$
;prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g1bwK/steps/h603.cosmo50cmb.3072g1bwK.'     + steps[0] + '.dir/h603.cosmo50cmb.3072g1bwK.'    + steps[0],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.'  + steps[1],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.'  + steps[1],$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h603.cosmo50cmb.3072g14HBWK.'  + steps[2],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/steps/h986.cosmo50cmb.3072g1bwK.' + steps[3] + '.dir/h986.cosmo50cmb.3072g1bwK.' + steps[3],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.' + steps[4] + '.dir/h986.cosmo50cmb.3072gs1MbwK.'  + steps[4],$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.dir/h986.cosmo50cmb.3072g14HBWK.'  + steps[5],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/steps/h516.cosmo25cmb.3072g3BWK.'     + steps[6] + '.dir/h516.cosmo25cmb.3072g3BWK.'    + steps[6],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.'   + steps[7] + '.dir/h516.cosmo25cmb.3072g1MBWK.'   + steps[7],$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[8] + '.dir/h516.cosmo25cmb.3072g14HBWK.'  + steps[8],$
                   prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3BWK/steps/h799.cosmo25cmb.3072g3BWK.'     + steps[9] + '.dir/h799.cosmo25cmb.3072g3BWK.'    + steps[9],$ 
;                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3bwdK/steps/h799.cosmo25cmb.3072g3bwdK.'   + steps[9] + '.dir/h799.cosmo25cmb.3072g3bwdK.'  + steps[9],$
                   prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g1MBWK/steps/h799.cosmo25cmb.3072g1MBWK.'   + steps[10] + '.dir/h799.cosmo25cmb.3072g1MBWK.'  + steps[10], $
                   prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.' + steps[11] + '.dir/h799.cosmo25cmb.3072g14HBWK.' + steps[11]]
        pfiles =  [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g1bwK/h603.cosmo50cmb.3072g1bwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.param',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param',$
                   prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/h516.cosmo25cmb.3072g3bwK.param',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
                   prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param',$
                   prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3BWK/h799.cosmo25cmb.3072g3BWK.param',$
;                   prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3bwdK/h799.cosmo25cmb.3072g3bwdK.param',$
                   prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g1MBWK/h799.cosmo25cmb.3072g1MBWK.param', $
                   prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/h799.cosmo25cmb.3072g14HBWK.param']
        filebases = [$
;'h603.cosmo50cmb.3072g1bwK.'   + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2] + '.halo.1',$
                     'h986.cosmo50cmb.3072g1bwK.'   + steps[3] + '.halo.1',$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[4] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[5] + '.halo.1',$
                     'h516.cosmo25cmb.3072g3BWK.'   + steps[6] + '.halo.1',$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[7] + '.halo.1',$
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[8] + '.halo.1',$
                     'h799.cosmo25cmb.3072g3BWK.'   + steps[9] + '.halo.1',$
;                     'h799.cosmo25cmb.3072g3bwdK.'  + steps[9] + '.halo.1',$
                     'h799.cosmo25cmb.3072g1MBWK.'  + steps[10] + '.halo.1',$
                     'h799.cosmo25cmb.3072g14HBWK.' + steps[11] + '.halo.1']
        filenames = [$
;'h603.cosmo50cmb.3072g1bwK.'   + steps[0],$
                    'h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h603.cosmo50cmb.3072g14HBWK.' + steps[2],$
                     'h986.cosmo50cmb.3072g1bwK.'   + steps[3],$
                     'h986.cosmo50cmb.3072gs1MbwK.' + steps[4],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[5],$
                     'h516.cosmo25cmb.3072g3BWK.'   + steps[6],$
                     'h516.cosmo25cmb.3072g1MBWK.'  + steps[7],$
                     'h516.cosmo25cmb.3072g14HBWK.' + steps[8],$
                     'h799.cosmo25cmb.3072g3BWK.'   + steps[9],$
;                     'h799.cosmo25cmb.3072g3bwdK.'  + steps[9],$
                     'h799.cosmo25cmb.3072g1MBWK.'  + steps[10],$
                     'h799.cosmo25cmb.3072g14HBWK.' + steps[11]]
        outext = ['_h603_'+steps[0] + '_warm',$
                  '_h603_'+steps[1] + '_met',$
                  '_h603_'+steps[2] + '_H2',$
                  '_h986_'+steps[3] + '_warm',$
                  '_h986_'+steps[4] + '_met',$
                  '_h986_'+steps[5] + '_H2',$
                  '_h516_'+steps[6] + '_warm',$
                  '_h516_'+steps[7] + '_met',$
                  '_h516_'+steps[8] + '_H2',$
                  '_h799_'+steps[9] + '_warm',$
                  '_h799_'+steps[9] + '_met',$
                  '_h799_'+steps[10] + '_H2']
        useH2 = [0,0,1,0,0,1,0,0,1,0,0,1]
        halos = [1,1,1,1,1,1,1,1,1,1,1,1]
        halos_str = strtrim(halos,2)
        distunits = [50000.,50000.,50000.,50000.,50000.,50000.,25000.,25000.,25000.,25000.,25000.,25000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16,1.84793e16,1.84793e16,1.84793e16,2.310e15,2.310e15,2.310e15,2.310e15,2.310e15,2.310e15]

        filternums = [13,13,13,13,13,13,13,13,13,13,13,13] ;number of the filters extension in broadband.fits
        cameras =    [15,15,15,15,15,15,15,15,15,15,15,15] ;number of the CAMERA-1-BROADBAND extension in broadband.fits 
        rotateAngle = [5,5,5,5,5,5,5,5,5,5,5,5] ;key for IDL to say how the Sunrise image must be rotated by to match the HI cube
        recenter =   [0,0,0,0,0,0,0,0,0,0,0,0]
        center = [[0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0],$
                  [0,0]]
        angle = [45.0,45.0,45.0,45.0,45.0,45.0,45.0,45.0,45.0,45.0,45.0,45.0] ;angle the galaxy is at
        intq = [33,33,33,38,33,33,33,32,33,33,33,33] ; number of the spectrum extension is mcrx.fits
        resizefactor = [2,2,2,2,2,2,1,1,1,1,1,1] ; resize factor for the jpeg image
        bands        = [5,5,5,5,5,5,5,5,5,5,5,5] ; number for bands (2MASS H)

        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8
        maxdistance = 8;10        
        yrange_SFH = [0,7]

        haloids = [[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]]
 ;       vfinal =  [102.159, 122.812, 111.321, 96.2218, 109.535, 101.855, 57.2425, 52.5520, 56.3969, 52.7254, 53.9138, 50.9242]
;        vfinal =  [104.054, 140.760, 110.402, 98.2964, 138.677, 98.8127, 60.5869, 59.5379, 65.1188, 57.6719, 60.0327, 54.4270]
        vfinal =  [103.925, 131.990, 111.032, 97.2769, 130.544, 100.847, 60.7539, 59.8370, 65.1547, 57.8003, 60.3383, 54.4358]
;        vfinal = [57.6719,60.0327,54.4270]
        colors = [120,fgcolor,254,120,fgcolor,254,120,fgcolor,254,120,fgcolor,254]
        obscolor = 50
        linestyles = [2,0,0,2,0,0,2,0,0,2,0,0]
        psym = [4,6,16,4,6,16,4,6,16,4,6,16]
        ctables = [0,0,39,0,0,39,0,0,39,0,0,39]
        IF KEYWORD_SET(outplot) THEN thicks = [5,5,8,5,5,8,5,5,8,5,5,8] ELSE thicks = [1,1,3,1,1,3,1,1,3,1,1,3]
;        keys   = ['no Metals, no ' + textoidl('H_2'),'Metals, no ' + textoidl('H_2'),  textoidl('H_2')]
        keys   = ['Primordial','Metals',textoidl('H_2')]
    END
    17: BEGIN
        steps = ['00512','00512','00512']
        dir =    ['/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[0],$
                  prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir']
       files =   ['/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.' + steps[0] + '/h986.cosmo50cmb.3072g1bwK.' + steps[0],$
                  prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.dir/h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.' + steps[2] + '.dir/h986.cosmo50cmb.3072g14HBWK.' + steps[2]]
        pfiles = ['/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.param',$
                  prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param']
        filebases = ['h986.cosmo50cmb.3072g1bwK.'   + steps[0] + '.halo.1',$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1] + '.halo.1',$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[2] + '.halo.1']
        filenames = ['h986.cosmo50cmb.3072g1bwK.'   + steps[0],$
                     'h603.cosmo50cmb.3072gs1MbwK.' + steps[1],$
                     'h986.cosmo50cmb.3072g14HBWK.' + steps[2]]
        outext = ['warm',$
                  'met',$
                  'H2']
        useH2 = [0,0,1]
        halos = [1,1,1]
        halos_str = strtrim(halos,2)
        distunits = [50000.,50000.,50000.]
        massunits  = [1.84793e16,1.84793e16,1.84793e16]

        filternums = [13,13,13] ;[13,14,13]
;        cameras = [14,16,14]    ;[14,15,14]
;        rotateAngle = [3,5,3]
;        angle = [90.0,90.0,90.0]
;        intq = [38,38,33]
;FO        cameras = [14,16,14]    ;[14,15,14]
;        rotateAngle = [3,3,3]       
        cameras = [15,15,15] ;45
        rotateAngle = [5,5,5] 
        recenter = [0,0,0]
        center = [[0,0],$
                  [0,0],$
                  [0,0]]
        angle = [45.0,45.0,45.0]
        intq = [38,33,33]
        resizefactor = [2,2,2]
        bands        = [5,5,5]

        yrange_vcirc = [0,250]
        yrange_photo = [26,15]
        maxdistance_photo = 8 ;16
        maxdistance = 8;10
        yrange_SFH = [0,7]
        xrange_fb = [0,8]
        yrange_fb = [0,1.2]
;        normalize = 1
;        yrange_fb = [0,0.03]
        yrange_fbcum = [0,0.3]

;orig satellites:       17, 12                
;satellities (metal): 
;at least 20,000 DM particles
        haloids = [[1, 2, 4, 9, 11, 13, 14],$
                   [1, 2, 3, 4,  5,  7, 10],$
                   [1, 2, 3, 9, 11, 14, 17]]
        vfinal = [96.2218, 109.535, 101.855]

        colors = [120,fgcolor,254]      
;        colors = [50,120,245]
;        thicks = [6,6,6]
        IF KEYWORD_SET(outplot) THEN thicks = [5,5,8] ELSE thicks = [1,1,3]
        linestyles = [2,0,0]
        psym = [4,6,16]
        ctables = [0,0,39]
;        ctables = [39,39,39]
        keys   = ['no Metals, no ' + textoidl('H_2'),'Metals, no ' + textoidl('H_2'),  textoidl('H_2')] 
       label = 'h986'
    END
ENDCASE
n = N_ELEMENTS(files)

dotpos = strpos(filenames,'.',/reverse_search)
fileprefixes = filenames
FOR i = 0, n - 1 DO fileprefixes[i] = strmid(filenames[i],0,dotpos[i])

dirhead = dir
FOR i = 0, n - 1 DO BEGIN
    cosmopos = strsplit(dir[i],'cosmo',/regex)
    dirheadend = strpos(dir[i],'/',cosmopos[n_elements(cosmopos) -2])
    dirhead[i] = strmid(dir[i],0,dirheadend + 1)
ENDFOR

a = fltarr(n)
FOR i = 0, n - 1 DO BEGIN
    rtipsy,files[i] + '.halo.1',h,g,d,s,/justhead
    a[i] = h.time
ENDFOR

IF 0 AND KEYWORD_SET(color) THEN BEGIN
    colors = [120,fgcolor,254]
    linestyles = [2,0,0]
    ctables = [0,0,39]
    obscolor = 100

;    colors = [fgcolor,254]
;    linestyles = [0,0]
;    ctables = [0,39]
;    !p.charsize=2.0
;    !x.charsize=2.0             ;2.25
;    !y.charsize=2.0             ;2.25
    IF KEYWORD_SET(outplot) THEN thicks = [5,5,8] ELSE  thicks = [1,1,3]
ENDIF

IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if NOT keyword_set(ctables) then ctables = [39,39,39]
    IF NOT keyword_set(obscolor) THEN obscolor = fgcolor
    if NOT keyword_set(colors) then  colors  = (findgen(n) + 1)*254/n else colors = colors
    IF NOT KEYWORD_SET(psym) THEN psym = fltarr(n) + 4
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT KEYWORD_SET(symsizes) THEN symsizes = fltarr(n) + 1.5
ENDIF ELSE BEGIN
    loadct,0    
    if NOT keyword_set(ctables) then ctables = [0,0,0]
    IF NOT keyword_set(obscolor) THEN obscolor = 100
    If NOT keyword_set(colors) then  colors = (fltarr(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(broadband)) + 5
    IF NOT KEYWORD_SET(psym) THEN  psym = (findgen(n)+2)*2
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2 ;thicks = (findgen(n) + 1)*6/n - 1
    IF NOT KEYWORD_SET(linestyles) THEN linestyles = REVERSE(findgen(n)*2) 
    IF NOT KEYWORD_SET(symsizes) THEN symsizes = fltarr(n) + 1.5
ENDELSE

IF keyword_set(outplot) THEN BEGIN
    xsize = 18                  ;10*n
    ysize = 12                  ;12
;    mxTitSize = 1.5
;    mxTitOffset = 2
ENDIF ELSE BEGIN
    xsize = 400                 ;*n
    ysize = 350                 ;475
;    mxTitSize = 1.5
;    mxTitOffset = 1
ENDELSE

nk = N_ELEMENTS(keys)
lb_keys = strarr(nk + 2)
lb_thicks = intarr(nk + 2)
lb_color = intarr(nk + 2) + bgcolor
lb_linestyle = intarr(nk + 2) + 1
lb_symbols = intarr(nk + 2) + 3
lb_symsizes = intarr(nk + 2) + 1

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
    matches0 = matchHalos([files[0] + '.amiga.stat',files[1] + '.amiga.stat'],minmass = 10^(8.5))
;    matches0 = matchHalos([files[6] + '.amiga.stat',files[5] + '.amiga.stat'],minmass = 10^(8.5))
    IF n gt 2 THEN matches1 = matchHalos([files[1] + '.amiga.stat',files[2] + '.amiga.stat'],minmass = 10^(8.5))
ENDIF

;--------------------------- Twins analysis -------------------------
IF plot_twinsanalysis THEN BEGIN
    twinsanalysis,outplot = outplot
ENDIF

;multiHaloGuo ----------------------------------------------------------------------------------------------------------------------------------
;/astro/users/christensen/code/IDL/HIcubes/multiHaloGuo
IF plot_multiHaloGuo THEN BEGIN
    IF KEYWORD_SEt(outplot) THEN multiHaloGuo,files,key = keys,yrange = [4,10.5],haloids = haloids,halomasses = halomasses,psym = psym,xrange = [8.5,12],color = colors, obscolor = obscolor, ctables = ctables, outplot = outplot + '_MstarMhalo.eps' $
                            ELSE multiHaloGuo,files,key = keys,yrange = [4,10.5],haloids = haloids,halomasses = halomasses,psym = psym,xrange = [8.5,12],color = colors, obscolor = obscolor, ctables = ctables

;   IF KEYWORD_SET(outplot) THEN device, filename= outplot + '_MstarMhaloRatio.eps',bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch ELSE stop
;   plot,alog10(halomasses[*,0]),alog10(halomasses[*,1]),psym = 5,xtitle = 'Log(Mass Halo)',ytitle = 'Log(Mstar_H2/Mstar_{no H2})',yrange = [-1,1],xrange = [9,11.5]
;   oplot,alog10(halomasses[*,0]),alog10(halomasses[*,1])
;   IF KEYWORD_SET(outplot) THEN device,/close
    IF KEYWORD_SEt(outplot) THEN multiHaloGuo,files,key = keys,yrange = [4,10.5],haloids = haloids,halomasses = halomasses,psym = psym,xrange = [8.5,12],/kcorrect,color = colors, obscolor = obscolor, ctables = ctables,outplot = outplot + '_MstarMhalo_Obs.eps' $
                            ELSE multiHaloGuo,files,key = keys,yrange = [4,10.5],haloids = haloids,halomasses = halomasses,psym = psym,xrange = [8.5,12],/kcorrect,color = colors, obscolor = obscolor, ctables = ctables
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
IF plot_vcirc THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN
        IF keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,verbose = verbose,outfile = outplot + '_vcirc.eps',color = colors $
                              ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,verbose = verbose,outfile = outplot + '_vcirc.eps'
;        IF keyword_set(color) THEN vcirc,files[0] + '.halo.' + halos_str[0],massunits[0],distunits[0],thicks = thicks[0],linestyle = linestyles[0],yrange = yrange_vcirc,maxdistance = maxdistance,vfinal = vfinal,type = 'all',ctables =ctables[0],verbose = verbose,outfile = outplot + '_vcirc0.eps',color = colors[0]
;        IF keyword_set(color) THEN vcirc,files[1] + '.halo.' + halos_str[1],massunits[1],distunits[1],thicks = thicks[1],linestyle = linestyles[1],yrange = yrange_vcirc,maxdistance = maxdistance,vfinal = vfinal,type = 'all',ctables =ctables[1],verbose = verbose,outfile = outplot + '_vcirc1.eps',color = colors[1]
;        IF keyword_set(color) THEN vcirc,files[2] + '.halo.' + halos_str[2],massunits[2],distunits[2],thicks = thicks[2],linestyle = linestyles[2],yrange = yrange_vcirc,maxdistance = maxdistance,vfinal = vfinal,type = 'all',ctables =ctables[2],verbose = verbose,outfile = outplot + '_vcirc2.eps',color = colors[2]
    ENDIF ELSE BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,verbose = verbose,color = colors   $
                              ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,verbose = verbose
    ENDELSE
    print,vfinal
ENDIF


IF 0 THEN BEGIN
IF plot_vcirc THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'dark',outfile = outplot + '_vcircdark.eps',color = colors $
                              ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'dark',outfile = outplot + '_vcircdark.eps'
    ENDIF ELSE BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'dark',color = colors $
                              ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'dark'
    ENDELSE
ENDIF


IF plot_vcirc THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'star',outfile = outplot + '_vcircstar.eps',color = colors $
                              ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'star',outfile = outplot + '_vcircstar.eps'
    ENDIF ELSE BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,color = colors,type = 'star' $
                              ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'star'
    ENDELSE
ENDIF

IF plot_vcirc THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'gas',outfile = outplot + '_vcircgas.eps',color = colors $
                              ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'gas',outfile = outplot + '_vcircgas.eps'
    ENDIF ELSE BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'gas',color = colors $
                              ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,maxdistance = 15,vfinal = vfinal,type = 'gas'
    ENDELSE
ENDIF
ENDIF

;--------------------------- DM profile -------------
IF plot_DMprof THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN
        if keyword_set(color) THEN DMprof,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,outfile = outplot + '_DMprof.eps',color = colors,formatthick = formatthick $
                              ELSE DMprof,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,outfile = outplot + '_DMprof.eps',formatthick = formatthick
    ENDIF ELSE BEGIN
        if keyword_set(color) THEN DMprof,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,color = colors,formatthick = formatthick $
                              ELSE DMprof,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,maxdistance = maxdistance,vfinal = vfinal,type = 'all',label = label,ctables =ctables,formatthick = formatthick
    ENDELSE
ENDIF

;--------------------------- Sigma v Radius -----
IF plot_sigmar Then BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN 
    IF KEYWORD_SET(color) THEN $
      sigmar,files,pfiles,keys = keys,thicks = thicks,linestyle = linestyles,outplot = outplot + 'sigmar.eps', colors = colors, label = label,ctables =ctables ,intq = intq, camera = cameras,center = center ELSE $
      sigmar,files,pfiles,keys = keys,thicks = thicks,linestyle = linestyles,outplot = outplot + '_sigmar.eps',label = label,ctables =ctables,intq = intq, camera = cameras,center = center
      ENDIF ELSE BEGIN 
        IF KEYWORD_SET(color) THEN $
          sigmar,files,pfiles,keys = keys,thicks = thicks,linestyle = linestyles,label = label,ctables =ctables,intq = intq,camera = cameras,center = center,colors = colors ELSE $
          sigmar,files,pfiles,keys = keys,thicks = thicks,linestyle = linestyles,label = label,ctables =ctables,intq = intq,camera = cameras,center = center
    ENDELSE


ENDIF

;--------------------------- Photometric Profile ------------
IF plot_photometricProf THEN BEGIN
    camerasP = cameras*0 + 14
    IF KEYWORD_SET(outplot) THEN BEGIN 
        IF KEYWORD_SET(color) THEN photometricProf,files + '.' + halos_str + '/broadband.fits',keys = keys,thicks = thicks,linestyle = linestyles,filternums = filternums,cameras = camerasP,yrange = yrange_photo,bands = bands,maxdistance = maxdistance_photo,recenter = recenter,label = label,ctables =ctables,outplot = outplot + '_photoProf.eps',colors = colors,formatthick = formatthick $
                              ELSE photometricProf,files + '.' + halos_str + '/broadband.fits',keys = keys,thicks = thicks,linestyle = linestyles,filternums = filternums,cameras = camerasP,yrange = yrange_photo,bands = bands,maxdistance = maxdistance_photo,recenter = recenter,label = label,ctables =ctables,outplot = outplot + '_photoProf.eps',formatthick = formatthick
    ENDIF ELSE BEGIN 
        IF KEYWORD_SET(color) THEN photometricProf,files + '.' + halos_str + '/broadband.fits',keys = keys,thicks = thicks,linestyle = linestyles,filternums = filternums,cameras = camerasP,yrange = yrange_photo,bands = bands,maxdistance = maxdistance_photo,recenter = recenter,label = label,ctables =ctables,colors = colors,formatthick = formatthick $
                              ELSE photometricProf,files + '.' + halos_str + '/broadband.fits',keys = keys,thicks = thicks,linestyle = linestyles,filternums = filternums,cameras = camerasP,yrange = yrange_photo,bands = bands,maxdistance = maxdistance_photo,recenter = recenter,label = label,ctables =ctables,formatthick = formatthick
    ENDELSE
ENDIF

;-------------------------------- TF ----------------------
IF plot_tully_fisher_obs THEN BEGIN
;    vmax = fltarr(n)
;    vmaxobs = fltarr(n)
;    head = ' '    
;    FOR i = 0, n - 1 DO BEGIN
;        data = read_stat_struc_amiga(files[i] + '.amiga.stat')
;        ind = where(data.group eq halos[i])
;        vmax[i] = data[ind].vc
;        openr,1,dir[i] + '/linewidths.txt'
;        readf,1,head,format = '(A)'
;        readf,1,width,width2,vmaxline
;        close,1
;        vmaxobs[i] = vmaxline ;width2
;    ENDFOR
    gmass = tully_fisher_obs_gasmass(files + '.halo.' + halos_str, massunits,smass = smass)
 ;   IF KEYWORD_SET(outplot) THEN tully_fisher_obs,files + '.' + halos_str + '/broadband.fits',vmax,gmass,filternums = filternums,symbols = psym, key = keys, outfile = outplot + '_vtrue_' $
 ;   ELSE BEGIN
 ;       tully_fisher_obs,files + '.' + halos_str + '/broadband.fits',vmax,gmass,filternums = filternums,symbols = psym, key = keys
 ;       stop
 ;   ENDELSE

;;    IF KEYWORD_SET(outplot) THEN tully_fisher_obs,files + '.' + halos_str + '/broadband.fits',vmaxobs,gmass,filternums = filternums,symbols = psym, key = keys, outfile = outplot + '_vobs_' $
;;    ELSE BEGIN
;;        tully_fisher_obs,files + '.' + halos_str + '/broadband.fits',vmaxobs,gmass,filternums = filternums,symbols = psym, key = keys
;;        stop
;;    ENDELSE

;    IF KEYWORD_SET(outplot) THEN tully_fisher_obs_btf, files + '.' + halos_str + '/broadband.fits', vmaxobs, gmass, key = keys,smass_true = smass,velocities_true = vfinal,filternums = filternums,symbols = psym, outfile = outplot + '_vAll_' $
;    ELSE BEGIN
;        tully_fisher_obs_btf, files + '.' + halos_str + '/broadband.fits', vmaxobs, gmass, key = keys,smass_true = smass,velocities_true = vfinal,filternums = filternums,symbols = psym
;        stop
;    ENDELSE
    IF KEYWORD_SET(outplot) THEN BEGIN
        IF KEYWORD_SET(color) THEN tully_fisher_obs_btf, files + '.' + halos_str + '/broadband.fits', vfinal, gmass, key = keys,filternums = filternums,symbols = psym, outfile = outplot + '', symsizes = symsizes,obscolor = obscolor,ctables = ctables,color = colors $
                              ELSE tully_fisher_obs_btf, files + '.' + halos_str + '/broadband.fits', vfinal, gmass, key = keys,filternums = filternums,symbols = psym, outfile = outplot + '', symsizes = symsizes,obscolor = obscolor
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(color) THEN tully_fisher_obs_btf, files + '.' + halos_str + '/broadband.fits', vfinal, gmass, key = keys,filternums = filternums,symbols = psym, symsizes = symsizes,obscolor = obscolor,ctables = ctables, color = colors $
                              ELSE tully_fisher_obs_btf, files + '.' + halos_str + '/broadband.fits', vfinal, gmass, key = keys,filternums = filternums,symbols = psym, symsizes = symsizes,obscolor = obscolor
    ENDELSE
ENDIF


;----------------------------------- Density Radius ----------
IF plot_densRadius THEN BEGIN
    if keyword_set(outplot) THEN BEGIN
;        IF KEYWORD_SET(color) THEN $
;          densRadius,files,distunits,massunits,maxdistance = maxdistance_dens,label = label,outplot = outplot,color = colors,ctables = ctables ELSE $
          densRadius,files,distunits,massunits,maxdistance = maxdistance_dens,label = label,outplot = outplot
    ENDIF ELSE BEGIN
;        IF KEYWORD_SET(color) THEN $
;          densRadius,files,distunits,massunits,maxdistance = maxdistance_dens,label = label,color = colors,ctables = ctables ELSE $
          densRadius,files,distunits,massunits,maxdistance = maxdistance_dens,label = label
    ENDELSE
ENDIF

;------------------------------------ SFH ----------------------
IF plot_sfh THEN BEGIN
;    IF KEYWORD_SET(color) THEN BEGIN
;        IF KEYWORD_SET(outplot) THEN sfh_z,dir,files  + '.halo.' + halos_str,pfiles,keys = keys,linestyle = linestyles,thick = thicks,colors = colors,yrange = yrange_SFH,label = label,ctables =ctables,binsize = 1e8,xrange = [0,14],outplot = outplot + '_SFH.eps',formatthick = formatthick $
;                                ELSE sfh_z,dir,files  + '.halo.' + halos_str,pfiles,keys = keys,linestyle = linestyles,thick = thicks,colors = colors,yrange = yrange_SFH,label = label,ctables =ctables,binsize = 1e8,xrange = [0,14],formatthick = formatthick
;    ENDIF ELSE BEGIN
;        IF KEYWORD_SET(outplot) THEN sfh_z,dir,files  + '.halo.' + halos_str,pfiles,keys = keys,linestyle = linestyles,thick = thicks,colors = colors,yrange = yrange_SFH,label = label,ctables =ctables,binsize = 1e8,xrange = [0,14],outplot = outplot + '_SFH.eps',formatthick = formatthick $
;                                ELSE sfh_z,dir,files  + '.halo.' + halos_str,pfiles,keys = keys,linestyle = linestyles,thick = thicks,colors = colors,yrange = yrange_SFH,label = label,ctables =ctables,binsize = 1e8,xrange = [0,14],formatthick = formatthick
;     ENDELSE
    IF KEYWORD_SET(color) THEN BEGIN
        IF KEYWORD_SET(outplot) THEN sfh,files  + '.halo.' + halos_str,pfiles,key = keys,linestyle = linestyles,thick = thicks,colors = colors,yrange = yrange_SFH,label = label,ctables =ctables,binsize = 1e8,outplot = outplot,formatthick = formatthick,/bulge $
                                ELSE sfh,files  + '.halo.' + halos_str,pfiles,key = keys,linestyle = linestyles,thick = thicks,colors = colors,yrange = yrange_SFH,label = label,ctables =ctables,binsize = 1e8,formatthick = formatthick, /bulge
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(outplot) THEN sfh,files  + '.halo.' + halos_str,pfiles,key = keys,linestyle = linestyles,thick = thicks,colors = colors,yrange = yrange_SFH,label = label,ctables =ctables,binsize = 1e8,outplot = outplot,formatthick = formatthick, /bulge $
                                ELSE sfh,files  + '.halo.' + halos_str,pfiles,key = keys,linestyle = linestyles,thick = thicks,colors = colors,yrange = yrange_SFH,label = label,ctables =ctables,binsize = 1e8,formatthick = formatthick, /bulge
    ENDELSE
ENDIF

;--------------------------Global SK law ----------------------------
IF  plot_schmidtlaw_global_obs THEN BEGIN
    data=fltarr(2,N_ELEMENTS(dir))
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
        print,dir[i]
        data[*,i] = schmidtlaw_global_obs(dir[i],filenames[i],pfiles[i],useH2 = useH2[i],halo_str = halos_str[i],angle = angle[i],intq = intq[i], camera = cameras[i],center = center[*,i],verbose = verbose,/intrinsic);,/Halpha,tipsyfile = filebases[i]
    ENDFOR 

    if (KEYWORD_SET(outplot)) then begin
        device,filename = outplot + '_SK.eps',/encapsulated,/color,/times,ysize=xsize,xsize=xsize,bits_per_pixel= 8
    endif else begin
        set_plot,'x'
        window,1
    endelse
    IF KEYWORD_SET(color) THEN loadct,39
    xsigma = 10.0^(findgen(600)/100 - 1.) ;xsigmalow
    ysigma=2.5e-4*xsigma^1.4
    
    readcol,datafile_prefix + 'HIcubes/ks98.dat',name,D,logHI,logH2,logH,logSFR,tdyn,co,HI,halpha ;,format='(A9D)'
    readcol,datafile_prefix + 'HIcubes/uavpl.dat',logHdarf,logSFRdwarf
    loadct,ctables[0]
    plot,alog10(xsigma),alog10(ysigma),ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-0.5,2.5], yrange = [-4,-0.5], xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]"   
    if keyword_set(color) THEN loadct,39 else loadct,0
    oplot,logH,logSFR,psym = 2,color = obscolor,symsize = 2
    oplot,logHdarf,logSFRdwarf,psym = 1,color = obscolor,symsize = 2 
    IF KEYWORD_SET(keys) THEN BEGIN
        l_keys = lb_keys
        l_keys[0] = 'Kennicutt et al. 1998'
        l_color = lb_color
        l_color[0] = obscolor
        l_symbols = lb_symbols
        l_symbols[0]  = 2
        l_symsizes = lb_symsizes
        l_symsizes[0]  = 2   
        l_linestyle = lb_linestyle
        legend,l_keys,color = l_color,linestyle = l_linestyle,symsize = l_symsizes,psym = l_symbols,/right,/bottom,box = 0
    ENDIF
    IF KEYWORD_SET(keys) THEN BEGIN
        l_keys = lb_keys
        l_keys[1] = 'Roychowdhury et al. 2009'
        l_color = lb_color
        l_color[1] = obscolor
        l_symbols = lb_symbols
        l_symbols[1]  = 1
        l_symsizes = lb_symsizes
        l_symsizes[1]  = 2  
        l_linestyle = lb_linestyle
        legend,l_keys,linestyle = l_linestyle,color = l_color,symsize = l_symsizes,psym = l_symbols,/right,/bottom,box = 0
    ENDIF
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
        loadct,ctables[i]
        oplot,[alog10(data[0,i]),alog10(data[0,i])],[alog10(data[1,i]),alog10(data[1,i])],psym = symcat(psym[i]),color = colors[i],symsize = symsizes[i] 
        IF KEYWORD_SET(keys) AND i lt N_ELEMENTS(keys) THEN BEGIN
            l_keys = lb_keys
            l_keys[i + 2] = keys[i]
            l_color = lb_color
            l_color[i + 2] = colors[i]
            l_symbols = lb_symbols
            l_symbols[i + 2]  = psym[i]
            l_symsizes = lb_symsizes
            l_symsizes[i + 2]  = symsizes[i]            
            legend,l_keys,color = l_color,symsize = l_symsizes,psym = l_symbols,/right,/bottom,box = 0
        ENDIF
    ENDFOR
    oplot,[-0.2,0.2],[-1,-1]
    oplot,[0,0],[-0.8,-1.2]
;    IF KEYWORD_SET(keys) THEN legend,[keys,'Kennicutt 98'],psym = [psym,2],color = [colors,obscolor],/bottom,/right
    if (KEYWORD_SET(outplot)) then device,/close else stop
ENDIF

;-------------------------- Res SK law ----------------------------
IF plot_schmidtlaw_res_obs THEN BEGIN
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
        cd,dir[i]
        center_temp = [-1.0*center[0,i],center[1,i]] 
        schmidtlaw_res_obs,filenames[i],pfiles[i],useH2 = useH2[i],angle = angle[i],extno = cameras[i],center = center_temp,rotateAngle = rotateAngle[i],verbose = verbose
    ENDFOR
ENDIF

;------------- Resolved K-S law plot ----------------------
;~/code/HIcubes/schmidtlaw_res_obs.pro
res = '0.750000'
IF plot_schmidtlaw_res_obs_master_out THEN BEGIN
        sk_files_P = [prefix +  'dwarfs.g1bwk.schmidtlaw_res_obs0.750000.dat',$
                          prefix + 'spirals.g1bwk.schmidtlaw_res_obs0.750000.dat']
        IF KEYWORD_SET(outplot) THEN $
          schmidtlaw_res_obs_master_out,  sk_files_P,color = [50,50],  /contour,thick = [thicks[0],thicks[0]],label = keys[0],symbols = [14,14],outplot = outplot+'_bin_P_' ELSE BEGIN
            schmidtlaw_res_obs_master_out,sk_files_P,color = [50,50],  /contour,thick = [thicks[0],thicks[0]],label = keys[0],symbols = [14,14]
            stop
        ENDELSE

        sk_files_M = [prefix +  'dwarfs.g1MbwK.schmidtlaw_res_obs0.750000.dat',$
                          prefix + 'spirals.g1MbwK.schmidtlaw_res_obs0.750000.dat']
        IF KEYWORD_SET(outplot) THEN $
          schmidtlaw_res_obs_master_out,  sk_files_M,color = [130,130],/contour,thick = [thicks[0],thicks[0]],label = keys[1],symbols = [15,15],outplot = outplot+'_bin_M_' ELSE BEGIN
            schmidtlaw_res_obs_master_out,sk_files_M,color = [130,130],/contour,thick = [thicks[0],thicks[0]],label = keys[1],symbols = [15,15]
            stop
        ENDELSE

        sk_files_H = [prefix +  'dwarfs.g14HBWK.schmidtlaw_res_obs0.750000.dat',$
                          prefix + 'spirals.g14HBWK.schmidtlaw_res_obs0.750000.dat']
        IF KEYWORD_SET(outplot) THEN $
          schmidtlaw_res_obs_master_out,  sk_files_H,color = [254,254],/contour,thick = [thicks[0],thicks[0]],label = keys[2],symbols = [16,16],outplot = outplot+'_bin_H_' ELSE BEGIN
            schmidtlaw_res_obs_master_out,sk_files_H,color = [254,254],/contour,thick = [thicks[0],thicks[0]],label = keys[2],symbols = [16,16]
            stop
        ENDELSE


;    IF KEYWORD_SET(color) THEN $
;      schmidtlaw_res_obs_master_out,dir+"/schmidtlaw_res_obs"+res+".dat",outplot = outplot,thick = thicks,symbols = 1+fltarr(n),key = keys,symsize = symsize,color = (findgen(n)+1)*254.0/n ELSE $
;      schmidtlaw_res_obs_master_out,dir+"/schmidtlaw_res_obs"+res+".dat",outplot = outplot,thick = thicks,symbols = 1+fltarr(n),key = keys,symsize = symsize

    IF 0 THEN BEGIN
        loadct,39
        sk_files = [prefix + 'hAll.g1bwk.schmidtlaw_res_obs0.750000.dat',$
                    prefix + 'hAll.g1MbwK.schmidtlaw_res_obs0.750000.dat',$
                    prefix + 'hAll.g14HBWK.schmidtlaw_res_obs0.750000.dat']
        IF KEYWORD_SET(outplot) THEN $
          schmidtlaw_res_obs_master_out,sk_files,color = [50,130,254],/contour,thick = thicks,/multiframe,outplot = outplot+'_bin_',thick = 4+fltarr(n) ELSE $
          schmidtlaw_res_obs_master_out,sk_files,color = [50,130,254],/contour,thick = thicks,/multiframe
    ENDIF

    IF 0 THEN BEGIN
        IF KEYWORD_SET(color) THEN $
          schmidtlaw_res_obs_master_out,dir+"/schmidtlaw_res_obs"+res+".dat",outplot = outplot,thick = thicks,symbols = 1+fltarr(n),key = keys,symsize = symsize,color = (findgen(n)+1)*254.0/n ELSE $
          schmidtlaw_res_obs_master_out,dir+"/schmidtlaw_res_obs"+res+".dat",outplot = outplot,thick = thicks,symbols = 1+fltarr(n),key = keys,symsize = symsize
        sk_files = [prefix + 'hAll.g1bwk.schmidtlaw_res_obs0.750000.dat']
        IF KEYWORD_SET(outplot) THEN $
          schmidtlaw_res_obs_master_out,sk_files,color = [50],/contour,outplot = outplot+'_CH_bin_',thick = 4+fltarr(n) ELSE $
          schmidtlaw_res_obs_master_out,sk_files,color = [50],/contour,thick = 4+fltarr(n)
        
        IF KEYWORD_SET(color) THEN $
          schmidtlaw_res_obs_master_out,dir+"/schmidtlaw_res_obs"+res+".dat",outplot = outplot,thick = thicks,symbols = 1+fltarr(n),key = keys,symsize = symsize,color = (findgen(n)+1)*254.0/n ELSE $
          schmidtlaw_res_obs_master_out,dir+"/schmidtlaw_res_obs"+res+".dat",outplot = outplot,thick = thicks,symbols = 1+fltarr(n),key = keys,symsize = symsize
        sk_files = [prefix + 'hAll.g1MbwK.schmidtlaw_res_obs0.750000.dat']
        IF KEYWORD_SET(outplot) THEN $
          schmidtlaw_res_obs_master_out,sk_files,color = [130],/contour,outplot = outplot+'_met_bin_',thick = 4+fltarr(n) ELSE $
          schmidtlaw_res_obs_master_out,sk_files,color = [130],/contour,thick = 4+fltarr(n)
                
        IF KEYWORD_SET(color) THEN $
          schmidtlaw_res_obs_master_out,dir+"/schmidtlaw_res_obs"+res+".dat",outplot = outplot,thick = thicks,symbols = 1+fltarr(n),key = keys,symsize = symsize,color = (findgen(n)+1)*254.0/n ELSE $
          schmidtlaw_res_obs_master_out,dir+"/schmidtlaw_res_obs"+res+".dat",outplot = outplot,thick = thicks,symbols = 1+fltarr(n),key = keys,symsize = symsize
        sk_files = [prefix + 'hAll.g14HBWK.schmidtlaw_res_obs0.750000.dat']
        IF KEYWORD_SET(outplot) THEN $
          schmidtlaw_res_obs_master_out,sk_files,color = [254],/contour,outplot = outplot+'_bin_',thick = 4+fltarr(n) ELSE $
          schmidtlaw_res_obs_master_out,sk_files,color = [254],/contour,thick = 4+fltarr(n)
    ENDIF
    IF NOT KEYWORD_SET(color) THEN loadct,0
ENDIF

;---------------- FFT -----------------
;/astro/users/christensen/code/IDL/HIcubes/fft_result.pro
IF plot_fft_result THEN BEGIN
    if keyword_set(outplot) THEN BEGIN
         IF KEYWORD_SET(color) THEN $
           fft_result,dir,keys = keys,formatthick = formatthick,label = label,symbols = psym,symsizes = symsizes,outplot = outplot,ctables = ctables,color = colors ELSE $
           fft_result,dir,keys = keys,formatthick = formatthick,label = label,symbols = psym,symsizes = symsizes,outplot = outplot
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(color) THEN $
          fft_result,dir,keys = keys,formatthick = formatthick,label = label,symbols = psym,symsizes = symsizes,ctables = ctables,color = colors ELSE $
          fft_result,dir,keys = keys,formatthick = formatthick,label = label,symbols = psym,symsizes = symsizes
    ENDELSE
ENDIF

;----------------- Surface Den -----------------
;/astro/users/christensen/code/IDL/MolecH/gas_surface_den.pro
IF plot_gas_surface_den THEN BEGIN
    IF KEYWORD_SET(outplot) THEN outfile = outplot + outext
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
;        IF useH2[i] eq 1 THEN BEGIN
            IF KEYWORD_SET(outplot) THEN gas_surface_den,dir[i] + '/' + filenames[i],useH2 = useH2[i],outplot = outfile[i],range = 12.0/2,color = 1 ELSE gas_surface_den,dir[i] + '/' + filenames[i],useH2 = useH2[i],range = 24.0/2,color = 1
;        ENDIF            
    ENDFOR
ENDIF

;---------------------------- Picture -------------------
IF plot_make_jpeg THEN BEGIN
    cameraFO = cameras - 1
    cameraEO = cameras + 1
    IF NOT KEYWORD_SET(outplot) THEN outfile = files ELSE outfile = outplot + outext
    FOR i = 0, n - 1 DO BEGIN
        print,files[i]
        cd,files[i] + '.' + halos_str[i]
        center_temp = [-1.0*center[0,i],center[1,i]] 
        make_jpeg,bands = [4,3,2],cam = cameraFO[i],outfile = outfile[i] + '_FO.jpeg',range = 24,resizefactor = resizefactor[i],center = center_temp,rotateAngle = rotateAngle[i]
        make_jpeg,bands = [4,3,2],cam = cameraEO[i],outfile = outfile[i] + '_EO.jpeg',range = 24,resizefactor = resizefactor[i],center = center_temp;,rotateAngle = rotateAngle[i]
    ENDFOR
ENDIF

;------------------------ Tracking the Mass ---------------
IF plot_track_mass THEN BEGIN
    track_mass_plot,dir,filenames,key = keys,linestyle = linestyles,thick = thicks
ENDIF

;------------------------ Cumlative Feedback ------------
IF plot_fbcum THEN BEGIN
    dirtop = dir
    base = filenames
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
        splitdir = strsplit(dir[i],'/')
        IF strpos(dir[i],'/steps') EQ -1 THEN dirtop[i] = strmid(dir[i],0,splitdir[N_ELEMENTS(splitdir)-1]) $
        ELSE dirtop[i] = strmid(dir[i],0,splitdir[N_ELEMENTS(splitdir)-2])
        splitbase = strsplit(filenames[i],'.')
        base[i] = strmid(filenames[i],0,splitbase[N_ELEMENTS(splitbase) - 1] - 1)
    ENDFOR
;    outflow_plots,dirtop,outplot = outplot, keys = keys, color = colors, thicks = thicks, linestyles = linestyles,label = label,ctables = ctables,yrange_fbcum = yrange_fbcum
    eject_plots,dirtop,outplot = outplot, keys = keys, color = colors, thicks = thicks,label = label,ctables = ctables,yrange_fbcum = [0,0.1],formatthick = formatthick,linestyles = linestyles;,/unscale
;    eject_plots,dirtop,outplot = outplot, keys = keys, color = colors, thicks = thicks,label = label,ctables = ctables,yrange_fbcum = [0,0.1],formatthick = formatthick,molecularH = molecularH,linestyles = linestyles,/massloading;,/unscale;/unscale;,/massloading
ENDIF

;------------------------- Radius where feedback comes from ---------
IF plot_fbrad THEN BEGIN
    lb_keys = strarr(n)
    lb_thicks = intarr(n)
    lb_color = intarr(n) + bgcolor
    lb_linestyle = intarr(n)

    IF KEYWORD_SET(outplot) AND KEYWORD_SET(normalize)     THEN device,filename = outplot + '_fbrad_norm.eps',/encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8
    IF KEYWORD_SET(outplot) AND NOT KEYWORD_SET(normalize) THEN device,filename = outplot + '_fbrad.eps',     /encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8

    FOR i = 0, N_ELEMENTS(files)-1 DO BEGIN
        splitdir = strsplit(dir[i],'/')
        IF strpos(dir[i],'/steps') EQ -1 THEN cd,strmid(dir[i],0,splitdir[N_ELEMENTS(splitdir)-1]) $
          ELSE cd,strmid(dir[i],0,splitdir[N_ELEMENTS(splitdir)-2])
        loadct,ctables[i]
        splitbase = strsplit(filenames[i],'.')
        base = strmid(filenames[i],0,splitbase[N_ELEMENTS(splitbase) - 1] - 1)

;        if i eq 1 THEN stop
        totalmass = 1
;        eject_gas_character,dirhead[i],ejecthistory,expellhistory,totalmass = totalmass
        ejecthistory =  mrdfits(dirhead[i] + 'grp1.eject_disk.fits',1)
        expellhistory = mrdfits(dirhead[i] + 'grp1.expell_disk.fits',1)
        IF 0 THEN BEGIN
         IF i EQ 0 AND     KEYWORD_SET(NORMALIZE) THEN histogramp,SQRT(ejecthistory.x*ejecthistory.x + ejecthistory.y*ejecthistory.y),weight = ejecthistory.mass,xrange = xrange_fb,/normalize,xtitle = 'Radius [kpc]',ytitle = '1/M dM/dr',                         max = xrange_fb[1],/nodata,nbins = 100
;histogramp,SQRT(ejecthistory.x*ejecthistory.x + ejecthistory.y*ejecthistory.y),weight = totalmass,xrange = xrange_fb,/normalize,xtitle = 'Radius [kpc]',ytitle = '1/M dM/dr',                         max = xrange_fb[1],/nodata,nbins = 100
        IF i EQ 0 AND NOT KEYWORD_SET(NORMALIZE) THEN histogramp,SQRT(ejecthistory.x*ejecthistory.x + ejecthistory.y*ejecthistory.y),xrange = xrange_fb,       xtitle = 'Radius [kpc]',ytitle = 'dN/dr',nbins = 100,max = xrange_fb[1],title = label,yrange = [0,4500],/nodata ;,xmargin = [16,3]

        IF     KEYWORD_SET(NORMALIZE) THEN histogramp,SQRT(ejecthistory.x*ejecthistory.x + ejecthistory.y*ejecthistory.y),weight = ejecthistory.mass,/overplot,normalize = normalize,color = colors[i],thick = thicks[i],max = xrange_fb[1],linestyle = linestyles[i],nbins = 100
;histogramp,SQRT(ejecthistory.x*ejecthistory.x + ejecthistory.y*ejecthistory.y),weight = totalmass/1e8,/overplot,normalize = normalize,color = colors[i],thick = thicks[i],max = xrange_fb[1],linestyle = linestyles[i],nbins = 100
        IF NOT KEYWORD_SET(NORMALIZE) THEN histogramp,SQRT(ejecthistory.x*ejecthistory.x + ejecthistory.y*ejecthistory.y),/overplot,                      color = colors[i],thick = thicks[i],max = xrange_fb[1],linestyle = linestyles[i],nbins = 100

        IF     KEYWORD_SET(NORMALIZE) THEN histogramp,SQRT(expellhistory.x*expellhistory.x + expellhistory.y*expellhistory.y),weight = ejecthistory.mass,/normalize,/overplot,color = colors[i],thick = thicks[i],max = xrange_fb[1],linestyle = 2,nbins = 100
        IF NOT KEYWORD_SET(NORMALIZE) THEN histogramp,SQRT(expellhistory.x*expellhistory.x + expellhistory.y*expellhistory.y),/overplot,                  color = colors[i],thick = thicks[i],max = xrange_fb[1],linestyle = 2,nbins = 100
    ENDIF ELSE BEGIN
        histogram_eject  = weighted_histogram(SQRT(ejecthistory.x*ejecthistory.x + ejecthistory.y*ejecthistory.y),    weight = ejecthistory.mass, nbins = 100,min = 0,max = xrange_fb[1],locations = locations)
        histogram_expell = weighted_histogram(SQRT(expellhistory.x*expellhistory.x + expellhistory.y*expellhistory.y),weight = expellhistory.mass,nbins = 100,min = 0,max = xrange_fb[1],locations = locations)        
        IF i EQ 0 THEN plot,locations,histogram_expell/histogram_eject,xtitle = 'Radius [kpc]',/nodata,yrange = [0,1]
        oplot,locations,histogram_expell/histogram_eject,color = colors[i],thick = thicks[i],linestyle = linestyles[i]
    ENDELSE
 
        IF KEYWORD_SET(keys) AND i lt N_ELEMENTS(keys) THEN BEGIN
            l_keys = lb_keys
            l_keys[i] = keys[i]
            l_thicks = lb_thicks
            l_thicks[i] = thicks[i]
            l_color = lb_color
            l_color[i] = colors[i]
            l_linestyle = lb_linestyle
            l_linestyle[i] = linestyles[i]
            legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/top,/right,box = 0
        ENDIF
        print,MAX(weighted_histogram(SQRT(ejecthistory.x*ejecthistory.x + ejecthistory.y*ejecthistory.y),weight = totalmass/1e8,max = xrange_fb[1],nbins = 100))

        stop

        undefine,totalmass
        undefine,ejecthistory
        undefine,expellhistory
    ENDFOR    
    IF KEYWORD_SET(outplot) THEN device,/close
ENDIF

;------------------------- metallicity of feedback gas ---------
IF plot_fbz THEN BEGIN
    lb_keys = strarr(n)
    lb_thicks = intarr(n)
    lb_color = intarr(n) + bgcolor
    lb_linestyle = intarr(n)
    zsolar = 0.0130215

    IF KEYWORD_SET(outplot) THEN device,filename = outplot + '_fbrad.eps',     /encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8

    FOR i = 0, N_ELEMENTS(files)-1 DO BEGIN
        splitdir = strsplit(dir[i],'/')
        IF strpos(dir[i],'/steps') EQ -1 THEN cd,strmid(dir[i],0,splitdir[N_ELEMENTS(splitdir)-1]) $
          ELSE cd,strmid(dir[i],0,splitdir[N_ELEMENTS(splitdir)-2])
        loadct,ctables[i]
        splitbase = strsplit(filenames[i],'.')
        base = strmid(filenames[i],0,splitbase[N_ELEMENTS(splitbase) - 1] - 1)

        totalmass = 1
;        eject_gas_character,dirhead[i],ejecthistory,expellhistory,totalmass = totalmass
        ejecthistory =  mrdfits(dirhead[i] + 'grp1.eject_disk.fits',1)
        expellhistory = mrdfits(dirhead[i] + 'grp1.expell_disk.fits',1)
        xrange_z= [0,0.4]
        IF i EQ 0 THEN histogramp,ejecthistory.metallicity/zsolar,xtitle = 'Z/Z'+sunsymbol(),ytitle = 'dN/dr',/nodata,nbins = 100,title = label,max = xrange_z[1],xrange = xrange_z ;,xmargin = [16,3]
        histogramp,               ejecthistory.metallicity/zsolar,/overplot,color = colors[i],thick = thicks[i],max = xrange_z[1],linestyle = linestyles[i],nbins = 100
        histogramp,               ejecthistory.metallicity/zsolar,/overplot,color = colors[i],thick = thicks[i],max = xrange_z[1],linestyle = 2,nbins = 100

        IF KEYWORD_SET(keys) AND i lt N_ELEMENTS(keys) THEN BEGIN
            l_keys = lb_keys
            l_keys[i] = keys[i]
            l_thicks = lb_thicks
            l_thicks[i] = thicks[i]
            l_color = lb_color
            l_color[i] = colors[i]
            l_linestyle = lb_linestyle
            l_linestyle[i] = linestyles[i]
            legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/top,/right,box = 0
        ENDIF
        print,MAX(weighted_histogram(SQRT(ejecthistory.x*ejecthistory.x + ejecthistory.y*ejecthistory.y),weight = totalmass/1e8,max = xrange_fb[1],nbins = 100))

        undefine,totalmass
        undefine,ejecthistory
        undefine,expellhistory
    ENDFOR    
    IF KEYWORD_SET(outplot) THEN device,/close
    stop
ENDIF

;----------------------- Angular Momentum -------------------
IF plot_j_test THEN BEGIN
    IF keyword_set(outplot) THEN device, filename=outplot + '_jjtot.eps', /encapsulated, /color,/times,ysize = ysize,xsize = xsize,bits_per_pixel = 8  $
    ELSE window,0;,ysize = ysize,xsize = xsize
    s =indgen(n_elements(h))*0.1
    xsi = 1.0-1.25*(1.0-(1.25-1.0)*alog(1.25/(1.25-1.0)))
    ps = xsi*1.25*(1.25-1.0)/(xsi*s+1.25-1.0)^2.
    loadct,39

    lb_keys = strarr(n)
    lb_thicks = intarr(n)
    lb_color = intarr(n) + bgcolor
    lb_linestyle = intarr(n)  
 
    FOR i = 0, n -1 DO BEGIN
        cd,dir[i]
        loadct,ctables[i]
        rmax = opticalRadii(filename = filenames[i],center = center[*,i],extno = cameras[i],verbose = verbose)
        print,filenames[i],rmax
;posx=[.15,.95] 
        ytickname=[' ','0.2',' ','0.6',' ','1.0'] 
;posy=[.25,.95]
        jdist, prefix=dir[i] + '/', filenames[i], h, hs, hg, d, s, g, fdisk, dmass, smass, gmass, rmax=rmax, /obs 
;        if i eq 0 then plot, indgen(n_elements(h))*0.05, dmass/max(dmass), xrange=[0,2.5], yrange=[0,1], pos=[posx[0],posy[0],posx[1],posy[1]], xtickname=xtickname, ytickname=ytickname,xtitle = textoidl('s = j/j_{tot}'),ytitle = 'P(s)',/nodata,title = label;, charthick=2
        if i eq 0 then plot, indgen(n_elements(h))*0.05, dmass/max(dmass), xrange=[0,2.5], yrange=[0,1], xtickname=xtickname,xtitle = 's = j/<j>',ytitle = 'P(s)/Max(P(s))',/nodata,title = label;, charthick=2
        oplot,indgen(n_elements(h))*0.05,dmass/max(dmass), color = colors[i],thick = thicks[i],linestyle = 2
        if n_elements(hs) lt n_elements(hg) then b = n_elements(hg) else b = n_elements(hs)
        oplot, indgen(b)*0.05, ((gmass + smass)/max(gmass + smass))*(fdisk/.2121), linestyle=0, color = colors[i],thick = thicks[i]
;        IF KEYWORD_SET(keys) THEN BEGIN
;        IF 0 THEN BEGIN
;            l_keys = lb_keys
;            l_keys[i] = keys[i]
;            l_thicks = lb_thicks
;            l_thicks[i] = thicks[i]
;            l_color = lb_color
;            l_color[i] = colors[i]
;            l_linestyle = lb_linestyle
;            l_linestyle[i] = 0 ;linestyle[i]
;            legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/right,box = 0
;        ENDIF
    ENDFOR
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = fltarr(n_elements(keys)),thick = thicks,ctables = ctables,/right,box = 0; ELSE 
    legend, ['Dark Matter', 'Disk Baryons'], linestyle=[2,0],box = 0,thick = [0,0] + thicks[0],pos = [2.5,0.35],/right;pos = [0.948205,0.428272 - 0.05] ;charthick = 2
    if keyword_set(outplot) then device, /close
ENDIF

;-------------------------- Contour of Star Forming Gas ----
IF plot_sfgas THEN BEGIN
    sfgas_contour,dirhead,fileprefixes,molecularH = useH2,outplot = outplot, keys = keys, colors = colors, ctables = ctables, formatthick = formatthick, label = label, thicks = thicks
;    IF keyword_set(outplot) THEN sfgas_contour,dirhead[2],fileprefixes[2],/molecularH,outplot = outplot + '2',/zdivide,min1 = 1.5,max1 = 3.5,min2 = 1, max2 = 3,nlevels = 15,thick =[2],maxlevel = 6e8 ELSE sfgas_contour,dirhead[2],fileprefixes[2],/molecularH,/zdivide,min1 = 1.5,max1 = 3.5,min2 = 1, max2 = 3,nlevels = 15,thick =[2],maxlevel = 6e8
ENDIF
stop
end
