;h603_h516_h986_h799evol_talk
pro twinsEvolPaper,formatthick = formatthick,outplot = outplot
formatplot,outplot = outplot,thick = formatthick
IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(outplot) THEN bgcolor = 255 ELSE bgcolor = 0

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN BEGIN 
   prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' 
   plotsdir = '~/Plots/'
ENDIF ELSE BEGIN
   prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
   plotsdir = '~/plots/'
ENDELSE

plot_densR = 0
plot_massEvol = 0
plot_fbcum = 0 ;checked
plot_fbmassload = 0 ;checked
plot_fbcumratio = 0 ;checked
plot_inoutflow = 1 ;checked
plot_sfr = 0
plot_resks = 0
plot_fbratio = 0
plot_eject_radius = 0
plot_eject_metallicity = 0

x = 3
CASE x OF
   1: BEGIN
      dir      = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/',$
                  prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']
   END
   2: BEGIN
      dir      = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']
   END
   3: BEGIN
      dir      = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']
      stepdir  = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir'] 
      basedir  = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/']
      files    = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512.halo.1',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
      pfiles   = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
      skfilename = stepdir + '/schmidtlaw_res_obs0.750000.dat'
       haloid = ['1','1']
;     haloid = ['2','2']
      maxdistance_photo = 6
      scale = [1,8]
      colors = [30,245]
      linestyle = [0,0]
      psym = [14,14]
      keys = ['Spiral'+textoidl(' Galaxy'),'Dwarf Galaxy']
      colors_sk = colors
      psym_sk   = psym
      keys_sk   = keys
      IF keyword_set(outplot) THEN outplot = plotsdir + 'twins/h603_h516_evol_talk'
   END
   4: BEGIN
      dir      = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/']
      stepdir  = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir'] 
      basedir  = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/']
      files    = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.halo.1',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
      pfiles   = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
      skfilename = stepdir + '/schmidtlaw_res_obs0.750000.dat'
      maxdistance_photo = 6
      scale = [1,8]
      colors = [30,245]
      linestyle = [0,0]
      psym = [14,14]
      keys = ['Spiral'+textoidl(' Galaxy'),'Dwarf Galaxy']
      colors_sk = colors
      psym_sk   = psym
      keys_sk   = keys
      IF keyword_set(outplot) THEN  outplot = plotsdir + 'twins/h986_h516_evol_talk'
   END
   5: BEGIN
      dir      = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/']
      stepdir  = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir',$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir'] 
      basedir  = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/',$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/']
      files    = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.halo.1',$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512.halo.1']
      pfiles   = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/h799.cosmo25cmb.3072g14HBWK.param']
      skfilename = stepdir + '/schmidtlaw_res_obs0.750000.dat'
      maxdistance_photo = 6
      scale = [1,8]
      colors = [30,245]
      linestyle = [0,0]
      psym = [14,14]
      keys = ['Spiral'+textoidl(' Galaxy'),'Dwarf Galaxy']
      colors_sk = colors
      psym_sk   = psym
      keys_sk   = keys
      IF keyword_set(outplot) THEN outplot = plotsdir + 'twins/h986_h799_evol_talk'
   END
   6: BEGIN
      dir      = [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']
      files    = [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
      pfiles   = [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
      scale    = [1]
      colors   = [245]
   END
   7: BEGIN
      dir      = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/']
      stepdir  = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir',$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir'] 
      basedir  = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/',$ 
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/']
      files    = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512.halo.1',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.halo.1',$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512.halo.1']
      pfiles   = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
                  prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param',$
                  prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param',$
                  prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/h799.cosmo25cmb.3072g14HBWK.param']
      skfilename =prefix +['schmidtlaw_res_obs0.750000.h986.h603.dat','schmidtlaw_res_obs0.750000.h516.h799.dat']
      maxdistance_photo = 6
      scale = [1,8,1,8]
      colors = [30,30,245,245]
;      colors = [245,30,245,30]
      linestyle = [0,2,0,2]
;      IF keyword_set(formatthick) THEN thicks = [8,3,8,3] ELSE thicks = [3,1,3,1]
      IF keyword_set(formatthick) THEN thicks = [8,8,8,8] ELSE thicks = [3,3,3,3] 
      psym = [14,16,14,16]
      keys = ['Spiral A','Dwarf A','Spiral B','Dwarf B']
      colors_sk = [50,120]
      colors_sk = [245,30]
      psym_sk   = [16,16]
      keys_sk   = ['Spirals','Dwarfs']
      IF keyword_set(outplot) THEN outplot = plotsdir + 'twins/h603_h516_h986_h799evol_talk'
   END
ENDCASE

IF plot_densR THEN BEGIN 
   kpcunit  = fltarr(n_elements(dir))
   msununit = fltarr(n_elements(dir))
   filesall = files
   FOR i = 0, n_elements(dir) - 1 DO BEGIN
      periodpos = strsplit(files[i],'.')
      filesall[i] = strmid(files[i],0,periodpos[n_elements(periodpos) - 2] - 1)
      units = tipsyunits(pfiles[i])
      kpcunit[i]  = units.lengthunit
      msununit[i] = units.massunit
   ENDFOR
   densRadius,filesall,kpcunit,msununit,maxdistance = maxdistance_photo,outplot = outplot
ENDIF

IF plot_massEvol THEN BEGIN
   mass_evol,dir,files,pfiles,scale = scale,linestyle = linestyle,thicks = thicks,colors = colors,psym = psym,keys = keys,outplot = outplot,formatthick = formatthick
   stop
ENDIF

IF plot_fbcum THEN BEGIN
;    eject_plots,basedir,outplot = outplot,keys = keys,colors = colors,thicks = thicks,formatthick = formatthick,linestyles = linestyle,/redshift
    eject_plots,basedir,outplot = outplot,keys = keys,colors = colors,thicks = thicks,formatthick = formatthick,linestyles = linestyle,/redshift,/unscale
    stop
ENDIF

IF plot_fbmassload THEN BEGIN
   eject_plots,basedir,outplot = outplot,keys = keys,colors = colors,thicks = thicks,formatthick = formatthick,linestyles = linestyle,/redshift,/massloading
   stop
ENDIF

IF plot_fbcumratio THEN BEGIN
   eject_plots,basedir,outplot = outplot,keys = keys,colors = colors,thicks = thicks,formatthick = formatthick,linestyles = linestyle,/redshift,/ratio
   stop
ENDIF

IF plot_inoutflow THEN BEGIN
   plot_inflow_outflow,basedir,haloid = haloid,scale = scale,outplot = outplot,keys = keys,colors = colors,thicks = thicks,formatthick = formatthick,linestyle = linestyle
   stop
ENDIF

IF plot_sfr THEN BEGIN
    sfr_twins,files,pfiles,scale = scale,linestyle = linestyle,outplot = outplot,keys = keys,colors = colors,thick = thicks,formatthick = formatthick
    stop
ENDIF

IF plot_resks THEN BEGIN
   schmidtlaw_res_obs_master_out,skfilename,outplot = outplot,color = colors_sk,symbols = psym_sk,key = keys_sk,symsize = symsize_sk,/scaleHe,/contour,formatthick = formatthick
ENDIF

IF plot_fbratio THEN BEGIN
    eject_radius,basedir,files,scalem = scale,scaler = scale^(1.0/3),outplot = outplot,keys = keys,colors = colors,thicks = thicks,formatthick = formatthick,linestyles = linestyle,/ratio
ENDIF

IF plot_eject_radius THEN BEGIN
   eject_radius,basedir,files,scalem = scale,scaler = scale^(1.0/3),outplot = outplot,keys = keys,colors = colors,thicks = thicks,formatthick = formatthick,linestyles = linestyle
ENDIF

IF plot_eject_metallicity THEN BEGIN
   eject_metallicity,basedir,files,scale = scale,outplot = outplot,keys = keys,colors = colors,thicks = thicks,formatthick = formatthick,linestyles = linestyle
ENDIF
END
