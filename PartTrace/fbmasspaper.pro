PRO fbmasspaper,formatthick = formatthick,outplot = outplot,colors = colors
formatplot,outplot = outplot,thick = formatthick
;IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
;IF KEYWORD_SET(outplot) THEN bgcolor = 255 ELSE bgcolor = 0

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN BEGIN 
   prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' 
   plotsdir = '~/Plots/'
ENDIF ELSE BEGIN
   prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
   plotsdir = '~/plots/'
ENDELSE

plot_reeject_r = 1
plot_eject_v_mass = 0
plot_times_ejected = 0
plot_times_expelled = 0
plot_times_accrdisk = 0
x = 1

CASE x OF
   1: BEGIN
      dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
      file799 = 'h799.cosmo25cmb.3072g14HBWK'
      dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
      file516 = 'h516.cosmo25cmb.3072g14HBWK'
      dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
      file986 = 'h986.cosmo50cmb.3072g14HBWK'
      dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
      file603 = 'h603.cosmo50cmb.3072g14HBWK'
      dir277 = prefix + 'h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/'
      file277 = 'h277.cosmo50cmb.3072g14HMbwK'
      dirs    = [dir799,dir516,dir516,dir986,dir986,dir986,dir986,dir603,dir603,dir603,dir277,dir277]
      files    = [file799,file516,file516,file986,file986,file986,file986,file603,file603,file603,file277,file277]
      haloid = ['1'   ,'1'   ,'2'   ,'1'   ,'2'   ,'3'   ,'9'   ,'1'   ,'2'    ,'3'  ,'1'   ,'2']
;      dirs    = [dir799,dir516,dir516,dir986,dir986,dir986,dir986,dir603,dir603,dir603,dir277]
;      files    = [file799,file516,file516,file986,file986,file986,file986,file603,file603,file603,file277]
;      haloid = ['1'   ,'1'   ,'2'   ,'1'   ,'2'   ,'3'   ,'9'   ,'1'   ,'2'    ,'3'  ,'2']
   END
ENDCASE

IF plot_eject_v_mass  THEN eject_v_mass, dirs,files,halo = haloid,colors = colors,outplot = outplot

IF plot_times_ejected THEN times_ejected,dirs,files,halo = haloid,colors = colors,outplot = outplot
IF plot_times_expelled THEN times_ejected,dirs,files,halo = haloid,colors = colors,outplot = outplot,/expelled
IF plot_times_accrdisk THEN times_ejected,dirs,files,halo = haloid,colors = colors,outplot = outplot,/accrdisk

IF plot_reeject_r THEN plot_half_eject,dirs,files,halo = haloid,colors = colors,outplot = outplot,/normalize,/stellarmass
stop
END
