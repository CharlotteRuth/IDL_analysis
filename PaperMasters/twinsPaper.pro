;.r opticalRadii¥
;.r densRadius
;.r photometricProf
;.r twinsPaper

;outplot = '~/plots/twins/twins_color'
;outplot = '~/plots/twins/twins2_color'
;outplot = '~/plots/twins/twins.talk_color'
;outplot = '~/Plots/twins/h603_h516_0_talk'
pro twinsPaper,outplot = outplot,color = color,verbose = verbose,thick = thick,psym = psym
spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

IF hostname EQ 'ozma' THEN idldir = '/home/christensen/Code/IDL/' ELSE idldir = '/astro/users/christensen/code/IDL/'
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL

plot_vcirc = 1
plot_photometricProf = 0
plot_densRadius = 0
plot_fft_result = 0
plot_gas_surface_den = 0
plot_make_jpeg = 0
plot_schmidtlaw_global_obs = 0
plot_schmidtlaw_res_obs = 0
plot_schmidtlaw_res_obs_master_out = 0
plot_j_test = 0

steps = ['00100','00512']
steps = ['00072','00512']
halos = [2,1]

steps = ['00512','00512']
halos = [1,1]

halos_str = strtrim(halos,2)
dir = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[0] + '.dir',$
       prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir']
files = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.' + steps[0] + '.dir/h603.cosmo50cmb.3072g14HBWK.' + steps[0],$
         prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.dir/h516.cosmo25cmb.3072g14HBWK.' + steps[1]]
pfiles = [prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param' ,$
          prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
filebases = ['h603.cosmo50cmb.3072g14HBWK.' + steps[0] + '.halo.' + halos_str[0], $
             'h516.cosmo25cmb.3072g14HBWK.' + steps[1] + '.halo.' + halos_str[1]]
filenames = ['h603.cosmo50cmb.3072g14HBWK.' + steps[0], $
             'h516.cosmo25cmb.3072g14HBWK.' + steps[1]]
;keys = ['L'+textoidl('^* Galaxy, z = 3.4'),'Dwarf Galaxy, z = 0']
useH2 = [1,1]
distunits = [50000.,25000]
massunits  = [1.84793e16,2.310e15]
filternums = [13,13]
cameras = [14,14]
center = [[0,0],$
          [0,0]]
yrange_vcirc = [0,120];75]
yrange_photo = [26,17]
maxdistance_photo = 4
maxdistance = 6
rotateAngle = [3,3]
resizefactor = [4,4]
;resizefactor = [8,8]
boxside = [4,3]
nocolorbar = [0,1]
intq = [33,32]
outext = ['h603','h516']
colors = [50,245]
colors = [30,245]
;colors = [60,220]

n = N_ELEMENTS(files)
formatplot,outplot = outplot,thick = thick
IF KEYWORD_SET(outplot) THEN BEGIN
    fgcolor = 0 
    thickbase = 4
ENDIF ELSE BEGIN
    fgcolor = 255
    thickbase = 2
ENDELSE
IF KEYWORD_SET(thick) THEN BEGIN
    thick = 8
ENDIF
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    obscolor = fgcolor
    if NOT keyword_set(colors) then  colors  = (findgen(n) + 1)*240/n else colors = colors
    IF NOT KEYWORD_SET(psym) THEN psym = fltarr(n) + 4
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + thickbase
    IF NOT KEYWORD_SET(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT keyword_Set(symsizes) THEN symsizes = fltarr(n) + 2
ENDIF ELSE BEGIN
    loadct,0    
    obscolor = 100
    colors = (fltarr(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(broadband)) + 5
    IF NOT KEYWORD_SET(psym) THEN  psym = (findgen(n)+2)*2
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + thickbase ;thicks = (findgen(n) + 1)*6/n - 1
    IF NOT KEYWORD_SET(linestyles) THEN linestyles = REVERSE(findgen(n)*2) 
    IF NOT keyword_Set(symsizes) THEN symsizes = fltarr(n) + 2
ENDELSE

;--------------------------- Vcirc -------------
IF plot_vcirc THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,outfile = outplot + '_vcirc.eps',maxdistance = 8,vfinal = vfinal,color = colors,type = 'all' ELSE vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc,outfile = outplot + '_vcirc.eps',maxdistance = 8,vfinal = vfinal,type = 'all' 
    ENDIF ELSE BEGIN
        if keyword_set(color) THEN vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = 8,vfinal = vfinal,color = colors,type = 'all'  ELSE   vcirc,files + '.halo.' + halos_str,massunits,distunits,keys = keys,thicks = thicks,linestyle = linestyles,yrange = yrange_vcirc, maxdistance = 8,vfinal = vfinal,type = 'all' 
    ENDELSE
    print,vfinal
ENDIF

;--------------------------- Photometric Profile ------------
IF plot_photometricProf THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN 
    IF KEYWORD_SET(color) THEN photometricProf,files + '.' + halos_str + '/broadband.fits',keys = keys,thicks = thicks,linestyle = linestyles,filternums = filternums,cameras = cameras,yrange = yrange_photo,outplot = outplot + '_photoProf.eps',bands = fltarr(N_ELEMENTS(files)) + 5,maxdistance = maxdistance_photo, recenter = recenter,colors = colors ELSE photometricProf,files + '.' + halos_str + '/broadband.fits',keys = keys,thicks = thicks,linestyle = linestyles,filternums = filternums,cameras = cameras,yrange = yrange_photo,outplot = outplot + '_photoProf.eps',bands = fltarr(N_ELEMENTS(files)) + 5,maxdistance = maxdistance_photo, recenter = recenter
      ENDIF ELSE BEGIN 
        IF KEYWORD_SET(color) THEN photometricProf,files + '.' + halos_str + '/broadband.fits',keys = keys,thicks = thicks,linestyle = linestyles,filternums = filternums,cameras = cameras,yrange = yrange_photo,bands = fltarr(N_ELEMENTS(files)) + 5,maxdistance = maxdistance_photo,recenter = recenter,colors = colors ELSE photometricProf,files + '.' + halos_str + '/broadband.fits',keys = keys,thicks = thicks,linestyle = linestyles,filternums = filternums,cameras = cameras,yrange = yrange_photo,bands = fltarr(N_ELEMENTS(files)) + 5,maxdistance = maxdistance_photo,recenter = recenter
    ENDELSE
ENDIF

;----------------------------------- Density Radius ----------
;/astro/users/christensen/code/IDL/HIcubes/densRadius.pro
IF plot_densRadius THEN BEGIN
    if keyword_set(outplot) THEN BEGIN
        IF KEYWORD_SET(color) THEN densRadius,files,distunits,massunits,maxdistance = maxdistance_photo,outplot = outplot,color = colors,halos_str = halos_str ELSE densRadius,files,distunits,massunits,maxdistance = maxdistance_photo,outplot = outplot,halos_str = halos_str
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(color) THEN densRadius,files,distunits,massunits,maxdistance = maxdistance_photo,color = colors,halos_str = halos_str  ELSE densRadius,files,distunits,massunits,maxdistance = maxdistance_photo,halos_str = halos_str
    ENDELSE
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

;----------------- Surface Den -----------------
;/astro/users/christensen/code/IDL/MolecH/gas_surface_den.pro
IF plot_gas_surface_den THEN BEGIN
    IF KEYWORD_SET(outplot) THEN outfile = outplot + outext
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
        IF useH2[i] eq 1 THEN BEGIN
            IF KEYWORD_SET(outplot) THEN gas_surface_den,dir[i] + '/' + filenames[i],useH2 = useH2[i],outplot = outfile[i],color = 1,halo_str = halos_str[i],range = maxdistance_photo, boxside = boxside[i], nocolorbar = nocolorbar[i] ELSE gas_surface_den,dir[i] + '/' + filenames[i],useH2 = useH2[i],outplot = outplot,color = 1,halo_str = halos_str[i],range = maxdistance_photo, boxside = boxside[i], nocolorbar = nocolorbar[i]
            ENDIF
    ENDFOR
ENDIF

;---------------------------- Picture -------------------
IF plot_make_jpeg THEN BEGIN
    cameraFO = cameras
    cameraEO = cameras + 2
    IF NOT KEYWORD_SET(outplot) THEN outfile = files ELSE outfile = outplot + outext
    FOR i = 0, n - 1 DO BEGIN
        print,files[i]
        cd,files[i] + '.' + halos_str[i]
        center_temp = [-1.0*center[0,i],center[1,i]] 
        make_jpeg,bands = [4,3,2],scales = [5.0,5.0,5.0],cam = cameraFO[i],outfile = outfile[i] + '_FO.jpeg',rotateAngle = rotateAngle[i],center = center_temp,range = 2.00*maxdistance_photo,resizefactor = resizefactor[i];range
        make_jpeg,bands = [4,3,2],scales = [5.0,5.0,5.0],cam = cameraEO[i],outfile = outfile[i] + '_EO.jpeg',range = 2.0*maxdistance_photo,resizefactor = resizefactor[i];,rotateAngle = rotateAngle,center = center_temp,range = 24;range
    ENDFOR
ENDIF
;---------------------------------- SK -global ----------------
IF  plot_schmidtlaw_global_obs THEN BEGIN
    data=fltarr(2,N_ELEMENTS(dir))
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
        print,dir[i]
        data[*,i] = schmidtlaw_global_obs(dir[i],filenames[i],pfiles[i],useH2 = useH2[i],/Halpha,intq = intq[i], camera = cameras[i],center = center[*,i],verbose = verbose,halo_str = halos_str[i]) ;,tipsyfile = files[i]);,tipsyfile = filebases[i]
    ENDFOR 

    if (KEYWORD_SET(outplot)) then begin
        device,filename = outplot + '_SK.eps',/color,bits_per_pixel= 8,/times,ysize=18,xsize=18,xoffset =  2,yoffset =  2
    endif else begin
        set_plot,'x'
        window,1
    endelse
    IF keyword_set(color) THEN loadct,39
    xsigma = 10.0^(findgen(600)/100 - 1.) ;xsigmalow
    ysigma=2.5e-4*xsigma^1.4
    
    readcol,'~/code/HIcubes/ks98.dat',name,D,logHI,logH2,logH,logSFR,tdyn,co,HI,halpha ;,format='(A9D)'
    readcol,'~/code/HIcubes/uavpl.dat',logHdarf,logSFRdwarf
    plot,alog10(xsigma),alog10(ysigma),ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-0.5,2.5], yrange = [-4,0], xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]"   
    oplot,logH,logSFR,psym = 2,color = obscolor
    oplot,logHdarf,logSFRdwarf,psym = 1,color = obscolor,symsize = 2 
    FOR i = 0, N_ELEMENTS(dir) - 1 DO oplot,[alog10(data[0,i]),alog10(data[0,i])],[alog10(data[1,i]),alog10(data[1,i])],psym = symcat(psym[i]),color = colors[i],symsize = symsizes[i]
    oplot,[-0.2,0.2],[-1,-1]
    oplot,[0,0],[-0.8,-1.2]
    IF KEYWORD_SET(keys) THEN legend,[keys,'Kennicutt 98'],psym = [symcat(psym),2],color = [colors,obscolor],/bottom,/right
    if (KEYWORD_SET(outplot)) then device,/close else stop
ENDIF

;-------------------------- Res SK law Prep ----------------------------
IF plot_schmidtlaw_res_obs THEN BEGIN
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
        cd,dir[i]
        center_temp = [-1.0*center[0,i],center[1,i]] 
        schmidtlaw_res_obs,filenames[i],pfiles[i],useH2 = useH2[i],extno = cameras[i],center = center_temp,rotateAngle = rotateAngle[i],verbose = verbose,halo_str = halos_str[i]
    ENDFOR
ENDIF

;------------- Resolved K-S law plot ----------------------
;~/code/HIcubes/schmidtlaw_res_obs.pro
res = '0.750000'
IF plot_schmidtlaw_res_obs_master_out THEN BEGIN
;    IF KEYWORD_SET(color) THEN $
;      schmidtlaw_res_obs_master_out,dir+"/schmidtlaw_res_obs"+res+".dat",outplot = outplot,thick = thicks,symbols = 1+fltarr(n),key = keys,symsize = symsize,color = (findgen(n)+1)*240.0/n ELSE $
;      schmidtlaw_res_obs_master_out,dir+"/schmidtlaw_res_obs"+res+".dat",outplot = outplot,thick = thicks,symbols = 1+fltarr(n),key = keys,symsize = symsize

    loadct,39
    sk_files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00072.dir/schmidtlaw_res_obs0.750000.dat',$
                '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/schmidtlaw_res_obs0.750000.dat']
    IF KEYWORD_SET(outplot) THEN $
      schmidtlaw_res_obs_master_out,sk_files,color = colors,outplot = outplot+'_bin_',thick = 4+fltarr(n),true = true,symbols = psym ELSE $
      schmidtlaw_res_obs_master_out,sk_files,color = colors,thick = 4+fltarr(n),true = true,symbols = psym

    IF NOT KEYWORD_SET(color) THEN loadct,0
ENDIF

;---------------J Dist ------------------------------
IF plot_j_test THEN BEGIN
    if keyword_set(outplot) THEN BEGIN
;       formatplot,/outplot
        xsize = 18 ;10*n
        ysize = 12 ;12
        mxTitSize = 1.5
        mxTitOffset = 2
        device, filename=outplot + '_jjtot.eps', /encapsulated, /color,/times,ysize = ysize,xsize = xsize,bits_per_pixel = 8 
    ENDIF ELSE BEGIN
        xsize = 400 ;*n
        ysize = 350 ;475
        mxTitSize = 1.5
        mxTitOffset = 1
        window,0,ysize = ysize,xsize = xsize
    ENDELSE
;    s =indgen(n_elements(h))*0.1
;    xsi = 1.0-1.25*(1.0-(1.25-1.0)*alog(1.25/(1.25-1.0)))
;    ps = xsi*1.25*(1.25-1.0)/(xsi*s+1.25-1.0)^2.
    loadct,39
   
    FOR i = 0, n -1 DO BEGIN
        cd,dir[i]
        rmax = opticalRadii(filename = filenames[i],center = center[*,i],extno = cameras[i],halo_str = halos_str[i])
        print,rmax
        posx=[.15,.95] 
        ytickname=[' ','0.2',' ','0.6',' ','1.0'] 
        posy=[.25,.95]
        jdist, prefix=dir[i] + '/', filenames[i], h, hs, hg, d, s, g, fdisk, dmass, smass, gmass, rmax=rmax, /obs 
        if i eq 0 then plot, indgen(n_elements(h))*0.05, dmass/max(dmass), xrange=[0,2.5], yrange=[0,1], pos=[posx[0],posy[0],posx[1],posy[1]], xtickname=xtickname, ytickname=ytickname,xtitle = textoidl('s = j/j_{tot}'),ytitle = 'P(s)',/nodata,title = label;, charthick=2
        oplot,indgen(n_elements(h))*0.05,dmass/max(dmass), color = colors[i],thick = thicks[i],linestyle = 2
        if n_elements(hs) lt n_elements(hg) then b = n_elements(hg) else b = n_elements(hs)
        oplot, indgen(b)*0.05, (gmass/max(gmass)+smass/max(smass))*(fdisk/.2121), linestyle=0, color = colors[i],thick = thicks[i]
    ENDFOR
;    legend, ['DM', 'Observable'], linestyle=[0,2], /top, /right ;charthick = 2
    if keyword_set(outplot) then device, /close
ENDIF

end
