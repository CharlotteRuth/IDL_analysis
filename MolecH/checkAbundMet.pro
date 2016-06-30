
pro checkAbundMet_master, outplot = outplot,color = color
spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/12M/' ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/12M/'
dir = prefix + 'Disk_Collapse_1e6/'
files = ['Disk_Collapse_1e6.00100.scalez1.00010','Disk_Collapse_1e6.00100.scalez0.3.00010','Disk_Collapse_1e6.00100.scalez0.1.00010']
;iles = ['Disk_Collapse_1e6.00100.scalez1.00020','Disk_Collapse_1e6.00100.scalez0.3.00020','Disk_Collapse_1e6.00100.scalez0.1.00020']
msol_per_sysmass = 2.362e5
kpc_per_syslength = 1.0
title = ['Z' + sunsymbol(),'0.3 Z' + sunsymbol(),'0.1 Z' + sunsymbol()]
if keyword_set(outplot) AND NOT keyword_set(color) Then outfile = outplot ;'/astro/store/student-scratch1/christensen/MolecH/12M/Disk_Collapse_'
if keyword_set(outplot) AND keyword_set(color) Then outfile = outplot ;'/astro/store/student-scratch1/christensen/MolecH/12M/Disk_Collapse_color'

checkAbundMet_Obs,dir + files,msol_per_sysmass,kpc_per_syslength,title = title,outplot = outfile,color = color
;checkAbundMet,dir + files,msol_per_sysmass,kpc_per_syslength,title = title,outplot = outfile,color = color
end

pro checkAbundMet,files,msol_per_sysmass,kpc_per_syslength,outplot = outplot,title = title,color = color
n = N_ELEMENTS(files)

formatplot,outplot = outplot

IF KEYWORD_SET(outplot) THEN BEGIN
    fgcolor = 0
    l_charsize = 0.75
    set_plot,'ps'
    device,filename = outplot+'metabundSD.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
ENDIF ELSE BEGIN
    l_charsize = 1.0
    fgcolor = 255
    set_plot,'x'
    window,1,xsize = 712,ysize = 392
ENDELSE
IF NOT KEYWORD_SET(title) THEN title = strind(n)
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
dens_convert =  msol_per_sysmass*gm_per_msol*amu_per_gm/kpc_per_syslength^3/cm_per_kpc^3 ;massunit*amu_per_gm/molec_weight/lengthunit^3 ;Converts to grams/cm^3

;!p.multi = [1,n,0]
multiplot,[n,1];mxtitle = textoidl(' \Sigma_{HI + H_{2}} [amu/cm^2]'),mxTitSize = 1.5,mxTitOffset = 2;,/doxaxis,/doyaxis
FOR i = 0, n - 1 DO BEGIN
    rtipsy,files[i],h,g,d,s
    readarr,files[i]+'.H2',h,H2,/ascii,part = 'gas'
    readarr,files[i]+'.HI',h,HI,/ascii,part = 'gas'
    readarr,files[i]+'.smoothlength',h,hl,/ascii,part = 'gas'

    gamu = g
    gamu.dens = g.dens*dens_convert/h.time/h.time/h.time
    hl = hl*cm_per_kpc*kpc_per_syslength*h.time

    IF i eq 0 THEN plot,gamu.dens*hl*(2.0*H2 + HI),2.0*H2/(2.0*H2 + HI),psym = 3,/ylog,/xlog,xrange = [1e19,1e23],yrange = [1e-6,1],title = title[i],ytitle = textoidl('H_2/H'),xtickname = [textoidl('10^{19}'),textoidl('10^{20}'),textoidl('10^{21}'),textoidl('10^{22}'),textoidl('10^{23}')],ytickname = [textoidl('10^{-6}'),textoidl('10^{-5}'),textoidl('10^{-4}'),textoidl('10^{-3}'),textoidl('10^{-2}'),textoidl('10^{-1}'),textoidl('1')] $
    ELSE plot,gamu.dens*hl*(2.0*H2 + HI),2.0*H2/(2.0*H2 + HI),psym = 3,/ylog,/xlog,xrange = [1e19,1e23],yrange = [1e-6,1],title = title[i],xtickname = [textoidl(' '),textoidl('10^{19}'),textoidl('10^{21}'),textoidl('10^{22}'),textoidl('10^{23}')]
    multiplot;,/doxaxis
ENDFOR

multiplot,/reset
IF KEYWORD_SET(outplot) THEN BEGIN
    device,/close
    device,filename = outplot+'metabund.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
    !p.multi = 0
ENDIF ELSE BEGIN
    stop
    window,1,xsize = 712,ysize = 392
    !p.multi = 0
ENDELSE

;!p.multi = [0,n,1]
multiplot,[n,1],mxtitle = textoidl('\rho [amu/cc]'),mxTitSize = 1.5,mxTitOffset = 2
FOR i = 0, n - 1 DO BEGIN
    rtipsy,files[i],h,g,d,s
    readarr,files[i]+'.H2',h,H2,/ascii,part = 'gas'
    readarr,files[i]+'.HI',h,HI,/ascii,part = 'gas'
;    readarr,files[i]+'.smoothlength',h,hl,/ascii,part = 'gas'

    gamu = g
    gamu.dens = g.dens*dens_convert/h.time/h.time/h.time
;    hl = hl*cm_per_kpc*kpc_per_syslength*h.time

    IF i eq 0 THEN plot,gamu.dens,2.0*H2/(2.0*H2 + HI),psym = 3,/xlog,xrange = [1e-1,1e4],yrange = [1e-6,1],title = title[i],ytitle = textoidl('H_2/H'),xtickname = [textoidl('10^{-1}'),textoidl('1'),textoidl('10^{1}'),textoidl('10^{2}'),textoidl('10^{3}'),textoidl('10^{4}')],ytickname = ['0.2','0.4','0.6','0.8','1.0'] $;,/ylog $
    ELSE plot,gamu.dens,2.0*H2/(2.0*H2 + HI),psym = 3,/xlog,xrange = [1e-1,1e4],yrange = [1e-6,1],title = title[i],xtickname = [textoidl(' '),textoidl('1'),textoidl('10^{1}'),textoidl('10^{2}'),textoidl('10^{3}'),textoidl('10^{4}')];,/ylog
    oplot,[100,100],[0,1]
    print,TOTAL(2.0*TOTAL(H2)/(2.0*TOTAL(H2) + TOTAL(HI)))
    multiplot
ENDFOR
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
multiplot,/reset

end


pro checkAbundMet_Obs,files,msol_per_sysmass,kpc_per_syslength,outplot = outplot,title = title,color = color
n = N_ELEMENTS(files)
spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '~/Code/Datafiles/MolecH/' ELSE prefix = '/astro/users/christensen/code/MolecH/'

formatplot,outplot = outplot
IF KEYWORD_SET(outplot) THEN BEGIN
    fgcolor = 0
    l_charsize = 0.75
    set_plot,'ps'
    device,filename = outplot+'metabundSDobs.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
ENDIF ELSE BEGIN
    l_charsize = 1.0
    fgcolor = 255
    set_plot,'x'
    window,0,xsize = 712,ysize = 392
ENDELSE
IF NOT KEYWORD_SET(title) THEN title = strind(n)
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
dens_convert =  msol_per_sysmass*gm_per_msol*amu_per_gm/kpc_per_syslength^3/cm_per_kpc^3 ;massunit*amu_per_gm/molec_weight/lengthunit^3 ;Converts to grams/cm^3
convert2msol = 1/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18

IF keyword_set(color) THEN BEGIN
    loadct,39
    obscolor = 240
    pointcolor = fgcolor
    datacolor = 50
ENDIF ELSE BEGIN
    loadct,0
    obscolor = fgcolor
    pointcolor = 100
    datacolor = fgcolor
ENDELSE


;************** Obs Data**********************
readcol,prefix + 'Wolfire08.dat',starw,namew,nhw,nh_erw,nh2w,nh2_erw,ncIw,ncI_erw,ncIIw,ncII_erw,avw,logfH2w,logfcIw,refw,format='A10,A8,F,F,F,F,F,F,F,F,F,F,F,I'
readcol,prefix + 'Gillmon06.dat',nameg,nh2g,nhg,refg,logfH2g,T01,Texc,format = 'A11,F,F,I,F,I,I,I'
multiplot,[n,1],mxtitle = textoidl('log \Sigma_{HI + H_{2}} [amu/cm^2]'),mxTitSize = 1.5,mxTitOffset = 2;,/doxaxis,/doyaxis
FOR i = 0, n - 1 DO BEGIN
;------------------- Tipsy Data
    rtipsy,files[i],h,g,d,s
    readarr,files[i]+'.H2',h,H2,/ascii,part = 'gas'
    readarr,files[i]+'.HI',h,HI,/ascii,part = 'gas'
    readarr,files[i]+'.smoothlength',h,hl,/ascii,part = 'gas'
    gamu = g
    gamu.dens = g.dens*dens_convert/h.time/h.time/h.time
    hl = hl*cm_per_kpc*kpc_per_syslength*h.time

;-------------------- Velocity Cube Data
    cubeHI = read_dencube_fits(files[i]+'.HI.fo.fits',headerHI)
    pix_areaHI = headerHI.CDELT1*headerHI.CDELT2*3.08568021d21*3.08568021d21
    cubeHI_surface_den = TRANSPOSE(cubeHI)/pix_areaHI*gm_per_msol*amu_per_gm
    vcubeHI = read_dencube_fits(files[i]+'.cubeHI.fits')/pix_areaHI*gm_per_msol*amu_per_gm/2.0

    cubeH2 = read_dencube_fits(files[i]+'.H2.fo.fits',headerH2)
    pix_areaH2 = headerH2.CDELT1*headerH2.CDELT2*3.08568021d21*3.08568021d21
    cubeH2_surface_den = TRANSPOSE(cubeH2)/pix_areaH2*gm_per_msol*amu_per_gm*2.0
    vcubeH2 = read_dencube_fits(files[i]+'.cubeH2.fits')/pix_areaHI*gm_per_msol*amu_per_gm

;Axes for the HI/H2 surface density
    xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1 + headerHI.CDELT1/2
    yaxes = (findgen(headerHI.NAXIS2) - headerHI.NAXIS1/2.0)*headerHI.CDELT2 + headerHI.CDELT1/2
    dxcell = headerHI.cdelt1
    areapc = 1e6*dxcell*dxcell
    areakpc = dxcell*dxcell
    minrange = headerHI.CRVAL1

    indplot = where(vcubeH2 gt 1e19)
    vcubeH = vcubeHI + 2.0*vcubeH2
    min1 = 19
    max1 = 23
    min2 = -6
    max2 = 0.4
    bin1 = 0.1                  ;(max1 - min1)/4.0/10.0
    bin2 = 0.1                  ;(max2 - min2)/6.0/10.0
    nxbins=FLOOR((max1-min1) / bin1)
    nybins=FLOOR((max2-min2) / bin2)
    xbins=(FINDGEN(nxbins)*bin1)+min1+bin1/2.0
    ybins=(FINDGEN(nybins)*bin2)+min2+bin2/2.0
    histv=HIST_2D(alog10(vcubeH),alog10(2.0*vcubeH2/(vcubeH)),min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)
;    histd=HIST_2D(alog10((cubeH2_surface_den + cubeHI_surface_den)*cos(angle)),alog10(cubeH2_surface_den/(cubeH2_surface_den + cubeHI_surface_den)),min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)
    indcool = where(g.tempg lt 1e4)
    histp=HIST_2D(alog10(gamu[indcool].dens*hl[indcool]*(2.0*H2[indcool] + HI[indcool])),alog10(2.0*H2[indcool]/(2.0*H2[indcool] + HI[indcool])),min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)

    IF (i eq 0) THEN BEGIN
        contour,histv,xbins,ybins,nlevels = 4 ,ytitle = textoidl('log H_2/H '),xrange = [19,23],yrange = [-6,0],xtickname = [textoidl('10^{19}'),textoidl('10^{20}'),textoidl('10^{21}'),textoidl('10^{22}'),textoidl('10^{23}')],ytickname = [textoidl('10^{-6}'),textoidl('10^{-5}'),textoidl('10^{-4}'),textoidl('10^{-3}'),textoidl('10^{-2}'),textoidl('10^{-1}'),textoidl('1')],title = title[i];,xtitle = textoidl('log(\Sigma_{HI + H_{2}}) [amu/cm^2]')
        contour,histp,xbins,ybins,nlevels = 20,color = pointcolor,levels = [5,25,45,65],/overplot
        contour,histv,xbins,ybins,nlevels = 4,/overplot,color = obscolor
        oplot,nhw,logfH2w,psym = 7,color = datacolor
        oplot,nhg,logfH2g,psym = 1,color = datacolor
;    plot,gamu.dens*hl*(2.0*H2 + HI),2.0*H2/(2.0*H2 + HI),psym = 3,/ylog,/xlog,xrange = [1e19,1e23],yrange = [1e-6,1],title = title[i],ytitle = textoidl('H_2/H'),xtickname = [textoidl('10^{19}'),textoidl('10^{20}'),textoidl('10^{21}'),textoidl('10^{22}'),textoidl('10^{23}')],ytickname = [textoidl('10^{-6}'),textoidl('10^{-5}'),textoidl('10^{-4}'),textoidl('10^{-3}'),textoidl('10^{-2}'),textoidl('10^{-1}'),textoidl('1')] $
      ENDIF ELSE BEGIN 
        contour,histv,xbins,ybins,nlevels = 4 ,xrange = [19,23],yrange = [-6,0],xtickname = [' ',textoidl('10^{20}'),textoidl('10^{21}'),textoidl('10^{22}'),textoidl('10^{23}')]
        contour,histp,xbins,ybins,nlevels = 4 ,xrange = [19,23],yrange = [-6,0],xtickname = [' ',textoidl('10^{20}'),textoidl('10^{21}'),textoidl('10^{22}'),textoidl('10^{23}')],min = 2,levels = [5,25,45,65],title = title[i];,xtitle = textoidl('log(\Sigma_{HI + H_{2}}) [amu/cm^2]')
;        contour,histp,xbins,ybins,nlevels = 4 ,xrange = [19,23],yrange = [-6,0],min = 2,levels = [5,25,45,65],color = obscolor
        contour,histp,xbins,ybins,nlevels = 4,color = pointcolor,/overplot,levels = [5,25,45,65]
        contour,histv,xbins,ybins,nlevels = 4 ,/overplot
;    plot,gamu.dens*hl*(2.0*H2 + HI),2.0*H2/(2.0*H2 + HI),psym = 3,/ylog,/xlog,xrange = [1e19,1e23],yrange = [1e-6,1],title = title[i],xtickname = [textoidl(' '),textoidl('10^{19}'),textoidl('10^{21}'),textoidl('10^{22}'),textoidl('10^{23}')]
    ENDELSE
    multiplot;,/doxaxis
ENDFOR

IF KEYWORD_SET(outplot) THEN BEGIN
    device,/close
    set_plot,'ps'
    device,filename = outplot+'metabundSDobs2.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
ENDIF  ELSE stop
multiplot,/reset
multiplot,[n,1],mxtitle = textoidl('log \Sigma_{HI + H_{2}} [amu/cm^2]'),mxTitSize = 1.5,mxTitOffset = 2;,/doxaxis,/doyaxis

;-----------------------------------------------------------
;In solar masses
FOR i = 0, n - 1 DO BEGIN
;------------------- Tipsy Data
    rtipsy,files[i],h,g,d,s
    readarr,files[i]+'.H2',h,H2,/ascii,part = 'gas'
    readarr,files[i]+'.HI',h,HI,/ascii,part = 'gas'
    readarr,files[i]+'.smoothlength',h,hl,/ascii,part = 'gas'
    gamu = g
    gamu.dens = g.dens*dens_convert/h.time/h.time/h.time
    hl = hl*cm_per_kpc*kpc_per_syslength*h.time

;-------------------- Velocity Cube Data
    cubeHI = read_dencube_fits(files[i]+'.HI.fo.fits',headerHI)
    pix_areaHI = headerHI.CDELT1*headerHI.CDELT2*3.08568021d21*3.08568021d21
    cubeHI_surface_den = TRANSPOSE(cubeHI)/pix_areaHI*gm_per_msol*amu_per_gm
    vcubeHI = read_dencube_fits(files[i]+'.cubeHI.fits')/pix_areaHI*gm_per_msol*amu_per_gm/2.0

    cubeH2 = read_dencube_fits(files[i]+'.H2.fo.fits',headerH2)
    pix_areaH2 = headerH2.CDELT1*headerH2.CDELT2*3.08568021d21*3.08568021d21
    cubeH2_surface_den = TRANSPOSE(cubeH2)/pix_areaH2*gm_per_msol*amu_per_gm*2.0
    vcubeH2 = read_dencube_fits(files[i]+'.cubeH2.fits')/pix_areaHI*gm_per_msol*amu_per_gm

;Axes for the HI/H2 surface density
    xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1 + headerHI.CDELT1/2
    yaxes = (findgen(headerHI.NAXIS2) - headerHI.NAXIS1/2.0)*headerHI.CDELT2 + headerHI.CDELT1/2
    dxcell = headerHI.cdelt1
    areapc = 1e6*dxcell*dxcell
    areakpc = dxcell*dxcell
    minrange = headerHI.CRVAL1

    indplot = where(vcubeH2 gt 1e19)
    vcubeH = vcubeHI + 2.0*vcubeH2
    min1 = -1
    max1 = 3
    min2 = -6
    max2 = 0.4
    bin1 = 0.1                  ;(max1 - min1)/4.0/10.0
    bin2 = 0.1                  ;(max2 - min2)/6.0/10.0
    nxbins=FLOOR((max1-min1) / bin1)
    nybins=FLOOR((max2-min2) / bin2)
    xbins=(FINDGEN(nxbins)*bin1)+min1+bin1/2.0
    ybins=(FINDGEN(nybins)*bin2)+min2+bin2/2.0
    histv=HIST_2D(alog10(vcubeH*convert2msol),alog10(2.0*vcubeH2/(vcubeH)),min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)
;    histd=HIST_2D(alog10((cubeH2_surface_den + cubeHI_surface_den)*cos(angle)),alog10(cubeH2_surface_den/(cubeH2_surface_den + cubeHI_surface_den)),min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)
    indcool = where(g.tempg lt 1e4)
    histp=HIST_2D(alog10(convert2msol*gamu[indcool].dens*hl[indcool]*(2.0*H2[indcool] + HI[indcool])),alog10(2.0*H2[indcool]/(2.0*H2[indcool] + HI[indcool])),min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)

    IF (i eq 0) THEN BEGIN
        contour,histv,xbins,ybins,nlevels = 4 ,ytitle = textoidl('log H_2/H '),xrange = [-1,3],yrange = [-6,0],xtickname = [textoidl('10^{-1}'),textoidl('1'),textoidl('10'),textoidl('10^{2}'),textoidl('10^{3}')],ytickname = [textoidl('10^{-6}'),textoidl('10^{-5}'),textoidl('10^{-4}'),textoidl('10^{-3}'),textoidl('10^{-2}'),textoidl('10^{-1}'),textoidl('1')],title = title[i],/nodata;,xtitle = textoidl('log(\Sigma_{HI + H_{2}}) [amu/cm^2]')
        contour,histp,xbins,ybins,nlevels = 20,color = pointcolor,levels = [5,25,45,65],/overplot
        contour,histv,xbins,ybins,nlevels = 4,/overplot,color = 30,thick = 3
        oplot,alog10(convert2msol*10^nhw),logfH2w,psym = 7,color = datacolor
        oplot,alog10(convert2msol*10^nhg),logfH2g,psym = 1,color = datacolor
      ENDIF ELSE BEGIN 
        contour,histp,xbins,ybins,nlevels = 4 ,xrange = [-1,3],yrange = [-6,0],xtickname = [' ',textoidl('1'),textoidl('10'),textoidl('10^{2}'),textoidl('10^{3}')],min = 2,levels = [5,25,45,65],title = title[i];,xtitle = textoidl('log(\Sigma_{HI + H_{2}}) [amu/cm^2]')
        contour,histp,xbins,ybins,nlevels = 4,color = pointcolor,/overplot,levels = [5,25,45,65]
    ENDELSE
    multiplot;,/doxaxis
ENDFOR

IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
multiplot,/reset

end
