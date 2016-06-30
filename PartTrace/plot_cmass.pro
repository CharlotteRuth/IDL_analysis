PRO plot_cmass,dirs,files,finalids,step1,step2,outplot = outplot
formatplot,outplot = outplot,thick = formatthick
;loadct,39
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 18
    thick = 2
    symsize = 2.5
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 400
    ysize = 400
    thick = 1
    symsize = 1.5
ENDELSE
device,decomposed = 0

colors = [30,120,60,254,100,220,160,80]
colors = [30,120,60,60,254,100,230,160,80,210]

yrange = [5e7,2e10]
xrange = [4e9,1e12]

IF keyword_set(outplot) THEN  device,filename = outplot + '_cm500.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,[0,0],[0,0],xtitle = textoidl('M_{vir} [M') + sunsymbol() + ']',ytitle = textoidl('M_{500} [M') + sunsymbol() + ']',psym = 4,/ylog,/xlog,yrange = yrange,xrange = xrange,/nodata
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   halodat = mrdfits(dirs[i] + '/grp' + finalids[i] + '.alignment.fits',1)
   readcol,dirs[i] + '/grp' + finalids[i] + '.mass_metal_track.dat',halo_t,time_t,z_t,mtot_t,mgas_t,mstar_t,mdark_t,metal_t,ox_t,fe_t,HI_t,H2_t,coldg_t,/silent 
   readcol,dirs[i] + '/grp' + finalids[i] + '.mass500.txt',haloids,time,z,gmass,smass,dmass,slope500,gmass1,smass1,dmass1,slope1,/silent  
   tmass = gmass + smass + dmass
   ind1 = where(halodat.file EQ files[i] + '.' + step1[i] + '/' + files[i] + '.' + step1[i])
   ind2 = where(halodat.file EQ files[i] + '.' + step2[i] + '/' + files[i] + '.' + step2[i])
   oplot,[mtot_t[ind1],mtot_t[ind2]],[tmass[ind1],tmass[ind2]],thick = thick
   oplot,[mtot_t[ind1],mtot_t[ind1]],[tmass[ind1],tmass[ind1]],psym = symcat(4),color = colors[i],symsize = symsize
   oplot,[mtot_t[ind2],mtot_t[ind2]],[tmass[ind2],tmass[ind2]],psym = symcat(14),color = colors[i],symsize = symsize
ENDFOR
IF keyword_set(outplot) THEN  device,/close ELSE stop

;--------------------------------------
yrange = [5e7,2e10]
IF keyword_set(outplot) THEN  device,filename = outplot + '_cm1.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,[0,0],[0,0],xtitle = textoidl('M_{vir} [M') + sunsymbol() + ']',ytitle = textoidl('M_{cen} [M') + sunsymbol() + ']',psym = 4,/ylog,/xlog,yrange = yrange,xrange = xrange,/nodata
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   halodat = mrdfits(dirs[i] + '/grp' + finalids[i] + '.alignment.fits',1)
   readcol,dirs[i] + '/grp' + finalids[i] + '.mass_metal_track.dat',halo_t,time_t,z_t,mtot_t,mgas_t,mstar_t,mdark_t,metal_t,ox_t,fe_t,HI_t,H2_t,coldg_t,/silent 
   readcol,dirs[i] + '/grp' + finalids[i] + '.mass500.txt',haloids,time,z,gmass,smass,dmass,slope500,gmass1,smass1,dmass1,slope1,/silent 
   tmass1 = gmass1 + smass1 + dmass1
   ind1 = where(halodat.file EQ files[i] + '.' + step1[i] + '/' + files[i] + '.' + step1[i])
   ind2 = where(halodat.file EQ files[i] + '.' + step2[i] + '/' + files[i] + '.' + step2[i])
   print,dirs[i] + '/grp' + finalids[i],z[ind1]
   oplot,[mtot_t[ind1],mtot_t[ind2]],[tmass1[ind1],tmass1[ind2]],thick = thick
   oplot,[mtot_t[ind1],mtot_t[ind1]],[tmass1[ind1],tmass1[ind1]],psym = symcat(4),color = colors[i],symsize = symsize
   oplot,[mtot_t[ind2],mtot_t[ind2]],[tmass1[ind2],tmass1[ind2]],psym = symcat(14),color = colors[i],symsize = symsize
ENDFOR
oplot,[1.4e11,1.4e11],[3.3e9,3.3e9],psym = symcat(45),symsize = symsize
oplot,[3.1e11,3.1e11],[2.7e9,2.7e9],psym = symcat(46),symsize = symsize
oplot,[1.4e11,3.1e11],[3.3e9,2.7e9],thick = thick
IF keyword_set(outplot) THEN  device,/close ELSE stop

yrange = [-7,0]
;----------------------------------------
IF keyword_set(outplot) THEN  device,filename = outplot + '_s500.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,[0,0],[0,0],xtitle = textoidl('M_{vir} [M') + sunsymbol() + ']',ytitle = textoidl('\alpha_{500} [M') + sunsymbol() + ']',psym = 4,/xlog,yrange = yrange,xrange = xrange,/nodata
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   halodat = mrdfits(dirs[i] + '/grp' + finalids[i] + '.alignment.fits',1)
   readcol,dirs[i] + '/grp' + finalids[i] + '.mass_metal_track.dat',halo_t,time_t,z_t,mtot_t,mgas_t,mstar_t,mdark_t,metal_t,ox_t,fe_t,HI_t,H2_t,coldg_t,/silent 
   readcol,dirs[i] + '/grp' + finalids[i] + '.mass500.txt',haloids,time,z,gmass,smass,dmass,slope500,gmass1,smass1,dmass1,slope1,/silent 
   tmass = gmass + smass + dmass
   ind1 = where(halodat.file EQ files[i] + '.' + step1[i] + '/' + files[i] + '.' + step1[i])
   ind2 = where(halodat.file EQ files[i] + '.' + step2[i] + '/' + files[i] + '.' + step2[i])
   oplot,[mtot_t[ind1],mtot_t[ind2]],[slope500[ind1],slope500[ind2]],thick = thick
   oplot,[mtot_t[ind1],mtot_t[ind1]],[slope500[ind1],slope500[ind1]],psym = symcat(4),color = colors[i],symsize = symsize
   oplot,[mtot_t[ind2],mtot_t[ind2]],[slope500[ind2],slope500[ind2]],psym = symcat(14),color = colors[i],symsize = symsize
ENDFOR
IF keyword_set(outplot) THEN  device,/close ELSE stop

yrange = [-13,1]
IF keyword_set(outplot) THEN  device,filename = outplot + '_s1000.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,[0,0],[0,0],xtitle = textoidl('M_{vir} [M') + sunsymbol() + ']',ytitle = textoidl('\alpha_{1000} [M') + sunsymbol() + ']',psym = 4,/xlog,yrange = yrange,xrange = xrange,/nodata
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   halodat = mrdfits(dirs[i] + '/grp' + finalids[i] + '.alignment.fits',1)
   readcol,dirs[i] + '/grp' + finalids[i] + '.mass_metal_track.dat',halo_t,time_t,z_t,mtot_t,mgas_t,mstar_t,mdark_t,metal_t,ox_t,fe_t,HI_t,H2_t,coldg_t,/silent 
   readcol,dirs[i] + '/grp' + finalids[i] + '.mass500.txt',haloids,time,z,gmass,smass,dmass,slope500,gmass1,smass1,dmass1,slope1,/silent 
   tmass = gmass + smass + dmass
   ind1 = where(halodat.file EQ files[i] + '.' + step1[i] + '/' + files[i] + '.' + step1[i])
   ind2 = where(halodat.file EQ files[i] + '.' + step2[i] + '/' + files[i] + '.' + step2[i])
   oplot,[mtot_t[ind1],mtot_t[ind2]],[slope1[ind1],slope1[ind2]],thick = thick
   oplot,[mtot_t[ind1],mtot_t[ind1]],[slope1[ind1],slope1[ind1]],psym = symcat(4),color = colors[i],symsize = symsize
   oplot,[mtot_t[ind2],mtot_t[ind2]],[slope1[ind2],slope1[ind2]],psym = symcat(14),color = colors[i],symsize = symsize
   print,[mtot_t[ind1],mtot_t[ind2]],[slope1[ind1],slope1[ind2]]
ENDFOR
IF keyword_set(outplot) THEN  device,/close ELSE stop

END
