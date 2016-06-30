pro vpeak_plot,outplot = outplot
filebase = 'cosmo50cmb.256g2MbwK.00512'
cd,'/astro/store/nbody2/christensen/HaloSelect'

readcol,filebase+'_halodat_red.dat',group_r,mass_r,smass_r,smassfrac_r,gmassfrac_r,vpeak_r,xc_r,yc_r,zc_r,format='(I,F,F,F,F,F,F,F,F)'
readcol,filebase+'_halodat_blue.dat',group_b,mass_b,smass_b,smassfrac_b,gmassfrac_b,vpeak_b,xc_b,yc_b,zc_b,format='(I,F,F,F,F,F,F,F,F)'
ind_r = where(smass_r ge 1e10 AND smass_r le 5e10)
ind_b = where(smass_b ge 1e10 AND smass_b le 5e10)

loadct,39
formatplot,outplot = outplot
IF KEYWORD_SET(outplot) THEN device,filename = 'smass_vpeak.eps',/encapsulated,/color,/times,ysize=xsize,xsize=xsize,bits_per_pixel= 8
plot,smass_b[ind_b],vpeak_b[ind_b],xtitle = 'Stellar Mass [M'+sunsymbol()+']',ytitle = 'Peak Velocity [km/s]',/nodata,/xlog,xrange = [4e9,1e11],yrange = [100,350]
oplot,smass_b,vpeak_b,psym = symcat(6),color = 80
oplot,smass_r,vpeak_r,psym = symcat(9),color = 220
oplot,smass_b[ind_b],vpeak_b[ind_b],psym = symcat(15),color = 60,symsize = 2
oplot,smass_r[ind_r],vpeak_r[ind_r],psym = symcat(16),color = 254,symsize = 2
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop

IF KEYWORD_SET(outplot) THEN device,filename = 'vmass_vpeak.eps',/encapsulated,/color,/times,ysize=xsize,xsize=xsize,bits_per_pixel= 8
plot,mass_b,vpeak_b,xtitle = 'Virial Mass [M'+sunsymbol()+']',ytitle = 'Peak Velocity [km/s]',/nodata,/xlog,xrange = [1e11,2e12],yrange = [100,350]
oplot,mass_b,vpeak_b,psym = symcat(6),color = 80
oplot,mass_r,vpeak_r,psym = symcat(9),color = 220
oplot,mass_b[ind_b],vpeak_b[ind_b],psym = symcat(15),color = 60,symsize = 2
oplot,mass_r[ind_r],vpeak_r[ind_r],psym = symcat(16),color = 254,symsize = 2
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop

IF KEYWORD_SET(outplot) THEN device,filename = 'bfrac_vpeak.eps',/encapsulated,/color,/times,ysize=xsize,xsize=xsize,bits_per_pixel= 8
plot,smass_b/mass_b,vpeak_b,xtitle = 'Stellar Mass /Virial Mass',ytitle = 'Peak Velocity [km/s]',/nodata,xrange = [0,0.11],yrange = [100,350];,/xlog;
oplot,smass_b/mass_b,vpeak_b,psym = symcat(6),color = 80
oplot,smass_r/mass_r,vpeak_r,psym = symcat(9),color = 220
oplot,smass_b[ind_b]/mass_b[ind_b],vpeak_b[ind_b],psym = symcat(15),color = 60,symsize = 2
oplot,smass_r[ind_r]/mass_r[ind_r],vpeak_r[ind_r],psym = symcat(16),color = 254,symsize = 2
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
end
