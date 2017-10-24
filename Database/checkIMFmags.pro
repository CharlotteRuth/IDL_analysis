pro checkIMFmags,outplot = outplot
loadct,39
IF KEYWORD_SET(outplot) THEN BEGIN
   set_plot, 'ps'
   !P.THICK=1.5                 ;4
   !P.CHARTHICK=1.5             ;4
   !X.THICK=1.5                 ;4
   !Y.THICK=1.5                 ;4
   !p.charsize=1.2              ;1.8
   !p.font=0
ENDIF ELSE set_plot,'x'

cd,'/astro/net/nbody1/christensen/Sunrise/Get_halo_mags_Test/h258.cosmo50cmb.512g2bwK.00512/'
kroupa = MRDFITS('h258.cosmo50cmb.512g2bwK.00512.amiga_vir.halos.star99_MS_ab.Mv.fits',1)
ms = MRDFITS('h258.cosmo50cmb.512g2bwK.00512.amiga_vir.halos.star99_K_ab.Mv.fits',1)
FMT = 'X,X,X,X,X,F,X,F,F,F,X,X,X,X,X,X,X,X,X,X,A'
readcol,'/astro/net/nbody1/christensen/Sunrise/Get_halo_mags_Test/h258.cosmo50cmb.512g2bwK.00512/h258.cosmo50cmb.512g2bwK.00512.amiga.stat',F=FMT,vir,gas,star,dark,contam,/silent 
hstar = where(kroupa.u ne 0)
vir = vir[hstar]
star = star[hstar]
kroupa = kroupa[hstar]
ms = ms[hstar]

IF KEYWORD_SET(outplot) THEN device,filename='/astro/users/christensen/code/Database/MS_K_mags.ps',/COLOR,bits_per_pixel= 8,/times,ysize=3.5,xsize=5,/inch
plot,vir,kroupa.u,psym = 4,xtitle = 'Halo Mass [Msol]',ytitle = 'AB Magnitude',yrange = [-10,-25],title = 'h258.cosmo50cmb.512g2bwK.00512',/xlog
oplot,vir,kroupa.u,psym = 4,color = 30
oplot,vir,ms.u,psym = 5,color = 30

oplot,vir,kroupa.b,psym = 4,color = 70
oplot,vir,ms.b,psym = 5,color = 70

oplot,vir,kroupa.v,psym = 4,color = 120
oplot,vir,ms.v,psym = 5,color = 120

oplot,vir,kroupa.r,psym = 4,color = 180
oplot,vir,ms.r,psym = 5,color = 180

oplot,vir,kroupa.i,psym = 4,color = 220
oplot,vir,ms.i,psym = 5,color = 220

oplot,vir,kroupa.k,psym = 4,color = 240
oplot,vir,ms.k,psym = 5,color = 240
legend,['Kroupa IMF','Miller-Scalo IMF'],psym = [4,5]
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop

IF KEYWORD_SET(outplot) THEN device,filename='/astro/users/christensen/code/Database/MS_K_mags_diff.ps',/COLOR,bits_per_pixel= 8,/times,ysize=3.5,xsize=5,/inch
plot,[1e9,1e12],[0,0],xtitle = 'Halo Mass [Msol]',ytitle = 'Kroupa - Miller-Scalo [AB Mag]',yrange = [-0.6, 0.6],title = 'h258.cosmo50cmb.512g2bwK.00512',/xlog
oplot,vir,kroupa.u - ms.u,psym = 4,color = 30
oplot,vir,kroupa.b - ms.b,psym = 4,color = 70
oplot,vir,kroupa.v - ms.v,psym = 4,color = 120
oplot,vir,kroupa.r - ms.r,psym = 4,color = 180
oplot,vir,kroupa.i - ms.i,psym = 4,color = 220
oplot,vir,kroupa.k - ms.k,psym = 4,color = 240
legend,['U','B','V','R','I','K'],psym = [4,4,4,4,4,4],color = [30,70,120,180,220,240],/right,/top
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop

mags_bes_k = transpose([[kroupa.u],[kroupa.b],[kroupa.v],[kroupa.r],[kroupa.i]])
mags_bes_ms = transpose([[ms.u],[ms.b],[ms.v],[ms.r],[ms.i]])
mags_errs = fltarr(5,N_ELEMENTS(kroupa.u))+0.02
;mags_errs = [0.05, 0.02, 0.02, 0.02, 0.03] ;From the kcorrect website ;ABS(sdss_errs[*,i] - sdss_mags[*,i])
mags_ivar=1./mags_errs^2
dmod = 19.4576 ;distance modulus for z = 0 in this model
z = fltarr(N_ELEMENTS(kroupa.u))

mgy_k =(10.D)^(-(0.4D)*(mags_bes_k + dmod))
mgy_ivar_k = mags_ivar/(0.4*alog(10.)*mgy_k)^2.
kcorrect, mgy_k, mgy_ivar_k, z, kcorrect, mass = mass_k, mtol = mtol_k, absmag = absmag_k, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']

mgy_ms =(10.D)^(-(0.4D)*(mags_bes_ms + dmod))
mgy_ivar_ms = mags_ivar/(0.4*alog(10.)*mgy_ms)^2.
kcorrect, mgy_ms, mgy_ivar_ms, z, kcorrect, mass = mass_ms, mtol = mtol_ms, absmag = absmag_ms, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']

IF KEYWORD_SET(outplot) THEN device,filename='/astro/users/christensen/code/Database/MS_K_mass.ps',/COLOR,bits_per_pixel= 8,/times,ysize=3.5,xsize=5,/inch
plot,vir,mass_ms,psym = 4,/ylog,xtitle = 'Halo Mass [Msol]',ytitle = 'Stellar Mass [Msol]',yrange = [1e6,2e11],ystyle = 1,title = 'h258.cosmo50cmb.512g2bwK.00512',/xlog
x = findgen(700)/100 + 9
y = alog10(0.1257*((10^x/10^(11.36))^(-0.9147) + (10^x/10^(11.36))^(0.2485))^(-2.574)*10^x) ;guo 09, (3)
oplot,10^x,10^y
oplot,vir,mass_ms,psym = 4,color = 50
oplot,vir,mass_k,psym = 5,color = 150
;oplot,vir,star,psym = 6,color = 240
;legend,['Kcorrect from MS','Kcorrect from Kroupa','Amiga Stellar Mass'],color = [50,150,240],psym = [4,5,6],/left,/top
legend,['Kcorrect from MS','Kcorrect from Kroupa'],color = [50,150],psym = [4,5],/left,/top
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
END
