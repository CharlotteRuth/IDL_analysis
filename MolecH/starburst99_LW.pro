PRO starburst99_LW
formatplot
;see envelope_LW for the back of the envelope calculation
h = 6.626068e-27; ergs seconds
c = 3e18 ;A/s
h_ev = 4.1356668e-15 ;ev seconds
zsolar =  0.0130215
loadct,39

;LW radiation is between 12.24 eV and 13.51 eV
; i.e, 1013.64 A and 918.357 A
lambda1 = 1013.64
lambda2 = 918.357
lambda = 995.
metallicity = zsolar
E = h*c/lambda ;in ergs
time = 1e7
thick = 3
;--------------- Back of the Envelope Calculation ---------------------------
fit = [ -332.61118, 210.37616, -43.684908, 4.0379845, -0.14145357]
age_fit = 10^(findgen(100)*(14 - 6)/100.0 + 6.)
age_rad_fit = fit[0] + $
              fit[1]*alog10(age_fit) + $
              fit[2]*alog10(age_fit)*alog10(age_fit) + $ 
              fit[3]*alog10(age_fit)*alog10(age_fit)*alog10(age_fit) + $
              fit[4]*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)
window,0
plot,age_fit,age_rad_fit,/xlog,xrange = [1e6,1e11],yrange = [0,47],xtitle = 'Age [yr]',ytitle = 'LOG_10(I_LW per Stellar Mass) [photons s^-1 M_sol^-1]',thick = thick
;oplot,age_fit,age_rad_fit - alog10(2.0),linestyle = 1
oplot,[1e7,1e7],[0,50]
legend,['Back of the Envelope Fit'],linestyle = 0,/left,/bottom
;stop

;-------------- My Starburst99 Simulation ------------------------------------
IF 0 THEN BEGIN
readcol,'/astro/users/christensen/Gasoline_H2/SB99/output/SimpleLW.spectrum1',ageSB,waveSB,totSB,sspSB,nebulaSB
; units of year, angstrom, ergs/s/A/1e6 solar masses
mass = 1e6 ;solar masses
totSB = totSB - 6
sspSB = sspSB - 6
nebulaSB = nebulaSB - 6

indwSB = where(waveSB eq lambda)
temp = MIN(ABS(ageSB - time),indtSB_age)
indtSB = where(ageSB eq (ageSB[indtSB_age])[0])
totSB_phot = totSB - alog10(E) + alog10(lambda1 - lambda2) ;Total Emission
nebSB_phot = nebulaSB - alog10(E) + alog10(lambda1 - lambda2); Nebular Emission
sspSB_phot = sspSB - alog10(E) + alog10(lambda1 - lambda2) ;Simple Stellar Population Emission

oplot,ageSB[indwSB],nebSB_phot[indwSB],linestyle = 1, color = 210,thick = thick ;Nebular
oplot,ageSB[indwSB],sspSB_phot[indwSB],linestyle = 3, color = 210,thick = thick ;SSP
oplot,ageSB[indwSB],totSB_phot[indwSB],linestyle = 2, color = 240,thick = thick ;Total
legend,['SSP','Nebula','Total'],linestyle = [3,1,2],color = [210,210,240],/right,/bottom
;stop
ENDIF

;-------------- Sunrise SSP -----------------------------------------------------
dir = '/home/christensen/Storage1/UW/' ;'/astro/net/safe1/christensen/'
SSPaxes = mrdfits(dir + 'Sunrise/datafiles_v3/Patrik-imfKroupa-Zmulti-ml-subsampled.fits',4)
SSPdata = mrdfits(dir + 'Sunrise/datafiles_v3/Patrik-imfKroupa-Zmulti-ml-subsampled.fits',2) ;Watts per Meter, I think it is in log
SSPdata = SSPdata - 3 ;ergs / s / Angstrom (in log)
waveSR = SSPaxes.lambda*1e10 ;angstrom
temp = min(abs(waveSR - lambda),indwSR)
temp = min(abs(SSPaxes.metallicity - metallicity),indzSR)
temp = min(abs(SSPaxes.time - time),indtSR)
SSPLW_phot = SSPdata[*,indzSR,indwSR] - alog10(E) + alog10(lambda1 - lambda2)
indyoung = where(SSPaxes.TIME lt 1e7 AND SSPaxes.TIME ne 0)
indmed = where( SSPaxes.TIME gt 1e7 AND SSPaxes.TIME lt 1e9)
indold = where( SSPaxes.TIME gt 1e9 AND SSPaxes.TIME lt 1e10)
oplot,SSPaxes.TIME,SSPLW_phot,linestyle = 2, color = 80,thick = thick ;All ages
;oplot,SSPaxes[indyoung].TIME,SSPLW_phot[indyoung] - alog10(5.0),linestyle = 2, color = 50,thick = thick,psym = 2 ;Young
oplot,SSPaxes[indmed].TIME,SSPLW_phot[indmed],linestyle = 3, color = 50,thick = thick ;Medium
oplot,SSPaxes[indyoung].TIME,SSPaxes[indyoung].TIME*0 + SSPLW_phot[indmed[0]],linestyle = 3, color = 50,thick = thick ;Old
fitSR = POLY_FIT(alog10(SSPaxes[indmed].TIME),SSPLW_phot[indmed],5)
age_LWSR_fit = fitSR[0] + $
              fitSR[1]*alog10(age_fit) + $
              fitSR[2]*alog10(age_fit)*alog10(age_fit) + $ 
              fitSR[3]*alog10(age_fit)*alog10(age_fit)*alog10(age_fit) + $
              fitSR[4]*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)+ $
              fitSR[5]*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)
;oplot,age_fit,age_LWSR_fit,linestyle = 0,color = 50,thick = thick
;legend,['Sunrise SSP','Fit to Sunrise SSP'],color = [50,50],linestyle = [2,0],/right,/top
legend,['Sunrise SSP'],color = [50],linestyle = [2],/right,/top


fitSRold = POLY_FIT(alog10(SSPaxes[indold].TIME),SSPLW_phot[indold],2)
age_LWSRold_fit = fitSRold[0] + $
              fitSRold[1]*alog10(age_fit); + $
;              fitSR[2]*alog10(age_fit)*alog10(age_fit) + $ 
;              fitSR[3]*alog10(age_fit)*alog10(age_fit)*alog10(age_fit) + $
;              fitSR[4]*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)+ $
;              fitSR[5]*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)
oplot,age_fit,age_LWSRold_fit,linestyle = 0,color = 100,thick = thick

stop

;-------------- Sunrise Mappings ------------------------------------------------
PDRaxes = mrdfits(dir + 'Sunrise/datafiles_v3/Smodel.fits',1)
PDRdata = mrdfits(dir + 'Sunrise/datafiles_v3/Smodel.fits',2)
wavePDR = pdraxes.wave*1e5 ;angstrom
ind = where(wavePDR gt lambda2 AND wavePDR lt lambda1)
temp = min(abs(wavePDR - lambda),indwPDR)
temp = min(abs(PDRaxes.metal*ZSOLAR - metallicity),indzPDR)
PDRdata_select = PDRdata[indzPDR,*,*,*,indwPDR]

stop
;--------------- Spectrum
window,1
plot,waveSR,SSPdata[indtSR,indzSR,*],/xlog,xtitle = 'Wavelength [A]',ytitle = 'LOG(Energy) [ergs/s/A]',yrange = [1,40],xrange = [100,1e6]
oplot,[lambda,lambda],[0,100]
oplot,waveSR,SSPdata[indtSR,indzSR,*],color = 50
IF 0 THEN oplot,waveSB[indtSB],totSB[indtSB],color = 240
stop

loadct,0
formatplot,/outplot
device, filename='~/plots/LWfit.eps',/COLOR,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch
plot,alog10(SSPaxes.TIME),SSPLW_phot,xtitle = 'Log Age [yr]',ytitle = textoidl('I_{LW} / Stellar Mass [photons s^{-1} M')+sunsymbol()+textoidl('^{-1}]'),/nodata,xrange = [6,10],yrange = [30,47],xstyle = 1;,xtickv  = [6,7,8,9,10]
oplot,alog10(SSPaxes.TIME),SSPLW_phot,color = 120,thick = 6
fitregion = where(age_fit ge 1e7 AND age_fit le 1e9)
lowregion = where(age_fit lt 1e7)
highregion = where(age_fit gt 1e9)
age_LWSR_total = age_LWSR_fit
age_LWSR_total[fitregion] = age_LWSR_fit[fitregion]
age_LWSR_total[lowregion] = 0*age_LWSR_total[lowregion] + age_LWSR_fit[fitregion[0]]
age_LWSR_total[highregion] = 0
oplot,alog10(age_fit),age_LWSR_total,linestyle = 2,thick = 6
oplot,[7,7],[30,50],linestyle = 1,thick = 2
legend,['Starburst99 SSP','Fit'],color = [120,0],linestyle = [0,2],thick = 6,/right
device,/close
END



