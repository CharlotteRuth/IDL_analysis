pro plot_baryonTF
msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 1.69875e16 ; Solar mass in system units
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 50000.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70
loadct,39
;set_plot,'x'

filename1 = 'X2X2randv.dat'
filename2 = 'X3X3randv.dat'
filename3 = 'X5X5randv.dat'

readcol,filename1,r11, v11, sm11, bm11, r21, v21, sm21, bm21, r31, v31, sm31, bm31,format = 'F,F,F,F,F,F,F,F,F,F,F,F'
;bm2001 = bm2001 - sm2001*dMsolUnit + sm2001
;bm11 = bm11 - sm11*dMsolUnit + sm11
;bm21 = bm21 - sm21*dMsolUnit + sm21
readcol,filename2,r12, v12, sm12, bm12, r22, v22, sm22, bm22, r32, v32, sm32, bm32,format = 'F,F,F,F,F,F,F,F,F,F,F,F'

;bm2002 = bm2002 - sm2002*dMsolUnit + sm2002
;bm12 = bm12 - sm12*dMsolUnit + sm12
;bm22 = bm22 - sm22*dMsolUnit + sm22
readcol,filename3,r13, v13, sm13, bm13, r23, v23, sm23, bm23, r33, v33, sm33, bm33,format = 'F,F,F,F,F,F,F,F,F,F,F,F'
;bm2003 = bm2003 - sm2003*dMsolUnit + sm2003
;bm13 = bm13 - sm13*dMsolUnit + sm13
;bm23 = bm23 - sm23*dMsolUnit + sm23

logv11 = ALOG10(v11)
logbm11 = ALOG10(bm11)
logv12 = ALOG10(v12)
logbm12 = ALOG10(bm12)
logv13 = ALOG10(v13)
logbm13 = ALOG10(bm13)
logv21 = ALOG10(v21)
logbm21 = ALOG10(bm21)
logv22 = ALOG10(v22)
logbm22 = ALOG10(bm22)
logv23 = ALOG10(v23)
logbm23 = ALOG10(bm23)
logv31 = ALOG10(v31)
logbm31 = ALOG10(bm31)
logv32 = ALOG10(v32)
logbm32 = ALOG10(bm32)
logv33 = ALOG10(v33)
logbm33 = ALOG10(bm33)



bestfitx11 = findgen(10)*(MAX(logv11) - MIN(logv11))/10.0 + MIN(logv11)
bestfitx12 = findgen(10)*(MAX(logv12) - MIN(logv12))/10.0 + MIN(logv12)
bestfitx13 = findgen(10)*(MAX(logv13) - MIN(logv13))/10.0 + MIN(logv13)

bestfitx21 = findgen(10)*(MAX(logv21) - MIN(logv21))/10.0 + MIN(logv21)
bestfitx22 = findgen(10)*(MAX(logv22) - MIN(logv22))/10.0 + MIN(logv22)
bestfitx23 = findgen(10)*(MAX(logv23) - MIN(logv23))/10.0 + MIN(logv23)

bestfitx31 = findgen(10)*(MAX(logv31) - MIN(logv31))/10.0 + MIN(logv31)
bestfitx32 = findgen(10)*(MAX(logv32) - MIN(logv32))/10.0 + MIN(logv32)
bestfitx33 = findgen(10)*(MAX(logv33) - MIN(logv33))/10.0 + MIN(logv33)


fit11 = poly_fit(logv11,logbm11,1)
bestfity11 = bestfitx11*fit11[1] + fit11[0]
fit12 = poly_fit(logv12,logbm12,1)
bestfity12 = bestfitx12*fit12[1] + fit12[0]
fit13 = poly_fit(logv13,logbm13,1)
bestfity13 = bestfitx13*fit13[1] + fit13[0]

fit21 = poly_fit(logv21,logbm21,1)
bestfity21 = bestfitx21*fit21[1] + fit21[0]
fit22 = poly_fit(logv22,logbm22,1)
bestfity22 = bestfitx22*fit22[1] + fit22[0]
fit23 = poly_fit(logv23,logbm23,1)
bestfity23 = bestfitx23*fit23[1] + fit23[0]

fit31 = poly_fit(logv31,logbm31,1)
bestfity31 = bestfitx31*fit31[1] + fit31[0]
fit32 = poly_fit(logv32,logbm32,1)
bestfity32 = bestfitx32*fit32[1] + fit32[0]
fit33 = poly_fit(logv33,logbm33,1)
bestfity33 = bestfitx33*fit33[1] + fit33[0]


plot,[10,400],[7e5,1e12],title = 'Baryonic TF',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[10,400],yrange = [1e5,1e12]
oplot,v11,bm11,psym = 4,color = 50
oplot,10.^bestfitx11,10.^bestfity11,color = 50
oplot,v12,bm12,psym = 4, color = 100
oplot,10.^bestfitx12,10.^bestfity12,color = 100
oplot,v13,bm13,psym = 4, color = 240
oplot,10.^bestfitx13,10.^bestfity13,color = 240

oplot,10.^bestfitx21,10.^bestfity21,color = 50,linestyle = 1
oplot,10.^bestfitx22,10.^bestfity22,color = 100,linestyle = 1
oplot,10.^bestfitx23,10.^bestfity23,color = 240,linestyle = 1

oplot,10.^bestfitx31,10.^bestfity31,color = 50,linestyle = 2
oplot,10.^bestfitx32,10.^bestfity32,color = 100,linestyle = 2
oplot,10.^bestfitx33,10.^bestfity33,color = 240,linestyle = 2

legend,['X2X2 - R200','V200*20/220','V200*10/220','X3X3 - R200','V200*20/220','V200*10/220','X5X5 - R200','V200*20/220','V200*10/220','McGaugh 2005'],linestyle=[0,1,2,0,1,2,0,1,2,0],color = [50,50,50,100,100,100,100,240,240,240,0]
stop

end
