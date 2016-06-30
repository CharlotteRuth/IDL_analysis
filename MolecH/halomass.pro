PRO halomass,outfile = outfile
loadct,39
if (KEYWORD_SET(outfile)) then begin
    !Y.STYLE = 3
    !X.STYLE = 3
    !P.THICK = 3.5
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !P.charsize=1.0
    !x.charsize=2.25
    !y.charsize=2.25
    set_plot,'ps' 
endif else set_plot,'x'


prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
param = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/h516.cosmo25cmb.1536g1MBWK.param'
units = tipsyunits(param)

base0 = "h516.cosmo25cmb.1536g14HBWK"
base1 = "h516.cosmo25cmb.1536g3HBWK"
base2 = "h516.cosmo25cmb.1536g1MBWK"
base3 = "h516.cosmo25cmb.2304g14HBWK"
step0 = ['00024','00048','00072','00092','00096','00100','00103','00104','00108','00112','00116','00120','00128','00132','00136','00140','00144','00148','00152','00156','00160','00164','00167','00168','00172','00176','00180','00184','00188','00192','00196','00216','00240','00264','00288','00312','00336','00360','00384','00408','00432','00456','00480','00504','00512']

step1 = ['00030','00045','00060','00075','00090','00105','00120','00135','00150','00165','00180','00195','00210','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00468','00480','00492','00504','00512']

step2 = ['00024','00037','00046','00048','00061','00072','00084','00096','00120','00128','00144','00168','00192','00216','00227','00240','00264','00271','00288','00312','00328','00336','00360','00384','00406','00408','00432','00455','00456','00480','00504','00512']

step3 = ['00024','00036','00037','00046','00048','00060','00061','00072','00084','00096','00108','00120','00132','00144','00148','00160','00172','00184','00196','00212','00224','00232','00240','00244','00256','00268','00272','00276','00288','00300','00308','00312','00324','00336','00348','00360','00372','00384','00396','00406','00416','00428','00440','00448']

keys = ['new H2, lr','old H2','no H2']
colors= [240,100,50]
if (KEYWORD_SET(outfile)) then outfile = '/astro/store/student-scratch1/christensen/MolecH/results/' + 'h516.cosmo25cmb.1536g14HBWK.eps'

keys = ['new H2','no H2']
colors= [240,100,60]
if (KEYWORD_SET(outfile)) then outfile = '/astro/store/student-scratch1/christensen/MolecH/results/' + 'h516.cosmo25cmb.1536g14H_1MBWK.00512_mass.eps'

keys = ['new H2, mr','new H2, lr','old H2','no H2']
colors= [240,100,50,210]
if (KEYWORD_SET(outfile)) then outfile = '/astro/store/student-scratch1/christensen/MolecH/results/' + 'h516.cosmo25cmb.1536g14HBWK.eps'

cd,prefix

data0 = fltarr(N_ELEMENTS(step0),4)
time0 = step0
FOR i = 0, N_ELEMENTS(step0) - 1 DO BEGIN 
    filename0 = prefix + 'h516.cosmo25cmb.1536g/' + base0 + "/steps/"+base0+"."+step0[i]+".dir/" +  base0 + "." + step0[i]
    data = read_stat_struc_amiga(filename0+'.amiga.stat')
    rtipsy,filename0+'.halo.1.std',h,g,d,s
    time0[i] = MAX(s.tform*units.timeunit/1e9)
    data0[i,0] = data[0].m_tot
    data0[i,1] = data[0].m_gas
    data0[i,2] = data[0].m_star
    data0[i,3] = data[0].m_tot - data[0].m_gas - data[0].m_star
ENDFOR

data1 = fltarr(N_ELEMENTS(step1),4)
time1 = step1
FOR i = 0, N_ELEMENTS(step1) - 1 DO BEGIN 
    filename1 = prefix + 'h516.cosmo25cmb.1536g/' + base1 + "/steps/"+base1+"."+step1[i]+".dir/" +  base1 + "." + step1[i]
    data = read_stat_struc_amiga(filename1+'.amiga.stat')
    rtipsy,filename1+'.halo.1.std',h,g,d,s
    time1[i] = MAX(s.tform*units.timeunit/1e9)
    data1[i,0] = data[0].m_tot
    data1[i,1] = data[0].m_gas
    data1[i,2] = data[0].m_star
    data1[i,3] = data[0].m_tot - data[0].m_gas - data[0].m_star
ENDFOR

data2 = fltarr(N_ELEMENTS(step2),4)
time2 = step2
FOR i = 0, N_ELEMENTS(step2) - 1 DO BEGIN 
    filename2 = prefix + 'h516.cosmo25cmb.1536g/' + base2 + "/steps/"+base2+"."+step2[i]+".dir/" +  base2 + "." + step2[i]
    data = read_stat_struc_amiga(filename2+'.amiga.stat')
    rtipsy,filename2+'.halo.1.std',h,g,d,s
    time2[i] = MAX(s.tform*units.timeunit/1e9)
    data2[i,0] = data[0].m_tot
    data2[i,1] = data[0].m_gas
    data2[i,2] = data[0].m_star
    data2[i,3] = data[0].m_tot - data[0].m_gas - data[0].m_star
ENDFOR

data3 = fltarr(N_ELEMENTS(step3),4)
time3 = step3
FOR i = 0, N_ELEMENTS(step3) - 1 DO BEGIN 
    filename3 = prefix + 'h516.cosmo25cmb.2304g/' + base3 + "/steps/"+base3+"."+step3[i]+".dir/" +  base3 + "." + step3[i]
    data = read_stat_struc_amiga(filename3+'.amiga.stat')
    rtipsy,filename3+'.halo.1.std',h,g,d,s
    time3[i] = MAX(s.tform*units.timeunit/1e9)
    data3[i,0] = data[0].m_tot
    data3[i,1] = data[0].m_gas
    data3[i,2] = data[0].m_star
    data3[i,3] = data[0].m_tot - data[0].m_gas - data[0].m_star
ENDFOR

if (KEYWORD_SET(outfile)) then device,filename = outfile,/color,bits_per_pixel= 8,/times,xsize = 32,ysize = 12,xoffset =  2,yoffset =  2  else stop

!p.multi = [0,2,1]
plot,time0, data0[*,3], psym = -2, xstyle = 1, xrange = [0,MAX([time0,time2])], xtitle = 'Time [Gyr]', ytitle = 'Mass [M'+sunsymbol()+']'
oplot,time0,data0[*,3],psym = -2,color = colors[0]
oplot,time0,data0[*,1] + data0[*,2],psym = -1,color = colors[0]

oplot,time1,data1[*,3],psym = -2,color = colors[1]
oplot,time1,data1[*,1] + data1[*,2],psym = -1,color = colors[1]

oplot,time2,data2[*,3],psym = -2,color = colors[2]
oplot,time2,data2[*,1] + data2[*,2],psym = -1,color = colors[2]

oplot,time3,data3[*,3],psym = -2,color = colors[3]
oplot,time3,data3[*,1] + data3[*,2],psym = -1,color = colors[3]
;legend,keys,color = colors,linestyle = fltarr(N_ELEMENTS(colors))
;legend,['Dark Matter','Baryonic Mass'],color = [0,0],psym = [-2,-1],/right

plot,time0,data0[*,1],psym = -4,xstyle = 1, xrange = [0,MAX([time0,time2])], xtitle = 'Time [Gyr]', ytitle = 'Mass [M'+sunsymbol()+']' 
oplot,time0,data0[*,1],psym = -4,color = colors[0]
oplot,time0,data0[*,2],psym = -5,color = colors[0]

oplot,time1,data1[*,1],psym = -4,color = colors[1]
oplot,time1,data1[*,2],psym = -5,color = colors[1]

oplot,time2,data2[*,1],psym = -4,color = colors[2]
oplot,time2,data2[*,2],psym = -5,color = colors[2]

oplot,time3,data3[*,1],psym = -4,color = colors[3]
oplot,time3,data3[*,2],psym = -5,color = colors[3]
legend,['Gas','Stars'],color = [0,0],psym = [-4,-5],/right
!p.multi = 0 
if (KEYWORD_SET(outfile)) then device,/close
END
