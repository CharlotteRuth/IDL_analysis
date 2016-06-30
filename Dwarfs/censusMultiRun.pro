pro censusMultiRun,mstr,mass
dir='fb'
;dir='nofb'
;dir='dESN.0.4'
length=5
;mass=13 ;ie 1e13
directories=['5E1R','1E2R', '1E3R','1E4R','1E5R','1E6R_Rok']
directories=['5E1R','1E2R', '1E3R','1E4R','1E5R']
resolutions=[50,100,1000,10000,100000,1000000]
resolutions=[50,100,1000,10000,100000]
;trials = ['1','2','3','4','5','6','7','8','9']
col=[64,0,254]
m=2.325e5
t=1e9
bin=1e8
outmassfrac=fltarr(length)
bofrac=fltarr(length)
starfrac=fltarr(length)
accfrac=fltarr(length)
notfrac=fltarr(length)
mtot=fltarr(length)
mtol=fltarr(length)
ms=fltarr(length)

;set_plot,'ps', /copy, /interpolate
device,filename=mstr+'censusMulti.eps', /color, bits_per_pixel=8
;DEVICE, DECOMPOSED = 0
loadct,39
!p.multi=[0,3,2]
for i=0,length-1 do begin
;  mstr=strtrim(mass)+'M'
;  mstr='11M'
  mvir=(10.^mass)/m
  rvir=(3.*mvir*m/(4.*!PI*200.*278*0.7^2))^(1./3.)
;  if (i lt 3) then begin trial = 9. 
;  endif else if (i lt 5) then begin  trial = 3.
;  endif else if (i lt 6) then begin trial = 1.
;  endif else begin trial = 0
;  endelse
  trial = 1.0
  s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
  g={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,DENS:0,TEMPG:0,H:0,ZMETAL:0,PHI:0}
  d={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,EPS:0,PHI:0}
  for t = 0, trial-1 do begin
            ssub={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
            gsub={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,DENS:0,TEMPG:0,H:0,ZMETAL:0,PHI:0}
            dsub={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,EPS:0,PHI:0}
            file=directories[i]+'/'+mstr+'/o'+mstr+'__'+trials[t]+'.00300'
       ;     print,file
            rtipsy,file,h,gsub,dsub,ssub
            ;Makes an array of all stars formed in all trials
            if (t ne 0) then concat_structs,s,ssub,ssum else ssum = ssub
            s = ssum
            if (t ne 0) then concat_structs,g,gsub,gsum else gsum = gsub
            g = gsum
            if (t ne 0) then concat_structs,d,dsub,dsum else dsum = dsub
            d = dsum
        ENDFOR

  gr = (g.x^2. + g.y^2. + g.z^2.)^(0.5)
  dr = (d.x^2. + d.y^2. + d.z^2.)^(0.5)
  sr = (s.x^2. + s.y^2. + s.z^2.)^(0.5)
  gv = (g.vx^2. + g.vy^2. + g.vz^2.)^(0.5)
  gv2 = (g.vx^2. + g.vy^2. + g.vz^2.)
; Just use z component since that's the way the wind is blowing
  gdot = g.vz*g.z               ;(g.vx*g.x + g.vy*g.y + g.vz*g.z)
  gdens = g.dens*9.4952922e-3
  gmr=gr
  be = 0.5*gv2 + g.phi          ;mvir/rvir
  vesc = (2*mvir/rvir)^(0.5)
  beind = where(be GT 0, nbeind, comp=boundind)
  accind = where (gdens GT 1e-2 AND g.tempg LT 3e4 AND gdot LE 0)
  outind = where(gdot GT 0 AND g.zmetal GT 0 AND be LE 0 AND g.vz GT 0.1*vesc, noutind)
  mstar=total(s.mass)
  mgas=total(g.mass)
  mdark=total(d.mass)
  ms[i]=mstar*m
  if (nbeind GT 0 ) then outmassfrac[i]=total(g[beind].mass)/(mstar+mgas)$
  else outmassfrac[i]=0 ;nbeind = Fraction of gas beyond mvir/rvir (blow away)
  if (noutind GT 0 ) then bofrac[i]=total(g[outind].mass)/(mstar+mgas)$
  else bofrac[i]=0 ;bofrac = Fraction of gas that will escape perminantly (blow out)
  starfrac[i]=mstar/(mstar+mgas) ;starfrac = fraction of gas as stars
  dummy =  total(accind)
  if dummy eq -1 then accfrac[i]=0 $
  else accfrac[i]=total(g[accind].mass)/(mstar+mgas); accfrac = fraction of gas that is falling toward the nucleus
  print,"accfrac: ",accfrac[i]
  print,"starfrac: ",starfrac[i]
  print,"outmassfrac: ",outmassfrac[i]
  notfrac[i]=1-(starfrac[i]+accfrac[i]+outmassfrac[i])
  mtot[i]=(mstar+mgas+mdark)*m
  mtol[i]=(mstar+mgas+mdark)/mstar
  ;endfor
endfor
 ;legend,['No FB','E!lSN!n=0.2','E!lSN!n=0.4'],line=[0,0,0],color=col,pos=[3.5,80],charsize=1.5
;ppsclose

fbind=[0,1,2,4]
!p.multi=[0,1,1]
zeroIND=where(outmassfrac GT 0)
print,zeroIND,total(outmassfrac[zeroIND])
print,outmassfrac[where(outmassfrac GT 0)]
print,min(outmassfrac[where(outmassfrac GT 0)])
if (total(outmassfrac[zeroIND]) GT -1) then zero=min(outmassfrac[where(outmassfrac GT 0)]) else zero=min(starfrac[where(starfrac GT 0)])
if (min(outmassfrac) EQ 0) then outmassfrac[where(outmassfrac EQ 0)] = zero

;paperps,file='ESN2census.eps',/encapsulated,/color,charsize=2.3
;paperps,file='nofbcensus.eps',/encapsulated,/color,charsize=2.3
;paperps,file='census.eps',/encapsulated,/color,charsize=2.3
plot,resolutions,outmassfrac,xrange=[5e1,1e6],/xlog,xstyle=1, yrange=[zero,1],/ylog,ystyle=1, $
/nodata,xtit='Number of Baryon Particles',ytit='Baryonic Mass Fraction', $
xmargin=[6.5,1.4],ymargin=[3.5,3.5],title='Mass Distribution at Different Resolutions for M='+mstr+sunsymbol()+' Galaxies'
resolutions=[50,resolutions,1000000]
;oplot,mtot,fltarr(n_elements(mtot))+1.,psym=-1
polyfill,resolutions,[zero,fltarr(n_elements(mtot))+1.,zero],color=120
print,[0.0,fltarr(n_elements(mtot))+1.,0.0]
;pstop
;oplot,mtot,starfrac+accfrac+outmassfrac,color=160,psym=-1
polyfill,resolutions,[zero,starfrac+accfrac+outmassfrac+bofrac,zero],color=20
print,[0,starfrac+accfrac+outmassfrac+bofrac,0]
;stop
polyfill,resolutions,[zero,starfrac+outmassfrac+bofrac,zero],color=64
print,[0,starfrac+outmassfrac+bofrac,0]
;stop

;oplot,mtot,starfrac+outmassfrac,color=64,psym=-1
polyfill,resolutions,[zero,outmassfrac+bofrac,zero],color=190
print,[0,outmassfrac+bofrac,0]
;stop
;oplot,mtot,outmassfrac,color=254,psym=-1

polyfill,resolutions,[zero,outmassfrac,zero],color=254
print,[0,outmassfrac,0]
print,resolutions
;stop

;plot,resolutions,outmassfrac,xrange=[0,1e5],/xlog,xstyle=1, $
;;/ylog,yrange=[zero,1],ystyle=1, $
;/nodata,/noerase,xtit='Number of Baryon Particles', $
;ytit='Baryonic Mass Fraction', xmargin=[6.5,1.4],ymargin=[3.5,0.7]
;xyouts,1e2,.01,'Blow Away',charsize=2
;xyouts,1e3,0.1,'Blow Out',charsize=2
;xyouts,1e4,0.5,'Stars',charsize=2
;xyouts,1e4,0.8,'Unaccreted',charsize=2

 ;legend,['No FB','E!lSN!n=0.2','E!lSN!n=0.4'],line=[0,0,0],color=col,/right,charsize=1.7
ppsclose

END
