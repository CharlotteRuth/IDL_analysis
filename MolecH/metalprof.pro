;
;
;
;
;
;
;
;
; create metalicity profiles
;
;
;
;



pro metalprof, file, list, gaslist, binsize, h=h, g=g, s=s, normalize=normalize

if not keyword_set(h) then begin
   rtipsy, file, h, g, d, s
   cofm, h, g, d, s, /pot
endif

set_plot, 'ps'
device, filename='metalprofiles.ps', /color, ysize = 4, xsize = 7, /inches

!p.multi = [0,2,1]

XSOLfe=0.117E-2
XSOLO=0.96E-2
XSOLH=0.706

; **********************************************************************
; get the metalicites of particles - note, you need the auxiliary files
; **********************************************************************

if keyword_set(normalize) then $
  normalize_metals, list, file, h, g, s, fehs, fehg $
else $
  metals, file, fehs, fehg, ofes, ofeg, h=h, g=g, s=s

;stop

fehs = fehs*(10^0.2)
fehg = fehg*(10^0.2)

; ************************
; make the radial profile 
; ************************

age = h.time - s.tform

age1 = where(age lt 0.01)
age2 = where(age gt 0.5 and age lt 3.0)
age3 = where(age gt 3.0)

rs = sqrt(s.x^2.+s.y^2.)
rg = sqrt(g.x^2.+g.y^2.)

ind = where(rg lt 15. and abs(g.z) lt 1.0) 
g = g[ind]
rg = rg[ind]

;
; gas metallicity profile
;

masshistg = hist1d(rg, g.mass, binsize = binsize)
meanfeg   = hist1d(rg, fehg*g.mass, binsize = binsize, $
                   obin = xg, stddev = stddevg)
meanfeg   = meanfeg/masshistg
meanfehg  = alog10(meanfeg)-alog10(XSOLFe/XSOLH)   

;
; all stars profile
;

masshists = hist1d(rs, s.mass, binsize = binsize)
meanfes   = hist1d(rs, fehs*s.mass, binsize = binsize, $
                   obin = xs, stddev = stddev_all)
meanfes   = meanfes/masshists
meanfehs  = alog10(meanfes)-alog10(XSOLFe/XSOLH)   

;
; young stars profile
;

masshists1 = hist1d(rs[age1], s[age1].mass, binsize = binsize)
meanfes1   = hist1d(rs[age1], fehs[age1]*s[age1].mass, $
                    stddev = stddev1, binsize = binsize, obin = xs1)
meanfes1   = meanfes1/masshists1
meanfehs1  = alog10(meanfes1)-alog10(XSOLFe/XSOLH)   

;
; intermediate age stars profile
;

masshists2 = hist1d(rs[age2], s[age2].mass, binsize = binsize)
meanfes2   = hist1d(rs[age2], fehs[age2]*s[age2].mass, $
                    stddev = stddev2, binsize = binsize, obin = xs2)
meanfes2   = meanfes2/masshists2
meanfehs2  = alog10(meanfes2)-alog10(XSOLFe/XSOLH)   

;
; old stars profile
;

masshists3 = hist1d(rs[age3], s[age3].mass, binsize = binsize)
meanfes3   = hist1d(rs[age3], fehs[age3]*s[age3].mass, $
                    stddev = stddev3, binsize = binsize, obin = xs3)
meanfes3   = meanfes3/masshists3
meanfehs3  = alog10(meanfes3)-alog10(XSOLFe/XSOLH)   

multiplot

plot, xg, meanfehg, title = 'Stars', $
  xtitle = 'R [kpc]', ytitle = '[Fe/H]', yrange = [-0.4,0.5], $
psym = 10, /nodata, xrange = [0,15]
oplot, xs, meanfehs, line = 5
oplot, xs1, meanfehs1, color = 50, psym = 10
oplot, xs2, meanfehs2, color = 150, psym = 10
oplot, xs3, meanfehs3, color = 250, psym = 10

legend, /right, ['all stars', 'age < 0.1 Gyr', '1 < age < 3 Gyr', 'age > 5 Gyr'], $
  linestyle = [5,0,0,0], color = [0,50,150,250], charsize = 1


;
; get the gas metallicity profiles at different times
;

readcol, gaslist, format='A60', files

legend_titles = strarr(n_elements(files))
legend_lines = fltarr(n_elements(files))
legend_color = fltarr(n_elements(files))

color = [250,150,50]

for i=0, n_elements(files)-1 do begin
    
    rtipsy, files[i], h, g, d, s
    cofm, h, g, d, s, /pot
    d=0

    metals, files[i], fehs, fehg, ofes, ofeg, gas, disk_only = 15., h=h,g=g,s=s
    fehg = fehg*(10^0.2)

    rg = sqrt(gas.x^2+gas.y^2)

    binsize = 0.5

    masshistg = hist1d(rg, gas.mass, binsize = binsize)
    meanfeg   = hist1d(rg, fehg*gas.mass, binsize = binsize, $
                       obin = xg, stddev = stddevg)
    meanfeg   = meanfeg/masshistg
    meanfehg  = alog10(meanfeg)-alog10(XSOLFe/XSOLH) 

    if(i eq 0) then begin
        multiplot
        
        plot, xg, meanfehg, xtitle = 'R [kpc]', title = 'Gas', $
          yrange=[-0.4,0.5], xrange = [0,15], /nodata
        oplot, xg, meanfehg, color = color[i], psym = 10
    endif else begin
        oplot, xg, meanfehg, color = color[i], psym = 10
    endelse

    legend_titles[i] = strtrim(string(h.time, format='(F5.1)'),2)+' Gyr'
    legend_color[i]  = color[i]

endfor

legend, legend_titles, color = legend_color, line = legend_lines, /right, charsize = 1

device, /close
set_plot, 'x'
!p.multi = 0

multiplot, /reset
multiplot, /default

end
