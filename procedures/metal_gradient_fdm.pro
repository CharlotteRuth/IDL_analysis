pro metal_gradient,nickname

  radial_binsize = 0.25;kpc
  X_sol = 0.0122; (Asplund,M 2009)
  ;fullrunname = fullname(nickname)
  ;establish directories
  ;set_machine,fullrunname,run_dir,data_dir
  ;check_dirs,run_dir,data_dir

  ;establish input files
  ;starlog = run_dir+fullrunname+'.starlog'
  ;satsidfile = data_dir+fullrunname+'.sats.trace_stars_mod.dat'
  ;check_files,starlog
 ; readcol, satsidfile, snapshots, currentid, origid,format='a,l,l,x',/silent

  ;The haloidfile may contain more than one group.  Limit to the one you
  ;want.
  ;grp_id = origid

  ;Units
 ;get_sim_units,snapshots[0],lunit,vunit,hnot,munit

  ;read simulation particle numbers by z=0
run_dir='/home/munshi/dragon-scratch/'
snapshots=strarr(1)
snapshots[0]='h937.cosmo25cmb.4096g1MbwK1C52.004096'
  z0file = run_dir+snapshots[0]
  rtipsy,z0file,h,g,d,s
  delvarx,g,d
grp_id=intarr(1)
grp_id[0]=3
  ;the haloid for each particle at that snapshot
  haloid = read_lon_array(run_dir+snapshots[0]+'.amiga.grp')
  ;the haloid for star particles only
  haloid = haloid[h.ngas+h.ndark:h.n-1]
  ;find the star particles in the halo you chose
  keep = where(haloid eq grp_id[0],/L64)

  ;the halo index for centering the stars
  readcol,run_dir+snapshots[0]+'.amiga.stat',hid,$
    format='l,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x',/silent

  hchose = where(hid eq  grp_id[0])

  ;comoving halo data
  rtipsy,run_dir+snapshots[0]+'.amiga.gtp',h1,g1,d1,sh
  delvarx,g1,d1
get_sim_units,snapshots[0],lunit,vunit,hnot,munit
  ;convert to physical units at z=0
  dsx0 = float(s[keep].x-sh[hchose].x)*lunit*1.00;kpc
  dsy0 = float(s[keep].y-sh[hchose].y)*lunit*1.00;kpc
  dsz0 = float(s[keep].z-sh[hchose].z)*lunit*1.00;kpc

  r0 = sqrt(dsx0^2+dsy0^2+dsz0^2)
  res = alog10(s[keep].metals*2./X_sol)

  make_hist,r0,rh0,rb0,10,0,0.1,ri

  mean_res = fltarr(n_elements(rb0))
  ;mean_time = fltarr(n_elements(rb0))

  ;bin time and metallicity as a function of radius
  for i=0,n_elements(rb0)-1 do begin
    if ri[i] ne ri[i+1] then begin
      res_here = res[ri[ri[i]:ri[i+1]-1]]
      res_good = res_here[where(res_here lt 0 and res_here gt -6,/L64)]
      mean_res[i] = avg(res_good)
      ;mean_time[i]= avg(gtt[ri[ri[i]:ri[i+1]-1]])
    endif
  endfor

  result = linfit(rb0[0:39],mean_res[0:39])
set_plot, 'ps'
device, filename='h937_h3_metgrad.eps', /encapsulated, /color
   ex = create_struct('charsize', 1.5, 'ytitle', '[Fe/H]', 'xtitle', 'Radial Distance [kpc]', 'thick', 4, 'charthick', 4, 'position', [.20,.20,.95,.95],'xthick',2,'ythick',2,'c_colors',[50,80,110,140])

  contour_plus, r0, res, xmin=0, xmax=5, ymin=-4, ymax=1, levels = [10,50,100,500], threshold =10, xbinsize=0.1, ybinsize=0.1, _EXTRA=ex

  oplot,rb0[0:49],result[1]*rb0[0:49]+result[0],thick=6
device, /close
set_plot, 'x'

stop

end
