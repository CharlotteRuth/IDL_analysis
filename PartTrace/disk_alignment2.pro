;Charlotte Christensen
;5/25/12
;Adds the virial radius to the halo data found in alignment.txt (from
;disk_alignment.pro) and outputs to a fits file.
;That fits file is used in gas_gal_disk_history to establish when the
;gas particles are within the virial radius.

PRO disk_alignment2,finalid = finalid,debug = debug
IF NOT keyword_set(finalid) THEN finalid = '1'
readcol,'grp' + finalid + '.alignment.txt',files,haloid,time,z,xc,yc,zc,vxc,vyc,vzc,xa,ya,za,format = '(A,I,F,F,F,F,F,F,F,F,F,F,F)'
halodat = {file: ' ', haloid: 0, time: 0D, z: 0D, xc: 0D, yc: 0D, zc: 0D, vxc: 0D, vyc: 0D, vzc: 0D, xa: 0D, ya: 0D, za: 0D, rvir: 0D, mvir: 0D}
;readcol,'grp' + finalid + '.alignment.txt',files,haloid,time,z,xc,yc,zc,xa,ya,za,format = '(A,I,F,F,F,F,F,F,F,F)'
;halodat = {file: ' ', haloid: 0, time: 0D, z: 0D, xc: 0D, yc: 0D, zc: 0D, xa: 0D, ya: 0D, za: 0D, rvir: 0D}

halodat = replicate(halodat,n_elements(files))
;rvir = fltarr(n_elements(files))
;files = files + '/' + files
mvir = fltarr(n_elements(files))

FOR i = 0, n_elements(files) - 1 DO BEGIN
;   amigastat = read_stat_struc_amiga(files[i] +'/' + files[i] + '.amiga.stat')
   amigastat = read_stat_struc_amiga(files[i] + '.amiga.stat')
   ind = where(amigastat.group EQ haloid[i])
   halodat[i].rvir = amigastat[ind].rvir
   print,files[i],haloid[i],halodat[i].rvir,amigastat[ind].m_tot
   mvir[i] = amigastat[ind].m_tot
   halodat[i].mvir = mvir[i]
   halodat[i].file = files[i]
   halodat[i].haloid = haloid[i]
   halodat[i].time = time[i]
   halodat[i].z = z[i]
   halodat[i].xc = xc[i]  
   halodat[i].yc = yc[i]
   halodat[i].zc = zc[i]
   halodat[i].vxc = vxc[i]
   halodat[i].vyc = vyc[i]
   halodat[i].vzc = vzc[i]
   halodat[i].xa = xa[i]  
   halodat[i].ya = ya[i]
   halodat[i].za = za[i] 
ENDFOR
IF keyword_set(debug) THEN BEGIN
   set_plot,'x'
   plot,time,mvir,psym = -4,xtitle = 'Time',ytitle = 'Virial Mass'
   stop
   plot,time,halodat.rvir,psym = -4,xtitle = 'Time',ytitle = 'Virial Radius'
   stop
ENDIF
mwrfits,halodat,'./grp' + finalid + '.alignment.fits',/create
END
