;Reads in a tipsy simualtion file and starlog file and then outputs an
;array file with the selected output from the starlog file

;Charlotte Christensen
;2/26/2018

;pfile -- name of paramfile
;tfile -- name of the tipsy file
;slfile -- name of the starlog file
;big -- set to 1 in using 64 bit format
;molecularH -- set to 1 if molecular hydrogen is included
;cosmo -- set to 1 if you are using a cosmological simulation

;Relies on idl4tipsy/writetipsyarray.pro, procedures/rtipsy.pro,
;idl4tipsy/read_tipsy_arr.pro,
;/home/christenc/Code/IDL/idlutils/goddard/pro/misc/match.pro, procedures/tipsyunits.pro

;Ex: write_starlog_arr,'h516.cosmo25cmb.2304g1HBWK.000156','h516.cosmo25cmb.2304g1HBWK.starlog','rhoform',/big,/molecularH

PRO write_starlog_arr,pfile,tfile,slfile,arr,big = big,molecularH = molecularH, cosmo = cosmo

  rtipsy,tfile,h,g,d,s
  units = tipsyunits(pfile,silent =1)
  read_tipsy_arr,tfile + '.iord',h,iord,type = 'long'
  sl = rstarlog(slfile,molecularH = molecularH,big = big)
  sl = sl[UNIQ(sl.iorderstar, SORT(sl.iorderstar))]
  match,sl.iorderstar,iord,ind_sl,ind_iord
  IF keyword_set(cosmo) THEN a = h.time ELSE a = 1
  IF arr EQ 'iordergas' THEN BEGIN
     array = intarr(n_elements(iord))
     array[ind_iord] = sl[ind_sl].iordergas
     writetipsyarray,tfile + '.' + arr, array, type = 'long'
  ENDIF ELSE BEGIN
     array = fltarr(n_elements(iord))
     CASE 1 OF
        arr EQ 'timeform': array[ind_iord] = sl[ind_sl].timeform*units.timeunit ;yrs
        arr EQ 'x': array[ind_iord] = sl[ind_sl].x*units.lengthunit*a ;kpc
        arr EQ 'y': array[ind_iord] = sl[ind_sl].y*units.lengthunit*a ;kpc
        arr EQ 'z': array[ind_iord] = sl[ind_sl].z*units.lengthunit*a ;kpc
        arr EQ 'vx': array[ind_iord] = sl[ind_sl].vx*units.vunit*a ;km/s
        arr EQ 'vy': array[ind_iord] = sl[ind_sl].vy*units.vunit*a ;km/s
        arr EQ 'vz': array[ind_iord] = sl[ind_sl].vz*units.vunit*a ;km/s
        arr EQ 'massform': array[ind_iord] = sl[ind_sl].massform*units.massunit ;Msol
        arr EQ 'rhoform': array[ind_iord] = sl[ind_sl].rhoform*units.rhounit/a^3 ;amu/cm^3
        arr EQ 'tempform': array[ind_iord] = sl[ind_sl].tempform ;K
        arr EQ 'H2form': array[ind_iord] = sl[ind_sl].H2form ;fraction
        ELSE: print,arr+' is not an available key in the starlog file'
     ENDCASE
     writetipsyarray,tfile + '.' + arr, array, type = 'float'
  ENDELSE
  
END
