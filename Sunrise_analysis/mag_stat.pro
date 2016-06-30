pro mag_stat, Xname, gal_list, scat

;if scat = 1, then scattering is on
;if scat is not specified, then scattering is off

if not keyword_set(scat) then scat = 0 
if (scat eq 0) then begin
;read in gal_list and sort in ascending order.  index is sorted array
readcol, gal_list, F='I',ind
results = bsort(ind, index)

;file to write to.  need to open it to write to before entering loop
outfile = Xname + '.unscat.mag'
openw, 1, outfile, width=300

;reading in the name of the filters to put into the top of the .mag file
readcol, 'filters', F='A', filters

trimmed = STRMID(filters,0,5)
trim_scat = '          ' + trimmed 
filter_header = [trim_scat]

printf, 1, '      id', filter_header

for i=0,n_elements(index)-1 do begin

direct = Xname + '_' + strn(index[i])
cd, direct

openr, 2, 'broadband_rf.fits', Error = err
close, 2

if (err eq 0) then begin

file = 'broadband_rf.fits'
filt = mrdfits(file,10, Error=err2)   


ab_mag0_nscat = (filt[*]).AB_MAG_NONSCATTER0

printf,1,index[i], ab_mag0_nscat 

cd, '../'
endif else begin
cd, '../'
endelse
 
endfor
close,1
endif else begin

;read in gal_list and sort in ascending order.  index is sorted array
readcol, gal_list, F='I',ind
results = bsort(ind, index)

;file to write to.  need to open it to write to before entering loop
outfile = Xname + '.scat.mag'
openw, 1, outfile, width=300

;reading in the name of the filters to put into the top of the .mag file
readcol, 'filters', F='A', filters

trimmed = STRMID(filters,0,5)
trim_scat = '          ' + trimmed 
filter_header = [trim_scat]

printf, 1, '      id', filter_header

for i=0,n_elements(index)-1 do begin

direct = Xname + '_' + strn(index[i])
cd, direct

openr, 2, 'broadband_rf.fits', Error = err
close, 2

if (err eq 0) then begin

file = 'broadband_rf.fits'
filt = mrdfits(file,10, Error=err2)   


ab_mag0_scat = (filt[*]).ab_mag0
printf,1,index[i],ab_mag0_scat 

cd, '../'
endif else begin

cd, '../'
endelse
 
endfor
close,1
endelse

end

