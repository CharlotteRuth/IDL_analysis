FUNCTION convertheader,header,max,min
s = strarr(28)
s[0] = "CTYPE1  = 'KPC-LIN'            /"
s[1] = "CUNIT1  =                  KPC /"
s[2] = "CRVAL1  =            "+strtrim(header.crval1,2)+" /"
s[3] = "CDELT1  =            "+strtrim(header.cdelt1,2)+" /" 
s[4] = "CRPIX1  =            "+strtrim(header.crpix1,2)+" /"
s[5] = "CTYPE2  = 'KPC-LIN'            /"
s[6] = "CUNIT2  =                  KPC /"
s[7] = "CRVAL2  =             "+strtrim(header.crval2,2)+" /"
s[8] = "CDELT2  =             "+strtrim(header.cdelt2,2)+" /"
s[9] = "CRPIX2  =             "+strtrim(header.crpix2,2)+" /"
s[10] = "CTYPE3  = 'VELO-W2V'           /"
s[11] = "CUNIT3  = 'KM/S'               /"
s[12] = "CRVAL3  =          "+strtrim(header.crval3,2)+" /"
s[13] = "CDELT3  =          "+strtrim(header.cdelt3,2)+" /"
s[14] = "CRPIX3  =          "+strtrim(header.crpix3,2)+" /"
;s[15] = "BLANK   =          "+strtrim(header.blank,2)+" /"
;s[16] = "BSCALE  =          "+strtrim(header.bscale,2)+" /"
;s[17] = "BZERO   =          "+strtrim(header.bzero,2)+" /"
s[18] = "BUNIT   = 'solMass'            /"
s[19] = "DATAMAX =           "+strtrim(max,2)+" /"
s[20] = "DATAMIN =           "+strtrim(min,2)+" /"
return,s
END
