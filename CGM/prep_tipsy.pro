;filebase = 'h603.cosmo50cmb.3072g14HBWK'
;haloid = '1'
;step = '00512'
;dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir'

;For standard cosmology, redshift to spatial converstion = 1/0.011674743077949343
PRO prep_tipsy,filebase,haloid,step,dir,angle = angle,smoothfile = smoothfile,maxr = maxr,remake = remake

IF NOT keyword_set(angle) THEN angle = [0,45.0,90.0]
IF NOT keyword_set(maxr) THEN maxr = 300 ;kpc
maxv = 300.
filename = filebase + '.' + step
pfile = '../../' + filebase + '.param'
units = tipsyunits(pfile)
outarray =['smoothlength','FeMassFrac','OxMassFrac','HI','H2','iord','HeI','HeII','hsm']
outarray =['FeMassFrac','OxMassFrac','HI','H2']

stat = read_stat_struc_amiga(filename + '.amiga.stat')
ind = where(stat.group EQ haloid)
mdark = string(strtrim((round(10*alog10(stat[ind].m_dark)))/10.0,2),format='(A4)')
mstar = string(strtrim((round(10*alog10(stat[ind].m_star)))/10.0,2),format='(A4)')
outfile = strtrim(haloid,2) + '_' + mdark + '_' + mstar + '_'
shortbase = (strsplit(filebase,'.',/extract))[0] + '.g' + (strsplit((strsplit(filebase,'.',/extract))[2],'g',/extract))[1]
LOSfile = 'LOS.' + shortbase + '.' + step + '.' + strtrim(haloid,2)
first = strpos(filebase,'cosmo') + 5
last = strpos(filebase,'cmb')
box = strmid(filebase,first,last - first)
h = 0.7
box = strtrim(box*h,2)
halofilename = filename + '.halo.' + strtrim(string(haloid),2) + '.' + strtrim(maxr,2)
IF (file_exist(halofilename + '.std') EQ 0) OR keyword_set(remake) THEN BEGIN
;Select the halo
    tipsysatshi,filename,haloid,units.lengthunit,units.massunit,cutout_rad = maxr,outarray = outarray
ENDIF
rtipsy,halofilename + '.std',h,g,d,s,/justhead
redshift = string(strtrim(round(100.0*(1/h.time - 1)/100.0),2),format='(A5)')
print,'Run ~/Code/smooth/smooth -s 32g -o ' + halofilename + ' hsmooth < ' + dir + filebase + '.' + step + '/'+ halofilename + '.std'
;print,'Run ~/tipsy_tools/smooth -s 32g -o ' + halofilename + ' hsmooth < ' + dir + filebase + '.' + step + '/'+ halofilename + '.std'
;stop ;Take out stop after having created the .hsm file
;Add the correct smoothing lengths to the file
fix_smooth,halofilename         ;Missing so temporarily commented out

;Create the auxilary files
IF (file_exist(halofilename + '.aux') EQ 0) OR keyword_set(remake) THEN $
  prep_tipsyaux,halofilename,pfile

;Rotate the tipsy files
FOR i = 0, n_elements(angle) - 1 DO $
  IF angle[i] NE 0 THEN $
;    IF (file_exist(halofilename + '.' + strtrim(fix(angle[i]),2)) EQ 0) OR keyword_set(remake) THEN $
      rotatetipsy,halofilename + '.std',angle[i]

;Create the Line-of-sight files
impact_param = [2.0,4.0,8.0,14.0,20.0,50.0,100.0,150.0,200.0,250.0,300.0]
nangle = 32
FOR i = 0, n_elements(angle) - 1 DO BEGIN
  makeLOS,LOSfile + '.' + strtrim(fix(angle[i]),2),outfile + strtrim(fix(angle[i]),2),impact_param,units.lengthunit,nangle
ENDFOR

;Create shell script
openw,1,halofilename + '.sh',1
FOR i = 0, n_elements(angle) - 1 DO BEGIN
    printf,1,"/home/christensen/Code/specexbin_merge/specexbin " $
           + dir + '/CGM/' + halofilename + '.' + strtrim(fix(angle[i]),2) + ' ' $
           + dir + '/CGM/' + LOSfile + '.' + strtrim(fix(angle[i]),2) + ' ' $
           + strtrim(redshift,2) + ' ' $
           + box $
           + ' 0.24 0.70 0.26 1 1 1' ; > ' $
;           + halofilename + '.' + strtrim(fix(angle[i]),2) + '.out'
ENDFOR
close,1

;Create the .tab file
c = 3e5 ;km/s
dz = maxv*(1 + redshift)/c
FOR i = 0, n_elements(angle) - 1 DO BEGIN
    openw,1,halofilename + '.' + strtrim(fix(angle[i]),2) + '.tab',1
    printf,1,dir + '/CGM/' + halofilename + '.' + strtrim(fix(angle[i]),2) + ' ' $
           + strtrim(redshift,2) + ' ' $
           + strtrim(-1*dz,2)  + ' ' $
           + strtrim(dz,2) + ' ' $
           + shortbase   
    close,1
ENDFOR

END

PRO prep_tipsy_master
  dir = '/nobackupp8/crchrist/MolecH/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
  filebase = 'h239.cosmo50cmb.3072g14HMbwK'
  step_arr = ['00048','00060','00084','00128','00228','00328','00424']
  haloid_arr = ['6','5','2','1','1','1','1']

  FOR i = 0, n_elements(step_arr) - 1 DO BEGIN
     cd,dir + filebase + '.' + step_arr[i]
     pwd
     prep_tipsy,filebase,haloid_arr[i],step_arr[i],dir,angle = ['0.0']
  ENDFOR
END
