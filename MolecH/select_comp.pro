;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.00512.dir'
;file = 'h603.cosmo50cmb.3072gs1MbwK'
;comp_i = 0
PRO select_comp,file,comp_i
comp = [1,2,3,4,5,0]
comp_st = ['disk','halo','bulge','thdisk','pb','extra']

rtipsy,file,h,g,d,s,/justhead
read_tipsy_arr,file + '.decomp',h,decomp,type = 'int'

FOR ic = 0, N_ELEMENTS(comp_i) - 1 DO BEGIN
    ind = where(decomp eq comp[comp_i[ic]]) + 1 

    openw,1,file + '.decomp.' + comp_st[comp_i[ic]] + '.mrk'
    printf,1,h.n,h.ngas,h.nstar
    FOR i = 0L, N_ELEMENTS(ind) - 1 DO printf,1,ind[i]
    close,1
ENDFOR
END
