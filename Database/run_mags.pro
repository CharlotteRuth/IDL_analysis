pro run_mags, run=run
cd,'/unsafe-scratch/volatile/sloebman/DARREN/'

if (run eq 256) then begin
    ;256 FILE
    tipsy    = 'cosmo50cmb.256g1MBWK.00048'
    i_mass   = 7.483359e-10     ;256
    ;snap     = 0.0341187        ;256
    print, tipsy
endif else begin
    if (run eq 768) then begin
        ;768 FILE
        tipsy    = 'cosmo50cmb.768pg1MBWK.00048'
        i_mass   = 2.771614e-11 ;768
        ;snap     = 0.0334894    ;768 max(s.tform)
        print, tipsy
    endif
endelse

;AMIGA FILES
statfile = tipsy + '.amiga.stat'
grpfile  = tipsy + '.amiga.grp'
gtpfile  = tipsy + '.amiga.gtp' 

;SAME FOR BOTH RUNS
m_unit   = 18.48e15
;l_unit   = 50000. ;KPC
l_unit   = 50. ;MPC

;GET HALO MAGS PARAMETERS
iset     = 0 ;KROUPA
hfinder  = 'amiga'
mag      = 'ab'

;print, tipsy
;print, mag
;print, ';-------------------'

get_halo_mags, m_unit, l_unit, $
               init_mass=i_mass, $
               iso_set=iset, $ ;snaptime=snap, $
               mag_sys=mag, /multiple, $
               halo_finder=hfinder
stop

end
