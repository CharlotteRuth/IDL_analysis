pro tempgcheck
;This runs through all of the IC galaxy files to see if any have
;tempg==0

r=['5E1R', '1E2R', '1E2R', '1E3R', '1E4R', '1E5R']
m=['7M', '8M', '9M', '10M', '11M', '12M', '13M']

for rct=0,5 do begin
    cd,r[rct]
    for mct=0,6 do begin
        cd,m[mct]
        file=m[mct]+'.std'
        rtipsy,file,h,g,d,s
        if min(g.tempg) eq 0 then print,'tempg=0 at '+r[rct]+'/'+m[mct]
        cd,'..'
    endfor
 cd,'..'
endfor
END
