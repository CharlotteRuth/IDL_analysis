pro test_ct
loadct,27

plot,[0,0],[0,1],xrange=[0,1],yrange=[0,255],ystyle=1
for i = 0, 255 do begin
    oplot,[0,1],[i,i],color = i,thick=2
endfor
oplot,[0.5,0.5][0,255],color=
end
