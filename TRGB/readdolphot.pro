function readdolphot,file

return,mrdfits(file,1,columns=['x','y','mag1_acs','mag2_acs'])
end
