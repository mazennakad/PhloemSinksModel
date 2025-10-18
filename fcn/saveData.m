function saveData(c,u,v,p,nu,dp,psip,psix,pt,maxU,cc,cm,psi_diff,cfl,filename)
% c sugar mass Kg/m2
% u axial velocity m/s
% v radial velocity m/s
% p dynamic pressure MPa
% nu dynamic viscosity Pa/(m s)
% dp dyanmic pressure drop MPa
% psip and psix average phloem and xylem water potential MPa
% pt total pressure MPa
% maxU is the average velocity at each hour m/s
% cc sugar concentration mol/m3
% cm is the average sugar concentration in mol/m3
% psi_diff is the difference in water potential MPa
% cfl is the courant number

writematrix(c,filename,'Sheet','Mass')
writematrix(u,filename,'Sheet','AxialV')
writematrix(v,filename,'Sheet','RadialV')
writematrix(p,filename,'Sheet','DynamicP') % Note this is the total pressure (dynamic + hydrostatic)
writematrix(nu,filename,'Sheet','DynamicV')
writematrix(dp,filename,'Sheet','PressureG')
writematrix(psip,filename,'Sheet','PWP')
writematrix(psix,filename,'Sheet','XWP')
writematrix(pt,filename,'Sheet','TotalP')
writematrix(maxU,filename,'Sheet','MaxU')
writematrix(cc,filename,'Sheet','Concentration')
writematrix(cm,filename,'Sheet','TotalMass')
writematrix(psi_diff,filename,'Sheet','DWP')
writematrix(cfl,filename,'Sheet','cfl')

end