function [SH SV P S]=rad_pat(strike,dip,rake,f1,f2,fin)

%f1=1; % Hz
%f2=3; % Hz
M=50000;
takeoff=120; % takeoff angles of receiver
% azr=(1:360);
azr=90-(1:360);
daz=azr(2)-azr(1);
fn=0;
for f=fin
    fn=fn+1;
    xn=0;
    for az=azr% station azimuth
        xn=xn+1;
        
        % Average radiation pattern
        t1=takeoff-(180/pi)*(pi/6)*(f-f1)/(f2-f1);
        t2=takeoff+(180/pi)*(pi/6)*(f-f1)/(f2-f1);
        az1=(180/pi)*(1.5*pi/3)*(f-f1)/(f2-f1);
        
        z=rand(M,1);
        n=rand(M,1);
        
        t=acos((1-z)*cosd(t1)+z*cosd(t2))*180/pi;
        p=az+az1*(0.5-n);
        
        [FP FSV FSH FS]=radpattern(strike,dip,rake,t,p);
        SH(xn)=mean(abs(FSH));
        SV(xn)=mean(abs(FSV));
        P(xn)=mean(abs(FP));
        S(xn)=mean(abs(FS));
    end
end
end

function [FP FSV FSH FS]=radpattern(strike,dip,rake,tko,az)
FP=cosd(rake).*sind(dip).*((sind(tko)).^2).*sind(2*(az-strike)) ...
    -cosd(rake).*cosd(dip).*sind(2*tko).*cosd(az-strike) ...
    +sind(rake).*sind(2*dip).*(cosd(tko).^2-((sind(tko).^2).*(sind(az-strike).^2))) ...
    +sind(rake).*cosd(2*dip).*sind(2*tko).*sind(az-strike);

FSV=sind(rake).*cosd(2*dip).*cosd(2*tko).*sind(az-strike) ...
    -cosd(rake).*cosd(dip).*cosd(2*tko).*cosd(az-strike) ...
    +0.5*cosd(rake).*sind(dip).*sind(2*tko).*sind(2*(az-strike)) ...
    -0.5*sind(rake).*sind(2*dip).*sind(2*tko).*(1+sind(az-strike).^2);

FSH=cosd(rake).*cosd(dip).*cosd(tko).*sind(az-strike) ...
    +cosd(rake).*sind(dip).*sind(tko).*cosd(2*(az-strike)) ...
    +sind(rake).*cosd(2*dip).*cosd(tko).*cosd(az-strike) ...
    -0.5*sind(rake).*sind(2*dip).*sind(tko).*sind(2*(az-strike));

FS=sqrt(FSH.^2+FSV.^2);
end