function nhat = calcnhat(w,h,hmin,A,wa,ha,c1,x1,na,BjerklienOpt)

switch BjerklienOpt
    case 1
        nhat=c1.*( w.*h./wa./ha ).^x1 .* na; %updated if either nau or x1u change
    case 2
        nhat=c1.*( h./ha ).^x1 .* na; %updated if either nau or x1u change
    case 3
        nhat=na .* (h-hmin+.1).^x1; %relative stage
    case 4
       nhat=na .* (A./w).^x1; %hydraulic depth power law
    case 5        
        %note: using x1=sigmad for now. see notes from july 17, little red notebook
        %would be better to have prior on average depth & average depth
        %variation from M&T
        as=1/2;
        ad=5/3;        
        Cd=x1./(A./w);
        Cs2=(Cd.^2+1) .^( (-ad/as).^2 )-1;  %this can get unrealistic
        nhat = na .* (1 + 1/2.* (as./(ad+as).*Cd.^2 + ad./(ad+as).*Cs2 ) );
end

return