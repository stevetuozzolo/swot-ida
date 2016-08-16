function nhat = calcnhat(w,h,hmin,wa,ha,c1,x1,na,BjerklienOpt)

switch BjerklienOpt
    case 1
        nhat=c1.*( w.*h./wa./ha ).^x1 .* na; %updated if either nau or x1u change
    case 2
        nhat=c1.*( h./ha ).^x1 .* na; %updated if either nau or x1u change
    case 3
        nhat=na .* (h-hmin+.1).^x1;
end

return