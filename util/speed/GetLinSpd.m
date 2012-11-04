function spd = GetLinSpd(x,y)

vx = dxdt(x);
vy = dxdt(y);

if length(Range(vx)) ~= length(Range(vy))
   if length(Range(vx)) < length(Range(vy))
      xr = Range(vx); yd = Data(vy);
      vy = tsd(xr,yd(1:length(xr)));
   else
       yr = Range(vy); xd = Data(vx);
       vx = tsd(yr,xd(1:length(yr)));
   end
end

spd = tsd(Range(vx),(sqrt(Data(vx).^2+Data(vy).^2))/2.9);