function Ham = energies2()
points = 20;
syms x y real;

assign = {1, 0.4, 0.06, 1, 0.05};
[a, pot, soc, tb, r] = assign{:};

Ham = [Hsoc(soc,x,y,a)+Hpot(pot), 0, Htb1(tb,x,y,a)+1i*Htb2(tb,x,y,a), -1i*Hr1(r,x,y,a)-Hr2(r,x,y,a)-Hr3(r,x,y,a)+1i*Hr4(r,x,y,a);... 
       0, -Hsoc(soc,x,y,a)+Hpot(pot), -1i*Hr1(r,x,y,a)+Hr2(r,x,y,a)-Hr3(r,x,y,a)-1i*Hr4(r,x,y,a), Htb1(tb,x,y,a)+1i*Htb2(tb,x,y,a);...
       Htb1(tb,x,y,a)-1i*Htb2(tb,x,y,a), 1i*Hr1(r,x,y,a)+Hr2(r,x,y,a)-Hr3(r,x,y,a)+1i*Hr4(r,x,y,a), -Hsoc(soc,x,y,a)-Hpot(pot), 0;... 
       1i*Hr1(r,x,y,a)-Hr2(r,x,y,a)-Hr3(r,x,y,a)-1i*Hr4(r,x,y,a), Htb1(tb,x,y,a)-1i*Htb2(tb,x,y,a), 0, Hsoc(soc,x,y,a)-Hpot(pot)];
Ham = simplify(Ham);

[X,Y] = meshgrid(0:2.3*pi/(points-1):2.3*pi);

band1 = zeros(points); band2 = zeros(points); band3 = zeros(points); band4 = zeros(points);

for i = 1:1:length(X)^2
    vals = subs(Ham,{x,y},{X(i),Y(i)});
    vals = eig(vals);
    vals = double(vals);
    vals = real(vals);
    vals = sort(vals);
    band1(i) = vals(1); band2(i) = vals(2); band3(i) = vals(3); band4(i) = vals(4);
end
band1 = real(band1); band2 = real(band2); band3 = real(band3); band4 = real(band4);

hold on

surf(X,Y,band1,blue(points));
surf(X,Y,band2,red(points));
surf(X,Y,band3,green(points));
surf(X,Y,band4,yellow(points));

legend('band1','band2','band3','band4')

end

function blue = blue(points)
blue = zeros(points);
blue(:,:,2) = zeros(points);
blue(:,:,3) = ones(points);
end
function red = red(points)
red = ones(points);
red(:,:,2) = zeros(points);
red(:,:,3) = zeros(points);
end
function green = green(points)
green = zeros(points);
green(:,:,2) = ones(points);
green(:,:,3) = zeros(points);
end
function yellow = yellow(points)
yellow = ones(points);
yellow(:,:,2) = ones(points);
yellow(:,:,3) = zeros(points);
end

%spin orbit coupling
function soc = Hsoc(soc,x,y,a)   
x = x*a/2;
y = sqrt(3)*a*y/2;
soc = soc*(2*sin(2*x)-4*sin(x)*cos(y));
end

%chemical potential 
function pot =  Hpot(lamda)
     pot = lamda;
end

%tight binding
function tb1 = Htb1(t,x,y,a)
    x = x*a/2;
    y = sqrt(3)*a*y/2;
    tb1 = t*(1+2*cos(x)*cos(y));
end
function tb2 = Htb2(t,x,y,a) 
    x = x*a/2;
    y = sqrt(3)*a*y/2;
    tb2 = -2*t*cos(x)*sin(y);
end

%Rashba term (tight binding and spin)
function r1 = Hr1(lamda,x,y,a)
    x = x*a/2;
    y = sqrt(3)*a*y/2;
    r1 =  lamda*(1-cos(x)*cos(y));
end

function r2 = Hr2(lamda,x,y,a)
    x = x*a/2;
    y = sqrt(3)*a*y/2;
    r2 = -sqrt(3)*lamda*sin(x)*sin(y);
end

function r3 = Hr3(lamda,x,y,a)
    x = x*a/2;
    y = sqrt(3)*a*y/2;
    r3 = -lamda*cos(x)*sin(y);
end

function r4 = Hr4 (lamda,x,y,a)
    x = x*a/2;
    y = sqrt(3)*a*y/2;
    r4 = sqrt(3)*lamda*sin(x)*cos(y);
end