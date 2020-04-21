function Ham = energies()
x = sym('x');
y = sym('y');

assign = {1, 0.1, 0.05, 1, 0.05};
[a, pot, soc, tb, r] = assign{:};

Ham = [Hsoc(soc,x,y,a)+Hpot(pot), 0, Htb1(tb,x,y,a)+1i*Htb2(tb,x,y,a), -1i*Hr1(r,x,y,a)-Hr2(r,x,y,a)-Hr3(r,x,y,a)+1i*Hr4(r,x,y,a);... 
       0, -Hsoc(soc,x,y,a)+Hpot(pot), -1i*Hr1(r,x,y,a)+Hr2(r,x,y,a)-Hr3(r,x,y,a)-1i*Hr4(r,x,y,a), Htb1(tb,x,y,a)+1i*Htb2(tb,x,y,a);...
       Htb1(tb,x,y,a)-1i*Htb2(tb,x,y,a), 1i*Hr1(r,x,y,a)+Hr2(r,x,y,a)-Hr3(r,x,y,a)+1i*Hr4(r,x,y,a), -Hsoc(soc,x,y,a)-Hpot(pot), 0;... 
       1i*Hr1(r,x,y,a)-Hr2(r,x,y,a)-Hr3(r,x,y,a)-1i*Hr4(r,x,y,a), Htb1(tb,x,y,a)-1i*Htb2(tb,x,y,a), 0, Hsoc(soc,x,y,a)-Hpot(pot)];
Ham = simplify(Ham);

e = eig(Ham);
e = simplify(e);

e1 = e(1); e2 = e(2); e3 = e(3); e4 = e(4);

[X,Y] = meshgrid(0:2*pi/50:2*pi);
e1 = subs(e1,x,X); e2 = subs(e2,x,X); e3 = subs(e3,x,X); e4 = subs(e4,x,X);


for i = 1:1:length(Y)
    for k = 1:1:length(Y)
        e1(i,k) = subs(e1(i,k),y,Y(i,k));
        e2(i,k) = subs(e2(i,k),y,Y(i,k));
        e3(i,k) = subs(e3(i,k),y,Y(i,k));
        e4(i,k) = subs(e4(i,k),y,Y(i,k));
    end
end
e1 = double(e1); e2 = double(e2); e3 = double(e3); e4 = double(e4);

surf(X,Y,e1);
surf(X,Y,e2);
surf(X,Y,e3);
surf(X,Y,e4);
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


















