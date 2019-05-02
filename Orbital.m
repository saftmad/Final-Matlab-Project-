
close all

% quantum numbers, can make whatever i want to change graph
n = 2; % principle quantum number 
l = 1;  % orbital/sublevl
m = 0;  % magnetic orientation/spin 

% plotting parameters, heisenberg uncertainty principle
probabilitydensity = 1e-5;
a = 1; % radius 

% angular part (Condon-Shortley phase) for angular movement in quantum
SphericalYlm = @(l, m, theta, phi) (-1)^m * sqrt((2 * l + 1) / (4 * pi) * factorial(l - abs(m)) / factorial(l + abs(m))) * AssociatedLegendre(l, m, cos(theta)) .* exp(1i * m * phi);
% give spherical harmonics in coordinates later 
% use associated legendre polynomial 
if (m < 0)
    Y = @(l, m, theta, phi) sqrt(2) * (-1)^m * imag(SphericalYlm(l, abs(m), theta, phi));
elseif (m == 0)
    Y = @(l, m, theta, phi) SphericalYlm(l, m, theta, phi);
else
    Y = @(l, m, theta, phi) sqrt(2) * (-1)^m * real(SphericalYlm(l, m, theta, phi));
end

% radial part function to find the acutal given radius 
R = @(n, l, r) sqrt((2 / (a * n))^3 * factorial(n - l - 1) / (2 * n * factorial(n + l))) .* exp(-r / (a * n)) .* (2 * r / (a * n)).^l * 1 / factorial(n - l - 1 + 2 * l + 1) .* AssociatedLaguerre(n - l - 1, 2 * l + 1, 2 * r / (a * n));

% wave function, differntial with respect to gamma 
psi = @(n, l, m, r, theta, phi) R(n, l, r) .* Y(l, m, theta, phi);

% setting the grid
border = 32;
accuracy = 100;
raster = linspace(-border, border, accuracy);
[x, y, z] = ndgrid(raster, raster, raster);

% converting rectangular to spherical 
r = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z ./ r);
phi = atan2(y, x);

% plot orbital,  - and + wave function phase
f = sign(psi(n, l, m, r, theta, phi));
isosurface(psi(n, l, m, r, theta, phi).^2, probabilitydensity,f);
colormap([0 0 1; 1 0.5 0])

% functions
function Anm = AssociatedLaguerre(n,m,x)
Anm = 0;
    for i = 0 : n
        Anm = Anm + factorial(m + n) * nchoosek(m + n, n - i) / factorial(i) * (-x).^i;
    end
end

function Alm = AssociatedLegendre(l,m,x)
Alm = 0;
    for r = 0 : floor(1/2 * l - 1/2 * abs(m))
        Alm = Alm + (-1)^r * nchoosek(l - 2 * r, abs(m)) * nchoosek(l, r) * ...
            nchoosek(2 * l - 2 * r, l) * x.^(l - 2 * r - abs(m));
    end
    Alm = (1 - x.^2).^(abs(m) / 2) .* (factorial(abs(m)) / 2^l * Alm);
end


