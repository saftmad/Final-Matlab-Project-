function PlotHydrogenMolecularOrbital
while(1) 
    clear; 
  
close all 

n = 4 
l = 1;  
m = 0;
 
 n = input('Enter the principal quantum number (shell #) n = '); 
 l = input ('Enter the orbital quantum number l = '); 
 while l > (n-1) 
     fprintf('l must be less than n \n'); 
     l = input('Enter the orbital quantum number l = ');
 end 
 
 m = input('Enter the magnetic quantum number m = '); 
 while m > l
     fprintf('m must be less than or equal to l \n'); 
     m = input('Enter the magnetic quantum number m = '); 
 end 
 
 
 
end 


probabilitydensity = 1e-5;
a = 1;  

R = @(n, l, r) sqrt((2 / (a * n))^3 * factorial(n - l - 1) / (2 * n * factorial(n + l))) .* exp(-r / (a * n)) .* (2 * r / (a * n)).^l * 1 / factorial(n - l - 1 + 2 * l + 1) .* AssociatedLaguerre(n - l - 1, 2 * l + 1, 2 * r / (a * n));

psi = @(n, l, m, r, theta, phi) R(n, l, r) .* Y(l, m, theta, phi);


raster = linspace(-30, 30, 100);
[x, y, z] = ndgrid(raster, raster, raster);

r = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z ./ r);
phi = atan2(y, x);

colors = sign(psi(n, l, m, r, theta, phi));
isosurface(psi(n, l, m, r, theta, phi).^2, probabilitydensity, colors);
colormap([0 0 1; 1 0.5 0])




