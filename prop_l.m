%parameters
lamda = 1.064 * 10.^(-3); %wavelength of light in mm 
w = 1 ; %beam waist in mm
k = 2 * pi / lamda; % wave number
N = 500;
L = 5; %in mm

%Rayleigh Range
zR = pi*w^2/lamda;

%Along z-axis
d = 3000;

%z-value list
Ez = [];

%range
x = linspace(-L,L,N);
y = linspace(-L,L,N);
z_values = linspace(-d,d,N);
[X, Y] = meshgrid(x, y);


%E(x,y,0)
E = exp(-(X.^2 + Y.^2) / w^2);

%E(kx,ky,0)
E_k = fftshift(fft2(E)); 

%Defining wavevectors
dx = x(2) - x(1);
dy = y(2) - y(1);
dkx = 2 * pi /(N*dx);
dky = 2 * pi /(N*dy);

kx = (-N/2:N/2-1)*dkx;
ky = (-N/2:N/2-1)*dky;
[KX,KY] = meshgrid(kx,ky);

for i = 1:length(z_values) 
    z = z_values(i);
    E_kz = E_k .* exp(-1j*z*(KX.^2+KY.^2)/(2*k)) .* exp(1j * k * z);
    E_z = ifft2(ifftshift(E_kz));
    newrow = E_z(N/2,:);
    Ez = [Ez; newrow];
end 

% Plotting the results
figure;
imagesc(z_values, x, abs(Ez').^2);
colormap turbo
colorbar;
xlabel('z (mm)');
ylabel('x (mm)');
title('Electric Field Distribution in xz-plane');

