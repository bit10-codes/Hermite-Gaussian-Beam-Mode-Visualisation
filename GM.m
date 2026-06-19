%parameters
lamda = 1.064 * 10.^(-3); %wavelength of light in mm 
w = 1 ; %beam waist in mm
k = 2 * pi / lamda; % wave number
N = 500;
L = 5; %in mm

%Check field and intensity at
lst = [0, 1500, 3500, 5000];

%range
x = linspace(-L,L,N);
y = linspace(-L,L,N);
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

for i = 1:length(lst)
    z = lst(i);
    E_kz = E_k .* exp(-1j*z*(KX.^2+KY.^2)/(2*k)) .* exp(1j * k * z);
    E_z = ifft2(ifftshift(E_kz));

    %Surface Plot
    figure;
    surf(X,Y, abs(E_z))
    shading interp
    colormap turbo
    zlim([0, 1])
    title(['Electric Field Profile at z = ', num2str(z/1000), ' m'])
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Electric Field Amplitude');
end


figure
hold on
for m = 1:length(lst)
    zl = lst(m);
    E_kzl = E_k .* exp(-1j*zl*(KX.^2+KY.^2)/(2*k)) .* exp(1j * k * zl);
    E_zl = ifft2(E_kzl);

    plot(x,abs(E_zl(N/2,:)), 'DisplayName',['z = ',num2str(zl/1000), 'm'])
    legend
end
title('Beam Profile at different z-planes')
xlabel('X (mm)');
ylabel('Electric Field Amplitude');