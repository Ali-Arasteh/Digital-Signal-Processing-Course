%%
% 1.2.1
n = -10:100;
Coefficent_Of_Input = [1 0.5];
Coefficent_Of_Output = [1 -1.8*cos(pi/8) 0.81];
Input = zeros(1,111);
Input(11) = 1;
h = filter(Coefficent_Of_Input,Coefficent_Of_Output,Input);
stem(n,h)
xlabel('n');
ylabel('h[n]');
title('Impulse response');
%%
% 1.2.2
n = -10:100;
h = ((0.5 - 1.933i)*(0.8315 + 0.3444i) .^ n + (0.5 + 1.933i)*(0.8315 - 0.3444i) .^ n) .* heaviside(n);
h(11) = h(11) + 1/2;
stem(n,h)
xlabel('n');
ylabel('h[n]');
title('Impulse response');
%%
% 1.3 1)
Phase_Of_R = roots([1 0 0.9]);
P1 = Phase_Of_R(1);
P2 = Phase_Of_R(2);
h0 = -0.033;
h1 = 0.6;
A = [1 1;P1 P2];
B = [h0 h1];
Coefficent = A\B';
Alpha = Coefficent(1);
Beta = Coefficent(2);
n = -10:100;
h = (Alpha * P1.^n + Beta * P2.^n) .* heaviside(n);
h(11) = h(11) + 1/3 + 1/2 * -0.033;
stem(n,h)
xlabel('n');
ylabel('h[n]');
title('Impulse response');
%%
% 1.3 2)
Phase_Of_R = roots([1 -1.8*cos(pi/8) 0.81]);
P1 = Phase_Of_R(1);
P2 = Phase_Of_R(2);
h0 = 1;
h1 = 0.5 + 1.8*cos(pi/8);
A = [1 1;P1 P2];
B = [h0 h1];
Coefficent = A\B';
Alpha = Coefficent(1);
Beta = Coefficent(2);
n = -10:100;
h = (Alpha * P1.^n + Beta * P2.^n) .* heaviside(n);
h(11) = h(11) + 1/2;
stem(n,h)
xlabel('n');
ylabel('h[n]');
title('Impulse response');
%%
% 2.1
%%
% 2.2.1 1)
L = 12;
n = 1:10*L;
r = ones(1,L);
R = DTFT(r,10*L);
Real_Part_Of_R = real(R);
Imaginary_Part_Of_R = imag(R);
Magnitude_Of_R = abs(R);
Omega = n*(2*pi/length(n)) - pi;
subplot(3,1,1)
plot(Omega,Real_Part_Of_R);
xlabel('Frequency');
ylabel('Real Part Of R');
title('Real part of fourier transform of r[n] with L = 12');
subplot(3,1,2)
plot(Omega,Imaginary_Part_Of_R);
xlabel('Frequency');
ylabel('Imaginary Part Of R');
title('Imaginary part of fourier transform of r[n] with L = 12');
subplot(3,1,3)
plot(Omega,Magnitude_Of_R);      
xlabel('Frequency');
ylabel('Magnitude Of R');
title('Magnitude of fourier transform of r[n] with L = 12');
%%
% 2.2.1 2)
L = 15;
n = 1:10*L;
r = ones(1,L);
R = DTFT(r,10*L);
Real_Part_Of_R = real(R);
Imaginary_Part_Of_R = imag(R);
Magnitude_Of_R = abs(R);
Omega = n*(2*pi/length(n)) - pi;
subplot(3,1,1)
plot(Omega,Real_Part_Of_R);
xlabel('Frequency');
ylabel('Real Part Of R');
title('Real part of fourier transform of r[n] with L = 15');
subplot(3,1,2)
plot(Omega,Imaginary_Part_Of_R);
xlabel('Frequency');
ylabel('Imaginary Part Of R');
title('Imaginary part of fourier transform of r[n] with L = 15');
subplot(3,1,3)
plot(Omega,Magnitude_Of_R);      
xlabel('Frequency');
ylabel('Magnitude Of R');
title('Magnitude of fourier transform of r[n] with L = 15');
%%
% 2.2 1)
L = 12;
n = 1:10*L;
r = ones(1,L);
R = DTFT(r,10*L);
Phase_Of_R = angle(R);
Unwrap_Phase_Of_R = unwrap(Phase_Of_R);
Omega = n*(2*pi/length(n)) - pi;
subplot(2,1,1)
plot(Omega,Phase_Of_R)
xlabel('Frequency');
ylabel('Phase Of R');
title('Phase of the fourier transform of r[n] with L = 12');
subplot(2,1,2)
plot(Omega,Unwrap_Phase_Of_R)
xlabel('Frequency');
ylabel('Unwrap Phase Of R');
title('Phase of the fourier transform of r[n] with L = 12');
%%
% 2.2 2)
L = 15;
n = 1:10*L;
r = ones(1,L);
R = DTFT(r,10*L);
Phase_Of_R = angle(R);
Unwrap_Phase_Of_R = unwrap(Phase_Of_R);
Omega = n*(2*pi/length(n)) - pi;
subplot(2,1,1)
plot(Omega,Phase_Of_R)
xlabel('Frequency');
ylabel('Phase Of R');
title('Phase of the fourier transform of r[n] with L = 15');
subplot(2,1,2)
plot(Omega,Unwrap_Phase_Of_R)
xlabel('Frequency');
ylabel('Unwrap Phase Of R');
title('Phase of the fourier transform of r[n] with L = 15');
%%
% 2.3
n = -499:500;
x1 = (sin(pi*n/10).^2)./((pi*n/10).^2);
x1(500) = 1;
x2 = (sin(pi*n/10))./((pi*n/10));
x2(500) = 1;
y1 = (sin(pi*n/5))./((pi*n/5));
y1(500) = 1;
y2 = zeros(1,1000);
y2(2:2:1000) = x2(251:1:750);
y3 = x2 .* sin(2*pi*0.3*n);
X1 = DTFT(x1,length(x1));
X2 = DTFT(x2,length(x2));
Y1 = DTFT(y1,length(y1));
Y2 = DTFT(y2,length(y2));
Y3 = DTFT(y3,length(y3));
Magnitude_Of_X1 = abs(X1);
Magnitude_Of_X2 = abs(X2);
Magnitude_Of_Y1 = abs(Y1);
Magnitude_Of_Y2 = abs(Y2);
Magnitude_Of_Y3 = abs(Y3);
Omega = n*(pi/500);
figure();
plot(Omega,Magnitude_Of_X1)
xlabel('Frequency');
ylabel('Magnitude Of X1');
title('Magnitude of fourier transform of x1[n]');
figure();
plot(Omega,Magnitude_Of_X2)
xlabel('Frequency');
ylabel('Magnitude Of X2');
title('Magnitude of fourier transform of x2[n]');
figure();
plot(Omega,Magnitude_Of_Y1)
xlabel('Frequency');
ylabel('Magnitude Of Y1');
title('Magnitude of fourier transform of y1[n]');
figure();
plot(Omega,Magnitude_Of_Y2)
xlabel('Frequency');
ylabel('Magnitude Of Y2');
title('Magnitude of fourier transform of y2[n]');
figure();
plot(Omega,Magnitude_Of_Y3)
xlabel('Frequency');
ylabel('Magnitude Of Y3');
title('Magnitude of fourier transform of y3[n]');
%%
function Output = DTFT(Input,n)
Output = fft(Input,n);
Output = fftshift(Output);
end
