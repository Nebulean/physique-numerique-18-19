d = load("g_fourier.out");

t = d(:,1).*1e2;
theta = d(:,2);
thetadot = d(:,3);

fourier_theta = abs(fft(theta));
fourier_thetadot = abs(fft(thetadot));

Fs = length(theta); % sampling
L = t(end); % longueur du signal

freq_theta = Fs*(0:(L/2))/L;
freq_thetadot = Fs*(0:(L/2))/L;

% On vire le dernier fourier
fourier_theta(end) = [];
fourier_thetadot(end) = [];

f=figure;
hold on;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 18);

plot(freq_theta, fourier_theta);
plot(freq_thetadot, fourier_thetadot);

xlabel("frequency [Hz]");
ylabel("fft($\theta$, $\dot{\theta}$)");

lgd = legend("$\theta$", "$\dot{\theta}$", 'Location', 'northwest');

grid on;

set(gca, 'XScale', 'log');
hold off;

saveas(f, 'graphs/g_fourier_non_chaotic','epsc')
