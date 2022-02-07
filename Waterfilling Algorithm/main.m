function [] = main ()

p.Nt = 3;
p.Nr = 3;
p.iterations = 1000;

SNRdB = -10 : 1 : 20;
SNR = 10.^(SNRdB/10);

openloop = zeros(p.iterations, length(SNRdB));
closedloop = zeros(p.iterations, length(SNRdB));

for idx = 1 : length(SNRdB)
    for jdx = 1 : p.iterations
        H = sqrt(1/2) * (randn(p.Nr, p.Nt) + 1i * randn(p.Nr, p.Nt));
        [~, SIGMA, ~] = svd(H);
        
        sigma = diag(SIGMA);
        noiselevel = zeros(1, rank(H));
        for kdx = 1 : rank(H)
            noiselevel(kdx) = p.Nt / (SNR(idx) .* sigma(kdx)^2);
        end
        
        wf = waterfilling (p, noiselevel);
        
        for ldx = 1 : rank(H)
            openloop(jdx, idx) = openloop(jdx, idx) + log(1 + (SNR(idx)/p.Nt) * sigma(ldx).^2 * wf(ldx)) / log(2);
            closedloop(jdx, idx) = closedloop(jdx, idx) + log(1 + SNR(idx)/p.Nt * sigma(ldx).^2) / log(2);
        end
        
    end
end

ergodiccapacity = mean(openloop);
spectraleff = mean(closedloop);

plot(SNRdB, ergodiccapacity, '-*', SNRdB, spectraleff, '--*', 'LineWidth', 1.5);
title('3 \times 3 MIMO System');
xlabel('SNR (dB)');
ylabel('Capacity / Spectral Efficiency');
legend('Closed Loop', 'Open Loop');
grid on

end