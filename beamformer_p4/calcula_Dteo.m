function [Dteo] = calcula_Dteo(w, f, d, Vprop, barrido)
    
    N = 7;
    D = zeros(1,length(barrido));

    for i=1:length(barrido)
        Dsum = 0;
        for n=1:N
            Dsum = Dsum + exp(2 * pi * 1j * f * n * d * cos(barrido(i)) / Vprop) * conj(w(f, n));
    
        D(i) = Dsum;
        end
    end
    
    Dsim = fliplr(D);
    Dteo = abs([D Dsim]);
    
end

