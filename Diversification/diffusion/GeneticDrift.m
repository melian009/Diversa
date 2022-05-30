#Octave-Matlab-code
for r = 1:10;
    N = 100;count = 0;G = 100;
    for i = 2:G;
        P(1,1) = 0.6;
        count = count + 1;
        X = binornd(N,P(count,1));#Return a matrix of random samples from the binomial distribution with parameters N and P, where N is the number of trials and P is the probability of success.
        P(i,1) = X/N;
    end
    n = 1:G;
    for j = 1:N;
           hold on
           t1 =unifrnd(0,1);
           t2 =unifrnd(0,1);
           t3 =unifrnd(0,1);
           h = plot(n,P,'k');
           set(h,'color',[t1 t2 t3]);
    hold on
    #plot(n,P,'k')
    end
end

