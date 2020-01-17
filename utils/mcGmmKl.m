function [kldiv] = mcGmmKl(gmm1, gmm2, n)
%MCGMMKL Computes the KL-divgergence between two gaussian mixture models [KL(gmm1||gmm2)] through a MonteCarlo approximation
%gmm1, gmm2: two gmdistribution objects
%n: number of monte carlo samples

%this does manual sampling
%chosen = roulette( gmm1.PComponents, n);
%mu = gmm1.mu(chosen,:);
%sigma = gmm1.Sigma(:,:,chosen);
%sample = mvnrnd(mu, sigma, n);
%
%let matlab generate the sample, just fixing the seed
s = RandStream.create('mt19937ar','seed',5489);
RandStream.setDefaultStream(s);
sample = random(gmm1, n);

denGMM1 = pdf(gmm1, sample);
denGMM2 = pdf(gmm2, sample);
kldiv = mean(log(denGMM1 ./ denGMM2));

end
