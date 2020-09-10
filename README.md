# BayesianStatistics
### Model of cell-to-cell interaction for analysis of tumor pathology images
**Team members**: Diana Isaeva, Alexandra Pershakova, Maria Vittoria Trussoni



## Abstract
The digital pathology imaging of tumor tissue produces massive data that can capture histological details in high resolution at a large scale. It enables detailed identification, analysis and classification up to an individual cell level thanks to the recent developments in computational techniques.

Our study focuses on the problem of modeling \textbf{spatial correlations} among three cells commonly observed in tumor pathology images: \textit{lymphocyte}, \textit{stromal}, and \textit{tumor}, and approaches the problem in a Bayesian framework, with the aim of modeling how the mark in a pattern might have been formed given the points, with interpretable underlying parameters.

The main goal of the study is to apply Markov Chain Monte Carlo (MCMC) sampling techniques, combined with the double Metropolis-Hastings (DMH) algorithm, to sample from the posterior distribution with an intractable normalizing constant. The additional aim is to study the statistical significance of the estimated parameters for the prediction of the patient survival time.
