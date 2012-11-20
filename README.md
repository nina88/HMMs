HMMs
====

### Why use HMMs?

`HMMs` fits hidden Markov models to amino acid sequences. The amino acid sequences are recoded according to charge and hydrophobicity. `HMMs` uses the power posterior method to determine how many segment types to look for. Gibbs sampling can then be used to determine where these segment types are. 

### Installing the package

The package can be installed directly from Github using the devtools R package.
  
    library(devtools)
    install_github("HMMs", username = "nina88")

### Using the package

You can load proteins into R  that are in FASTA format saved in a txt file using (`read_FASTA`) for use in this package or you can use data from the package as shown below. The data must be of class `hmm_FASTA` which is done by the `read_FASTA` oe `read_FASTA_string` function. These functions also convert the amino acid sequence into another sequence which depends on the properties of the amino acids. They are either converted according to hydrophobicity (f=2), charge (f=3), or both (f=4). You can combine proteins in order to provide more information for the analysis using the `combine_FASTA`. This process is shown below. 

    library(HMMs)
    f = 4
    data(TAF15)
    data(FUS)
    data(EWS)
    x1 = read_FASTA_string(TAF15,f)
    x2 = read_FASTA_string(FUS,f)
    x3 = read_FASTA_string(EWS,f)
    y = combine_FASTA(x1,x2)
    y = combine_FASTA(y,x3)

You can run a power posterior analysis my initialising a prior for P and lambda and a checkpointing file name, then running the analysis. You can repeat the analysis for a range of values of r (the number of segment types), for example r=2,..,10. This analysis gives you the log marginal likelihood for r. You can then use Bayes theorem to calculate the marginal posterior distribution for r.

    r = 3
    P.mat = array(1, c(f, f, r))
    prior = initialise_prior(P.mat, mu=0.99, s=0.01, r, f)
    cp = initialise_checkpoint_power("cp.Rdata")
    powerpost(N=40, prior, m=10000, r, y, burnin=2000, checkpoint = cp)
    

Once you have a chosen value for r from your marginal posterior distribution for r you can then use Gibbs sampling to find out where the segments types are and what the transition structures are.

    cp = initialise_checkpoint("cp.Rdata", hour=500)
    gibbs(y, iter=10000, prior, r, burnin=2000, thin=20, checkpoint = cp)

