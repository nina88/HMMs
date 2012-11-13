HMMs
====

# Why use HMMs?

`HMMs` fits hidden Markov models to amino acid sequences. The amino acid sequences are recoded according to charge and hydrophobicity. `HMMs` uses the power posterior method to determine how many segment types to look for. Gibbs sampling can then be used to determine where these segment types are. 

# Installing the package

The package can be installed from Github using the devtools R package.
  
    library(devtools)
    install_github("HMMs", username = "nina88")

# Using the package

You can load proteins into R  that are in FASTA format saved in a txt file for use in this package 

    library(HMMs)
    x1=read_FASTA("p53.txt",f)
    x2=read_FASTA("mdm2.txt",f)
    x3=read_FASTA("CBP.txt",f)
    y=combine_FASTA(x1,x2)
    y=combine_FASTA(y,x3)

You can run a power posterior analysis my initialising a prior for P and lambda and a checkpointing file name, then running the analysis.

    prior=initialise_prior(P.mat, mu, s, r, f)
    cp=initialise_checkpoint_power("cp.Rdata")
    powerpost(N, prior, m, r, y, burnin, checkpoint = cp)

Once you have a chosen value for r then you can use Gibbs sampling to find out where the segments types are and what the transition structures are.

    cp=initialise_checkpoint("cp.Rdata",hour)
    gibbs(y, iter, prior, r, burnin, thin, checkpoint = cp)


# Thanks




-----------

