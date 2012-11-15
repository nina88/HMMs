HMMs
====

# Why use HMMs?

`HMMs` fits hidden Markov models to amino acid sequences. The amino acid sequences are recoded according to charge and hydrophobicity. `HMMs` uses the power posterior method to determine how many segment types to look for. Gibbs sampling can then be used to determine where these segment types are. 

# Installing the package

The package can be installed from Github using the devtools R package.
  
    library(devtools)
    install_github("HMMs", username = "nina88")

# Using the package

You can load proteins into R  that are in FASTA format saved in a txt file using (read_FASTA) for use in this package or you can use data from the package as follows 

    library(HMMs)
    f=4
    data(TAF15)
    data(FUS)
    data(EWS)
    x1=read_FASTA_string(TAF15,f)
    x2=read_FASTA_string(FUS,f)
    x3=read_FASTA_string(EWS,f)
    y=combine_FASTA(x1,x2)
    y=combine_FASTA(y,x3)

You can run a power posterior analysis my initialising a prior for P and lambda and a checkpointing file name, then running the analysis.

    r=3
    P.mat=array(1,c(f,f,r))
    prior=initialise_prior(P.mat, mu=0.99, s=0.01, r, f)
    cp=initialise_checkpoint_power("cp.Rdata")
    powerpost(N=40, prior, m=10000, r, y, burnin=2000, checkpoint = cp)
    

Once you have a chosen value for r then you can use Gibbs sampling to find out where the segments types are and what the transition structures are.

    cp=initialise_checkpoint("cp.Rdata",hour=500)
    gibbs(y, iter=10000, prior, r, burnin=2000, thin=20, checkpoint = cp)


# Thanks




-----------

