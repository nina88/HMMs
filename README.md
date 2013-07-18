HMMs
====

### Why use HMMs?

`HMMs` fits hidden Markov models to amino acid sequences. The amino acid sequences are recoded according to charge and hydrophobicity. `HMMs` uses the power posterior method to determine how many segment types to look for. Gibbs sampling can then be used to determine where these segment types are. 

### Installing the package

The package can be installed directly from Github using the devtools R package.
  
    library(devtools)
    install_github("HMMs", username = "nina88")

### Using the package

You can load proteins into R  that are in FASTA format saved in a txt file using (`read_FASTA`) for use in this package or you can use data from the package as shown below. The data must be of class `hmm_FASTA` which is done by the `read_FASTA` oe `read_FASTA_string` function. These functions also convert the amino acid sequence into another sequence which depends on the properties of the amino acids. This means you can recode 20 amino acids into either 2, 3 or 4 items. They are either converted according to hydrophobicity (f=2), charge (f=3), or both (f=4). You can combine proteins in order to provide more information for the analysis using the `combine_FASTA`. This process is shown below. 

    library(HMMs)
    f = 4
    data(TAF15)
    data(FUS)
    x1 = read_FASTA_string(TAF15,f)
    x2 = read_FASTA_string(FUS,f)
    y = combine_FASTA(x1,x2)

You can run a power posterior analysis by firstly initialising the prior distributions. P.mat sets all the parameters of the Dirichlet distributions to be 1 in the priors for P. (Each row of P is given a dirichlet prior). This gives us a prior mean of changing between amino acids types to be 1/f for all segment types r. We have more prior information about lambda as we expect segment types to be fairly long so we set the prior mean of staying in the same segment type (mu) to be quite high and the standard deviation of staying in the same segment type to be quite low (s). We use the `initialise_prior` command to initialise the prior distributions for P and lambda. This command also requires the numer of segment types (r) and the number of amino acid types after recoding (f). Next if you want the code to checkpoint you need to initialise checkpointing by choosing a file name. This must be an `.Rdata` file. Finally you can run the power posterior analysis using the `powerpost` function. The powerposterior analysis has a temperature parameter which are values between 0 and 1, N is used to determine how many values for the temperature parameter are used. m is the number of iterations that are used per temperature parameter. y is the `hmm_FASTA` object. You must choose a burnin and you can checkpoint if you have initialised checkpointing or you can choose not to checkpoint by just missing out the command.  You can repeat the analysis for a range of values of r (the number of segment types), for example r=2,..,10. This analysis gives you the log marginal likelihood for r. You can then use Bayes theorem to calculate the marginal posterior distribution for r.

    r = 3
    P.mat = array(1, c(f, f, r))
    prior = initialise_prior(P.mat, mu=0.99, s=0.01)
    cp = initialise_checkpoint_power("cp.Rdata")
    powerpost(N=40, prior, m=10000, y, burnin=2000, checkpoint = cp)
    

Once you have a chosen value for r from your marginal posterior distribution for r you can then use Gibbs sampling to find out where the segments types are and what the transition structures are. You must initialise checkpointing and specify after how many iterations you would like checkpointing to happen (hour). You must again specify a burnin and how often you would line thinning to occur. iter is the number of iterations.

    cp = initialise_checkpoint("cp.Rdata", hour=500)
    gibbs(y, iter=10000, prior, burnin=2000, thin=20, checkpoint = cp)

