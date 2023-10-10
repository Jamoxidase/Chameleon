# Stealth

Stealth looks for under-represented Kmers in a genome file by using a statistical model to establish  an expected count of a Kmer and comparing that to what is actually there.

## Background
### Null Model derivation
If we assume complete independence of the any of the characters that make up a sequence K, we could approximate the expected count of K, E(K),  as:  
- Pr(K)=E(K)/N = Pr(k1)Pr(k2)Pr(k3)Pr(k4), and
- Pr(k) = c(k)/N

so:
- E(K) = N * Pr(k1)Pr(k2)Pr(k3)Pr(k4)

the above is a markov(0) approximation for the expected count (E) of our sequence K.

Rarely is this assumption of independence a good null model, due to many factors including codon bias which drives selection of 3-mers based on underlying amino acid frequencies. We can model this dependence using a markov model, where the probability associated with any symbol is conditioned on the preceding characters. For a kmer of size 4, we can condition our character probabilities using the preceding two bases without requiring use of the count of the specific 4-mer that we are trying to approximate(see note on information content). We need to consider edge conditions in our approximation. In general, we can use a model of order k-2, where k is the size of the kmer as long as we have sufficient data, therefore:

- Pr(K)=E(K)/N = Pr(k1)Pr(k2|k1)Pr(k3|k1k2)Pr(k4|k2k3) = c(k1)/N * c(k1k2)/c(k1) * c(k1k2k3)/c(k1k2) * c(k2k3k4)/c(k2k3)
 after cancelling, we get this:
- Pr(K) = 1/N * c(k1k2k3) * c(k2k3k4) / c(k2k3)

and we need a way to estimate parameters of some distribution so that we can estimate a probability of an actual observation of c(K).

If we think of this as a binomial, or a poisson, we can estimate  a mean and variance. 
- mu = N * P, which in our case is E(K), the expected count of K
- var = N * P * (1-P), which simplifies to E(K) * (1 - E(K)/N)


because we have so much data, and the probabilities are so small, we could just estimate, var = mu, which is the Poisson variance
we can then build Z-scores from these data and actual counts by reshaping to a standard normal by subtracting the mean from our observation and dividing my the SD ( sqrt(var)).
So:
- Z = (c(K) - E(K)) / sqrt(var)

In short, count all of the kmers of interest and all of their component parts down to size 1
For this version of Stealth, I deal with odd sized Kmers by assuming that the middle character is N (dont care), and then summing all counts that are part of the composite sequence, so:
- c(AGNCT) = sum(c(AGACT, c(AGGCT), c(AGCCT), c(AGTCT)))

## Installing Stealth
You need Python3, with scipy installed. Stealth will use sys, math, itertools and scipy.stats

## Running Stealth
A Fasta file containing your genome of interest is read through STDIN.
Options:
- -c --cutoff Zscore ( this needs to be a negative float. default -4.0)
- -p --pseudocount ( an integer that is added to every possible kMer. default 0)
- -P --Palindrome ( flag telling stealth to only print reverse-complement palindromes)
- -m --max (the largest Kmer size to test. Keep this at 8 or below, especially for small genomes. Default is 8)

from command line, an example would be:
- python3 stealth.v0 -P -c -6.0 -m 6 < Hp-SS1.fa >Hp-SS1.out

If there are multiple .fa records in the input file, all of them will be evaluated as a single (not concatenated) genome.

## Multiple Hypothesis testing
Stealth.v0 does not adjust for multiple hypothesis testing. If, for example, you want a cutoff of -6.0, you will need to calculate that p-value from a Z-table, divide it by 96, then 
you get a Z-score cutoff of -6.7. With just palindromic sequences, 96 comes from the 16 tests of 4mers, 16 tests of degenerate 5mers, and the 64 possible 6mers. Doing this automatically 
would be a really nice future add.

## Roadmap
This is the first Stealth release. The next release is nearly complete and will include a FDR estimation method for multiple hypothesis testing. It will also be based on Chi-square test of independence,
and will provide better visuals.

## Contributing
I am open to contributions.
Standard test cases are needed.

## Authors and acknowledgment
David Bernick is the original author.

## License
This work is protected with a UC Santa Cruz Noncommercial License. Please see the License document contained in this repository.

## Project status
Still in active development as needed.

## Modifications

```text
Stealth.V0.py
- Stealth.V0.py filename changed -> _stealth.py
- __all__ defined

__init__.py
- Created for the Stealth package module
- Imports everything from _stealth.py as defined by __all__

_stealthClass.py
- StealthV0 class to store Stealth analysis as a Python Object
```

