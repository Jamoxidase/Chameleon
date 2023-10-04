#!/usr/bin/env python3
import itertools, math
import scipy.stats as stats
import sys


class FastAreader:
    """Read fasta file from specified fname or STDIN if not given"""

    def __init__(self, fname=""):
        """contructor: saves attribute fname"""

        self.fname = fname
        self.fileH = None

    def doOpen(self):
        if self.fname == "":
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        header = ""
        sequence = ""

        with self.doOpen() as self.fileH:
            header = ""
            sequence = ""

            # skip to first fasta header
            line = self.fileH.readline()
            while not line.startswith(">"):
                line = self.fileH.readline()
            header = line[1:].rstrip()

            for line in self.fileH:
                if line.startswith(">"):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ""
                else:
                    sequence += "".join(line.rstrip().split()).upper()

        yield header, sequence


########################################################################
# CommandLine
########################################################################


class CommandLine:
    """
    Handle the command line, usage and help requests.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    """

    def __init__(self, inOpts=None):
        """
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        """

        import argparse

        self.parser = argparse.ArgumentParser(
            description="Program prolog - read a genome and return under-represented kMers",
            epilog="Thank you",
            add_help=True,  # default is True
            prefix_chars="-",
            usage="%(prog)s [options] -option1[default] <input >output",
        )

        self.parser.add_argument(
            "-m",
            "--max",
            nargs="?",
            default=8,
            type=int,
            action="store",
            help="max kMer size ",
        )
        self.parser.add_argument(
            "-c",
            "--cutoffProb",
            nargs="?",
            type=float,
            default=-4.0,
            action="store",
            #
            help="lowest Zscore allowed for printing ",
        )

        self.parser.add_argument(
            "-p",
            "--pseudoCount",
            type=int,
            default=1,
            action="store",
            help="pseudocount to be added to each count",
        )
        self.parser.add_argument(
            "-P",
            "--Palindrome",
            type=bool,
            nargs="?",
            const=True,
            default=False,
            action="store",
            help="set if only palindromes print",
        )

        self.parser.add_argument(
            "-v", "--version", action="version", version="%(prog)s 0.1"
        )
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


def rc(s):
    """Return the reverse complement of the sequence."""
    return s[::-1].translate(
        str.maketrans("ACGTN", "TGCAN")
    )  # N is needed for composites


class Genome:
    bases = "ACGT"

    def __init__(self, min=1, max=8, pseudo=1):
        # initialize all kMers to the pseudocount value
        self.counts = {}
        for k in range(min, max + 1):  # kMer size
            newKcounts = {
                "".join(kmer): pseudo
                for kmer in itertools.product(Genome.bases, repeat=k)
            }
            self.counts.update(newKcounts)

        self.min, self.max, self.pseudo = min, max, pseudo
        self.n = sum(self.counts.values())  # add the pseudos to initialize n

    def addSequence(self, seq):
        """Count Kmers in the genome, and build composite sums that will be needs for odd kMers."""
        for l in range(self.min, self.max + 1):
            for p in range(len(seq) - self.max):
                try:
                    self.counts[seq[p : p + l]] += 1
                except KeyError:
                    pass
        self.n += len(seq) - self.max
        # build the composite sequences for all of the odd ones and their components, then put them in counts
        composites = {}
        for s in self.counts:
            if (len(s) > 2) and (
                len(s) % 2 == 1
            ):  # it's odd, so build a composite sum by adding over all center bases.
                flank = len(s) // 2
                newSeq = s[:flank] + "N" + s[flank + 1 :]  # build the composite
                composites[newSeq] = sum(
                    (
                        self.counts[s[:flank] + b + s[(flank + 1) :]]
                        for b in Genome.bases
                    )
                )
                # build left and right k-1mers and central k-2mer
                composites[newSeq[:-1]] = sum(
                    (
                        self.counts[s[:flank] + b + s[(flank + 1) : -1]]
                        for b in Genome.bases
                    )
                )
                composites[newSeq[1:]] = sum(
                    (
                        self.counts[s[1:flank] + b + s[(flank + 1) :]]
                        for b in Genome.bases
                    )
                )
                composites[newSeq[1:-1]] = sum(
                    (
                        self.counts[s[1:flank] + b + s[(flank + 1) : -1]]
                        for b in Genome.bases
                    )
                )
        self.counts.update(composites)

    def E(self, s):
        """return expected count of s. P(ab) = P(a|b) * P(b). P(abc) = P(a|b) * P(b|c) * P(c)"""
        # P(abcd) = P(a) * P(b|a) * P(c|ab) * P(d|bc)
        # E(abcd) = c(a) * c(ab)/c(a) * c(abc)/c(ab) * c(bcd)/c(bc) *len(s) = c(abc) * c(bcd) / c(bc)
        # for odd length s, assumes that the composites have already been built

        try:
            return self.counts[s[:-1]] * self.counts[s[1:]] / self.counts[s[1:-1]]
        except KeyError:
            return self.counts[s]
        except ZeroDivisionError:
            return 0.0

    def pValue(self, s):
        return stats.binom.cdf(self.counts[s], self.n, self.E(s) / self.n)

    def Zscore(self, s):
        """calculate Zscore
        s = kMer sequence
        """
        mu = self.E(s)
        var = mu * (1 - mu / self.n)  # using binomial, estimate variance
        if var == 0.0:
            return 1000000000
        return (self.counts[s] - mu) / math.sqrt(var)

    def Evalue(self, s):
        return self.pValue(s) * self.n


########################################################################
# Main
# Here is the main program
#
#
########################################################################


def main(myCommandLine=None):
    """
    Report underRepresented kMers in a genome.

    STDIN: FastA file of genome
    Options:
        -c cutoff value in Zscores. Sequences with lower Zscore probs will be reported [-4.0]
        -p pseudoCount value if desired [0]
        -P boolean flag specifying that only RC Palindromes should be reported [False]
        -m maximum kMer size to analyze.[8]

    :param myCommandLine:
    :return:
    """
    cl = CommandLine()

    # get the genome
    sourceReader = FastAreader()

    # setup the Tetramer object using centerSequences as the reference
    thisGenome = Genome(1, cl.args.max, cl.args.pseudoCount)
    molecules = 0
    for head, seq in sourceReader.readFasta():
        thisGenome.addSequence(seq)
        molecules += 1
    print("N = {}, molecules = {}".format(thisGenome.n, molecules))

    outArray = []
    for seq, count in thisGenome.counts.items():
        # for odd len sequence, we just print the composite that is for the full sequence.
        if "N" in seq:  # this is a composite
            p = seq.find("N")
            if len(seq) != (p * 2 + 1):
                continue  # skip this one, it is a part of a composite
        elif len(seq) % 2 == 1:
            continue  # odd sequences print only if they are a composite

        outArray.append(
            [
                seq,
                count,
                thisGenome.E(seq),
                thisGenome.Evalue(seq),
                thisGenome.Zscore(seq),
            ]
        )
    outArray.sort(key=lambda e: (len(e[0]), -e[4]), reverse=True)
    for seq, count, E, eVal, Z in outArray:
        if (seq in thisGenome.counts) and (Z <= cl.args.cutoffProb):
            # if the -P flag was set, only print palindromes, otherwise print palindromes with a flag
            rcFlag = ""
            if rc(seq) != seq:  # this must not be a palindrome
                if cl.args.Palindrome:
                    continue  # dont print non palindromes
            elif (
                not cl.args.Palindrome
            ):  # it is a palindrome so if we are printing everything give it a flag
                rcFlag = "RC palindrome"

            print(
                "{0:8}\t{1:0d}\t{2:0.2f}\t{3:0.2e}\t{4:0.2f}\t{5:13}".format(
                    seq, count, E, eVal, Z, rcFlag
                )
            )


__all__ = ["Genome", "rc"]

if __name__ == "__main__":
    main()
