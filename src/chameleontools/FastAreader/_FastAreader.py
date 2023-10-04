from typing import Iterator
import sys


class FastAreader:
    """
    Define objects to read FastA files.

    Args:
    fname (str): file name (optional), default is None (STDIN)

    Usage:
    thisReader = FastAreader ('testTiny.fa')
    for head, seq in thisReader.readFasta():
        print (head,seq)
    """

    def __init__(self, fname=None):
        """
        Constructor: saves attribute fname
        Args:
        fname (str): file name (optional), default is None (STDIN)
        """
        self.fname = fname

    def doOpen(self):
        """
        Handle file opens, allowing STDIN.

        Returns:
        file: either sys.stdin or file handler for the input file
        """
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self) -> Iterator[tuple[str, str]]:
        """
        Read an entire FastA record and return the sequence header/sequence

        Yields:
        tuple: header and sequence as a tuple
        """
        header = ""
        sequence = ""

        with self.doOpen() as fileH:
            header = ""
            sequence = ""

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith(">"):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith(">"):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ""
                else:
                    sequence += line.strip()
        yield header, sequence


__all__ = ["FastAreader"]
