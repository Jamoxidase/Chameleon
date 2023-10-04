from ._stealth import Genome, rc as reverse_complement
from chameleontools.SeqParser import StealthGenome
from Bio.Seq import Seq


class StealthV0:
    def __init__(
        self,
        seq: Seq | str | list[str | Seq] | StealthGenome,
        zScore_cut: float,
        pseudo: int,
        kMax: int,
        kMin: int = 1,
    ) -> None:
        if type(seq) == StealthGenome:
            seq = seq.getGenome()
        self.output = self._runStealth(seq, zScore_cut, pseudo, kMax, kMin)

    def _runStealth(
        self, genome: str, zScore_cut: float, pseudo: int, kMax: int, kMin: int = 1
    ) -> list[tuple[str, float, bool]]:
        """
        runs Stealth analysis on input genome
        saves output to array -> [(motif,zScore)]
        """
        if type(genome) != list:
            genome = [genome]
        stealth = Genome(1, kMax, pseudo)
        for seq in genome:
            stealth.addSequence(str(seq))  # conversion of BioPython Seq() to str
        """Taken from original stealth.vo.py main() func"""
        ######V Written by David L. Bernick V#######
        outArray = []
        for seq, _ in stealth.counts.items():
            if len(seq) < kMin:
                continue
            # for odd  len sequence, we just print the composite that is for the full sequence.
            if "N" in seq:  # this is a composite
                p = seq.find("N")
                if len(seq) != (p * 2 + 1):
                    continue  # skip this one, it is a part of a composite
            elif len(seq) % 2 == 1:
                continue  # odd sequences print only if they are a composite
            ######^ Written by David L. Bernick ^#######
            zscore = stealth.Zscore(seq)
            if zscore <= zScore_cut:
                if seq in stealth.counts:
                    rc_flag = False
                    if reverse_complement(seq):
                        rc_flag = True
                    outArray.append((seq, zscore, rc_flag))
        return outArray

    def getOutput(self) -> list[tuple[str, float, bool]]:
        return self.output


__all__ = ["StealthV0"]
