from argparse import ArgumentParser
from pathlib import Path
import re

from configs import SOFTWARE, DATABASE
from utils import Cmds, Results


REF_FA = DATABASE.joinpath('cons_ref', 'sars-cov2.ref.fa')
BWA = SOFTWARE.joinpath('bin', 'bwa')
SAMTOOLS = SOFTWARE.joinpath('bin', 'samtools')


class MgiBwa:

    def __init__(self):

        ap = ArgumentParser()
        ap.add_argument("-id", "--in_dir", type=str, required=True)
        ap.add_argument("-od", "--out_dir", type=str, required=True)
        ap.add_argument("-c", "--cpu", type=int, default=8)
        ag = ap.parse_args()

        self.in_dir = Path(ag.in_dir).absolute()
        self.cpu = int(ag.cpu)
        self.out_dir = Path(ag.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.clean_fqs = list(self.in_dir.glob("*.fq.gz"))
        self.samples = [re.sub(r"\.fq\.gz$", "", i.name) for i in self.clean_fqs]
        self.out_bams = [self.out_dir.joinpath(f"{i}.bam") for i in self.samples]

    def bwa(self):
        """
        将reads比对到参考基因组上
        :return:
        """
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        Cmds(
            [
                f"{BWA} mem -t 8 {REF_FA} {clean_fq} | "
                f"{SAMTOOLS} view -@ 8 -h -F 4 > {out_bam}"
                for clean_fq, out_bam in zip(self.clean_fqs, self.out_bams)
            ],
            [
                self.sh_dir.joinpath(f"{sample}.sh")
                for sample in self.samples
            ]
        ).multi_run(self.cpu)
        Results(self.out_bams).check_exists().check_empty()
        return self


if __name__ == "__main__":
    MgiBwa().bwa()
