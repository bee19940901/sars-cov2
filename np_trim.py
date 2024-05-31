from argparse import ArgumentParser
from pathlib import Path
from re import sub

from configs import SOFTWARE, DATABASE
from utils import Cmds, Results


SAMTOOLS = SOFTWARE.joinpath("bin", "samtools")
IVAR = SOFTWARE.joinpath("bin", "ivar")
BED = DATABASE.joinpath('primers', 'mgiseq.primers.bed')


class NpTrim:

    def __init__(self):

        ap = ArgumentParser()
        ap.add_argument("-id", "--in_dir", type=str, required=True)
        ap.add_argument("-od", "--out_dir", type=str, required=True)
        ap.add_argument("-c", "--cpu", type=int, default=8)
        ag = ap.parse_args()

        self.cpu = int(ag.cpu)
        self.in_dir = Path(ag.in_dir).absolute()
        self.out_dir = Path(ag.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.in_bams = list(
            self.in_dir.glob("*.bam")
        )
        self.samples = [
            sub(r"\.bam$", "", i.name)
            for i in self.in_bams
        ]
        self.out_bams = [
            self.out_dir.joinpath(f"{i}.sorted.bam")
            for i in self.samples
        ]
        self.out_trimmed_fqs = [
            self.out_dir.joinpath(f"{i}.trimmed.fq")
            for i in self.samples
        ]

    def trim(self):
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        Cmds(
            [
                f"{IVAR} trim -b {BED} -i {in_bam} -e -q 0 | "
                f"{SAMTOOLS} sort -@ 8 -o {out_bam} && "
                f"{SAMTOOLS} index {out_bam} && "
                f"{SAMTOOLS} fastq -@ 8 {out_bam} > {out_trimmed_fq} \n "
                for in_bam, out_bam, out_trimmed_fq
                in zip(self.in_bams, self.out_bams, self.out_trimmed_fqs)
            ],
            [
                self.sh_dir.joinpath(f"{i}.sh")
                for i in self.samples
            ]
        ).multi_run(self.cpu)
        Results(self.out_bams).check_empty().check_exists()
        return self


if __name__ == "__main__":
    (
        NpTrim()
        .trim()
    )
