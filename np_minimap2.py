from argparse import ArgumentParser
from pathlib import Path
from re import sub

from utils import Cmds, Results
from configs import SOFTWARE, DATABASE


MINIMAP2 = SOFTWARE.joinpath("bin", "minimap2")
SAMTOOLS = SOFTWARE.joinpath("bin", "samtools")
REF = DATABASE.joinpath('cons_ref', 'sars-cov2.ref.fa')


class NpMinimap2:

    def __init__(self):
        ap = ArgumentParser()
        ap.add_argument("-id", "--in_dir", type=str, required=True)
        ap.add_argument("-od", "--out_dir", type=str, required=True)
        ap.add_argument("-c", "--cpu", type=int, default=8)
        ag = ap.parse_args()
        self.cpu = ag.cpu
        self.in_dir = Path(ag.in_dir).absolute()
        self.out_dir = Path(ag.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.clean_fqs = list(self.in_dir.glob("*.fq.gz"))
        self.samples = [sub(r"\.fq\.gz$", "", i.name) for i in self.clean_fqs]
        self.bams = [self.out_dir.joinpath(f'{i}.bam') for i in self.samples]

    def minimap2(self):
        Cmds(
            [
                f"{MINIMAP2} -ax map-ont -t 8 "
                f"{REF} {clean_fq} | "
                f"{SAMTOOLS} view -@ 8 -h -F 4 > {bam}"
                for clean_fq, bam in zip(self.clean_fqs, self.bams)
            ],
            [
                self.sh_dir.joinpath(f'{i}.sh') for i in self.samples
            ]
        ).multi_run(self.cpu)
        Results(self.bams).check_exists().check_empty()


if __name__ == "__main__":
    NpMinimap2().minimap2()
