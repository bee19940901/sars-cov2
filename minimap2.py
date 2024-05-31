from pathlib import Path
import re

from utils import Cmds
from configs import SOFTWARE, DATABASE


class Minimap2:

    def __init__(self,args):

        self.cpu = int(args.cpu)
        self.in_dir = Path(args.in_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()

        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True ,exist_ok=True)

        self.clean_fqs = list(self.in_dir.glob("*.fq.gz"))
        self.samples = [re.sub(r"\.fq\.gz$", "", i.name) for i in self.clean_fqs]
        self.bams = [self.out_dir.joinpath(f'{i}.sorted.bam') for i in self.samples]

    def minimap2(self):
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'minimap2')} -ax map-ont -t 8 "
                f"{DATABASE.joinpath('cons_ref', 'sars-cov2.ref.fa')} {clean_fq} | "
                f"{SOFTWARE.joinpath('bin', 'samtools')} view -@ 8 -h -F 4 | "
                f"{SOFTWARE.joinpath('bin', 'samtools')} sort -@ 8 -o {bam} && "
                f"{SOFTWARE.joinpath('bin', 'samtools')} index {bam}\n"
                for clean_fq, bam in zip(self.clean_fqs, self.bams)
            ],
            [
                self.sh_dir.joinpath(f'{i}.sh') for i in self.samples
            ]
        ).multi_run(self.cpu)


def main(args):
    Minimap2(args).minimap2()
