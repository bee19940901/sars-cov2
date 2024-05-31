"""
获取样本的测序深度
"""
import re
from pathlib import Path

from configs import SOFTWARE
from utils import Cmds, Results


class Depth:

    def __init__(self, args):
        self.cpu = int(args.cpu)
        self.in_dir = Path(args.in_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()

        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)

        self.in_bams = list(self.in_dir.glob("*.sorted.bam"))
        self.samples = [re.sub(r"\.sorted\.bam$", "", i.name) for i in self.in_bams]
        self.out_depths = [self.out_dir.joinpath(f"{s}.depth") for s in self.samples]

    def depth(self):
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'samtools')} depth -@ 8 -aa {in_bam} > {out_depth}\n"
                for in_bam, out_depth in zip(self.in_bams, self.out_depths)
            ],
            [
                self.sh_dir.joinpath(f"{s}.sh") for s in self.samples
            ]
        ).multi_run(self.cpu)
        Results(self.out_depths)
        return self


def main(args):
    Depth(args).depth()
