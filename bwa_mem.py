from pathlib import Path
import re

from configs import SOFTWARE, DATABASE
from utils import Cmds


class BwaMem:

    def __init__(self, args):
        """
        :param args: 命令行参数
        """
        self.bed_file = Path(args.bed_file).absolute() if args.bed_file else ""
        self.in_dir = Path(args.in_dir).absolute()
        self.cpu = int(args.cpu)
        self.out_dir = Path(args.out_dir).absolute()

        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)

        self.clean_fqs = list(self.in_dir.glob("*.fq.gz"))
        self.samples = [re.sub(r"\.fq\.gz$", "", i.name) for i in self.clean_fqs]
        self.sorted_bams = [self.out_dir.joinpath(f"{i}.sorted.bam") for i in self.samples]

    def bwa_mem(self):
        """
        将reads比对到参考基因组上
        :return:
        """
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'bwa')} mem -t 8 "
                f"{DATABASE.joinpath('cons_ref', 'sars-cov2.ref.fa')} {clean_fq} | "
                f"{SOFTWARE.joinpath('bin', 'ivar')} trim -b {DATABASE.joinpath('primers', 'mgiseq.primers.bed')} "
                f"-x 3 -m 30 | "
                f"{SOFTWARE.joinpath('bin', 'samtools')} view -@ 8 -h -F 4 | "
                f"{SOFTWARE.joinpath('bin', 'samtools')} sort -@ 8 -o {sorted_bam} && "
                f"{SOFTWARE.joinpath('bin', 'samtools')} index {sorted_bam} "
                if
                self.bed_file
                else
                f"{SOFTWARE.joinpath('bin', 'bwa')} mem -t 8 "
                f"{DATABASE.joinpath('cons_ref', 'sars-cov2.ref.fa')} {clean_fq} | "
                f"{SOFTWARE.joinpath('bin', 'samtools')} view -@ 8 -h -F 4 | "
                f"{SOFTWARE.joinpath('bin', 'samtools')} sort -@ 8 -o {sorted_bam} && "
                f"{SOFTWARE.joinpath('bin', 'samtools')} index {sorted_bam} "
                for clean_fq, sorted_bam in zip(self.clean_fqs, self.sorted_bams)
            ],
            [
                self.sh_dir.joinpath(f"{sample}.sh")
                for sample in self.samples
            ]
        ).multi_run(self.cpu)
        return self


def main(args):
    """
    :param args: 命令行参数
    :return:
    """
    BwaMem(args).bwa_mem()
