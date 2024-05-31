from pathlib import Path
import re

from utils import Results, Cmds
from configs import SOFTWARE, DATABASE


class SelectReads:

    def __init__(self, args):
        self.in_dir = Path(args.in_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.cpu = args.cpu

    def select_reads(self):
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        in_fqs = list(self.in_dir.glob("*.fq.gz"))
        samples = [re.sub(r"\.fq\.gz$", "", i.name) for i in in_fqs]
        mapped_fqs = [self.out_dir.joinpath(f"{sample}.mapped.fq") for sample in samples]
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'bwa')} mem -t 32 {DATABASE.joinpath('ncbi_dataset','ref.fa')} {in_fq} | "
                f"{SOFTWARE.joinpath('bin', 'samtools')} fastq -@ 32 -F 4 > {mapped_fq} "
                for in_fq, mapped_fq in zip(in_fqs, mapped_fqs)
            ],
            [
                self.sh_dir.joinpath(f"{sample}.sh")
                for sample in samples
            ]
        ).multi_run(self.cpu)
        Results(mapped_fqs).check_exists().check_empty()


def main(args):
    """
    主函数
    :param args: 命令行参数
    :return:
    """
    SelectReads(args).select_reads()
