"""
将reads比对到参考基因组上
获取consensus sequence
"""
from pathlib import Path

from Bio import SeqIO
import re


from utils import Cmds, Results
from configs import SOFTWARE


class IvarConsensus:

    def __init__(self, args):
        """
        :param args: 命令行参数
        """
        self.cpu = int(args.cpu)
        self.in_dir = Path(args.in_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()

        self.out_file = self.out_dir.joinpath("all.consensus.merged.fas")
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)

        self.sorted_bams = list(self.in_dir.glob("*.sorted.bam"))
        self.samples = [re.sub(r"\.sorted\.bam$", "", i.name) for i in self.sorted_bams]
        self.out_pres = [self.out_dir.joinpath(i) for i in self.samples]
        self.out_fas = [self.out_dir.joinpath(f"{i}.fa") for i in self.samples]
        self.merged_fas = self.out_dir.joinpath("all.consensus.merged.fas")

    def ivar_consensus(self):
        """
        使用ivar获取组装的基因组
        :return:
        """
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'samtools')} mpileup -aa -A -Q 0 -d 0 {sorted_bam} | "
                f"{SOFTWARE.joinpath('bin', 'ivar')} consensus -p {out_pre} -m 10 -n N -t 0.5\n"
                for sorted_bam, out_pre in zip(self.sorted_bams, self.out_pres)
            ],
            [
                self.sh_dir.joinpath(f"{sample}.sh")
                for sample in self.samples
            ]
        ).multi_run(self.cpu)
        Results(self.out_fas).check_exists().check_empty()
        return self

    def merge(self):
        with open(self.merged_fas, "w", encoding="utf-8") as fw:
            for s, f in zip(self.samples, self.out_fas):
                with open(f, "r", encoding="utf-8") as fr:
                    for record in SeqIO.parse(fr, "fasta"):
                        fw.write(f">{s}\n{record.seq}\n")
        return self


def main(args):
    """
    :param args: 命令行参数
    :return:
    """
    (
        IvarConsensus(args)
        .ivar_consensus()
        .merge()
    )
