"""
使用fastp对原始数据进行质控
"""
from pathlib import Path
import re

from utils import Cmds, Results
from configs import  SOFTWARE


class FastpQC:

    def __init__(self, args):
        """
        :param args: 命令行参数
        """
        self.cpu = args.cpu
        self.in_dir = Path(args.in_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.rename_fqs = list(self.in_dir.glob("*.fq.gz"))
        self.samples = [re.sub(r"\.fq\.gz$", "", i.name) for i in self.rename_fqs]
        self.out_jsons = [self.out_dir.joinpath(f"{i}.json") for i in self.samples]
        self.out_htmls = [self.out_dir.joinpath(f"{i}.html") for i in self.samples]
        self.clean_fqs = [self.out_dir.joinpath(f"{i}.fq.gz") for i in self.samples]

    def fastp(self):
        """
        使用fastp进行质控
        :return:
        """
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'fastp')} -w 8 -i {rename_fq} -o {clean_fq} -h {out_html} -j {out_json} "
                for rename_fq, clean_fq, out_html, out_json
                in zip(self.rename_fqs, self.clean_fqs, self.out_htmls, self.out_jsons)
            ],
            [
                self.sh_dir.joinpath(f"{sample}.sh")
                for sample in self.samples
            ]
        ).multi_run(self.cpu)
        Results(self.clean_fqs + self.out_jsons + self.out_htmls).check_exists().check_empty()
        return self


def main(args):
    """
    :param args: 命令行参数
    :return:
    """
    (
        FastpQC(args)
        .fastp()
    )
