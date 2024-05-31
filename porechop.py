from pathlib import Path
import re

from utils import Cmds, Results
from configs import SOFTWARE


class Porechop:

    def __init__(self, args):

        self.cpu = int(args.cpu)
        self.in_dir = Path(args.in_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()

        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)

        self.in_fqs = list(self.in_dir.glob("*.fq.gz"))
        self.samples = [re.sub(r"\.fq\.gz$", "", i.name) for i in self.in_fqs]
        self.clean_fqs = [self.out_dir.joinpath(f"{i}.fq.gz") for i in self.samples]

    def porechop(self):
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'porechop')} -t 8 --format fastq.gz -i {in_fq} -o {clean_fq}\n"
                for in_fq, clean_fq in zip(self.in_fqs, self.clean_fqs)
            ],
            [
                self.sh_dir.joinpath(f"{i}.sh") for i in self.samples
            ]
        ).multi_run(self.cpu)
        Results(self.clean_fqs)
        return self


def main(args):
    Porechop(args).porechop()
