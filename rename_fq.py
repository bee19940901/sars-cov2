from pathlib import Path

from utils import Cmds


class RenameFq:

    def __init__(self, args):
        self.cpu = int(args.cpu)
        self.in_file = Path(args.in_file).absolute()
        self.out_dir = Path(args.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)

    def rename_fq(self):
        cmds, shells = list(), list()
        with open(self.in_file, "r", encoding="utf-8") as fr:
            for line in fr.readlines()[1:]:
                if line := line.strip():
                    
                    sample_name, fq_path = line.split("\t")
                    
                    new_fq_path = self.out_dir.joinpath(f"{sample_name}.fq.gz")
                    cmds.append(f"cp {fq_path} {new_fq_path}\n")
                    shells.append(self.sh_dir.joinpath(f"{sample_name}.sh"))
        Cmds(cmds, shells).multi_run(self.cpu)


def main(args):
    RenameFq(args).rename_fq()
