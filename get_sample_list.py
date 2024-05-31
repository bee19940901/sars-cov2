from pathlib import Path
import re


class GetSampleList:

    def __init__(self, args):
        self.in_dir = Path(args.in_dir).absolute()
        self.raw_fqs = list(self.in_dir.glob("*.fq.gz"))
        self.raw_sample_names = [re.sub(r"\.fq\.gz$", "", i.name) for i in self.raw_fqs]
        self.samples = [re.search(r"^XG\d+-\d+", i).group(0) for i in self.raw_sample_names]
        self.out_file = Path(args.out_file)
        self.out_dir = self.out_file.parent
        self.out_dir.mkdir(parents=True, exist_ok=True)

    def get_sample_list(self):
        with open(self.out_file, "w", encoding="utf-8") as fw:
            fw.write("sample_name\tfastq_path\n")
            for raw_fq, sample in zip(self.raw_fqs, self.samples):
                fw.write(f"{sample}\t{raw_fq}\n")


def main(args):
    GetSampleList(args).get_sample_list()
