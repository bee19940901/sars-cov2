from pathlib import Path
import re
from Bio import SeqIO

from utils import Cmds
from configs import SOFTWARE, DATABASE


class Medaka:

    def __init__(self, args):

        self.cpu = int(args.cpu)
        self.in_dir = Path(args.in_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()

        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)

        self.merged_fas = self.out_dir.joinpath("all.consensus.merged.fas")
        self.in_fqs = list(self.in_dir.glob("*.fq.gz"))
        self.samples = [re.sub(r"\.fq\.gz$", "", i.name) for i in self.in_fqs]
        self.out_dirs = [self.out_dir.joinpath(i) for i in self.samples]
        self.out_consensus_fas = [i.joinpath("consensus.fasta") for i in self.out_dirs]
        self.out_hdfs = [i.joinpath("consensus_probs.hdf") for i in self.out_dirs]
        self.out_bams = [i.joinpath("calls_to_draft.bam") for i in self.out_dirs]
        self.out_raw_vcfs = [i.joinpath("consensus.vcf") for i in self.out_dirs]
        self.raw_vcfs = [self.out_dir.joinpath(f"{s}.raw.vcf") for s in self.samples]
        self.good_vcfs = [self.out_dir.joinpath(f"{s}.good.vcf") for s in self.samples]

    def medaka(self):
        Cmds(
            [
                f"export PATH={SOFTWARE.joinpath('bin')}:$PATH && "
                f"medaka_consensus -i {in_fq} "
                f"-d {DATABASE.joinpath('cons_ref', 'sars-cov2.ref.fa')} "
                f"-t 8 -m r941_min_hac_variant_g507 -o {out_dir} && "
                f"medaka variant {DATABASE.joinpath('cons_ref', 'sars-cov2.ref.fa')} {out_hdf} {out_raw_vcf} && "
                f"medaka tools annotate {out_raw_vcf} {DATABASE.joinpath('cons_ref', 'sars-cov2.ref.fa')} "
                f"{out_bam} {raw_vcf} && "
                f"bcftools view -i 'QUAL >= 10 && INFO/DP >= 10' {raw_vcf} > {good_vcf}\n"
                for in_fq, out_dir, out_hdf, out_raw_vcf, out_bam, raw_vcf, good_vcf
                in zip(self.in_fqs, self.out_dirs, self.out_hdfs, self.out_raw_vcfs,
                       self.out_bams, self.raw_vcfs, self.good_vcfs)
            ],
            [
                self.sh_dir.joinpath(f"{s}.sh") for s in self.samples
            ]
        ).multi_run(self.cpu)
        return self

    def merge(self):
        with open(self.merged_fas, "w", encoding="utf-8") as fw:
            for s, f in zip(self.samples, self.out_consensus_fas):
                with open(f, "r", encoding="utf-8") as fr:
                    for record in SeqIO.parse(fr, "fasta"):
                        fw.write(f">{s}\n{record.seq}\n")
        return self


def main(args):
    Medaka(args).medaka().merge()
