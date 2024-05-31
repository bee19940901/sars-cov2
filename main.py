from pathlib import Path
from argparse import ArgumentParser

from utils import Cmds, Counter
from configs import PYTHON3, TOOLS, DATABASE, SOFTWARE, PIPELINE, PLATFORM2, PLATFORM3


class Pipeline:

    def __init__(self):

        parser = ArgumentParser(prog="xiehe-sars-cov2", description="协和新冠分析流程主脚本")
        parser.add_argument("-nsl", "--np_sample_list", type=str, required=True, help="三代样本与原始数据对应关系")
        parser.add_argument("-msl", "--mgi_sample_list", type=str, required=True, help="二代样本与原始数据对应关系")
        parser.add_argument("-od", "--out_dir", type=str, required=True, help="输出目录")
        parser.add_argument("-c", "--cpu", type=int, default=16)
        parser.add_argument("-mmd", "--mgi_min_depth", type=int, default=30)
        parser.add_argument("-nmd", "--np_min_depth", type=int, default=30)
        args = parser.parse_args()

        self.cpu = int(args.cpu)
        self.mgi_min_depth = int(args.mgi_min_depth)
        self.np_min_depth = int(args.np_min_depth)
        self.np_sample_list = Path(args.np_sample_list).absolute()
        self.mgi_sample_list = Path(args.mgi_sample_list).absolute()
        self.out_dir = Path(args.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.np_dir = self.out_dir.joinpath("nanopore")
        self.mgi_dir = self.out_dir.joinpath("mgi")
        self.oa_dir = self.out_dir.joinpath("overall")
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.np_dir.mkdir(parents=True, exist_ok=True)
        self.mgi_dir.mkdir(parents=True, exist_ok=True)
        self.oa_dir.mkdir(parents=True, exist_ok=True)
        self.counter = Counter()

    def mgi_rename(self):
        """
        二代数据拷贝并改名
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {TOOLS} rename_fq "
                f"-if {self.mgi_sample_list} "
                f"-od {self.mgi_dir.joinpath('rename_fqs')} "
                f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_rename.sh")
            ]
        ).multi_run(1)
        return self

    def np_rename(self):
        """
        三代数据拷贝并改名
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {TOOLS} rename_fq "
                f"-if {self.np_sample_list} "
                f"-od {self.np_dir.joinpath('rename_fqs')} "
                f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_rename.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_fastp(self):
        """
        二代数据质控
        :return:
        """
        Cmds(
            [
                 f"{PYTHON3} {TOOLS} fastp_qc "
                 f"-id {self.mgi_dir.joinpath('rename_fqs')} "
                 f"-od {self.mgi_dir.joinpath('clean_fqs')} "
                 f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_fastp.sh")
            ]
        ).multi_run(1)
        return self

    def np_porechop(self):
        """
        三代数据质控
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {TOOLS} porechop "
                f"-id {self.np_dir.joinpath('rename_fqs')} "
                f"-od {self.np_dir.joinpath('clean_fqs')} "
                f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_porechop.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_reads_stats(self):
        """
        二代 数据统计
        :return:
        """
        Cmds(
            [
                 f"{PYTHON3} {TOOLS} stats_reads "
                 f"-id {self.mgi_dir.joinpath('clean_fqs')} "
                 f"-od {self.mgi_dir.joinpath('reads_stats')} "
                 f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_reads_stats.sh")
            ]
        ).multi_run(1)
        return self

    def np_reads_stats(self):
        """
        三代 数据统计
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {TOOLS} stats_reads "
                f"-id {self.np_dir.joinpath('clean_fqs')} "
                f"-od {self.np_dir.joinpath('reads_stats')} "
                f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_reads_stats.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_bwa(self):
        """
        二代数据比对
        :return:
        """
        Cmds(
            [
                 f"{PYTHON3} {PIPELINE.joinpath('mgi_bwa.py')} "
                 f"-id {self.mgi_dir.joinpath('clean_fqs')} "
                 f"-od {self.mgi_dir.joinpath('bams')} "
                 f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_bwa.sh")
            ]
        ).multi_run(1)
        return self

    def np_minimap2(self):
        """
        三代数据比对
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('np_minimap2.py')} "
                f"-id {self.np_dir.joinpath('clean_fqs')} "
                f"-od {self.np_dir.joinpath('bams')} "
                f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_minimap2.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_trim(self):
        """
        二代数据修剪
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('mgi_trim.py')} "
                f"-id {self.mgi_dir.joinpath('bams')} "
                f"-od {self.mgi_dir.joinpath('trimmed_bams')} "
                f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_trim.sh")
            ]
        ).multi_run(1)
        return self

    def np_trim(self):
        """
        二代数据修剪
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('np_trim.py')} "
                f"-id {self.np_dir.joinpath('bams')} "
                f"-od {self.np_dir.joinpath('trimmed_bams')} "
                f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_trim.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_coverage(self):
        Cmds(
            [
                 f"{PYTHON3} {TOOLS} coverage "
                 f"-id {self.mgi_dir.joinpath('trimmed_bams')} "
                 f"-od {self.mgi_dir.joinpath('coverages')} "
                 f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_coverage.sh")
            ]
        ).multi_run(1)
        return self

    def np_coverage(self):
        Cmds(
            [
                 f"{PYTHON3} {TOOLS} coverage "
                 f"-id {self.np_dir.joinpath('trimmed_bams')} "
                 f"-od {self.np_dir.joinpath('coverages')} "
                 f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_coverage.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_depth(self):
        Cmds(
            [
                 f"{PYTHON3} {TOOLS} depth "
                 f"-id {self.mgi_dir.joinpath('trimmed_bams')} "
                 f"-od {self.mgi_dir.joinpath('depths')} "
                 f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_depth.sh")
            ]
        ).multi_run(1)
        return self

    def np_depth(self):
        Cmds(
            [
                 f"{PYTHON3} {TOOLS} depth "
                 f"-id {self.np_dir.joinpath('trimmed_bams')} "
                 f"-od {self.np_dir.joinpath('depths')} "
                 f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_depth.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_ivar(self):
        Cmds(
            [
                 f"{PYTHON3} {TOOLS} ivar_consensus "
                 f"-id {self.mgi_dir.joinpath('trimmed_bams')} "
                 f"-od {self.mgi_dir.joinpath('consensus')} "
                 f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_ivar.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_bcftools(self):
        Cmds(
            [
                 f"{PYTHON3} {TOOLS} bcftools_call "
                 f"-id {self.mgi_dir.joinpath('trimmed_bams')} "
                 f"-od {self.mgi_dir.joinpath('mutations')} "
                 f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_bcftools.sh")
            ]
        ).multi_run(1)
        return self

    def np_medaka(self):
        Cmds(
            [
                 f"{PYTHON3} {TOOLS} medaka "
                 f"-id {self.np_dir.joinpath('clean_fqs')} "
                 f"-od {self.np_dir.joinpath('consensus_mutations')} "
                 f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_medaka.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_nextclade(self):
        """
        评估二代数据基因组的质量
        :return:
        """
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'nextclade')} run "
                f"--input-dataset {DATABASE.joinpath('sars-cov-2')} "
                f"--output-all {self.mgi_dir.joinpath('nextclade')} "
                f"{self.mgi_dir.joinpath('consensus', 'all.consensus.merged.fas')} \n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_nextclade.sh")
            ]
        ).multi_run(1)
        return self

    def np_nextclade(self):
        """
        评估三代数据基因组质量
        :return:
        """
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'nextclade')} run "
                f"--input-dataset {DATABASE.joinpath('sars-cov-2')} "
                f"--output-all {self.np_dir.joinpath('nextclade')} "
                f"{self.np_dir.joinpath('consensus_mutations', 'all.consensus.merged.fas')} \n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_nextclade.sh")
            ]
        ).multi_run(1)
        return self

    def overall_stats(self):
        """
        对二三代测序数据进行一个整体的统计分析
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_stats.py')} "
                f"-mcf {self.mgi_dir.joinpath('coverages', 'all.coverage.merged.tsv')} "
                f"-ncf {self.np_dir.joinpath('coverages', 'all.coverage.merged.tsv')} "
                f"-msf {self.mgi_dir.joinpath('reads_stats', 'all.stats.merged.tsv')} "
                f"-nsf {self.np_dir.joinpath('reads_stats', 'all.stats.merged.tsv')} "
                f"-od {self.oa_dir.joinpath('stats')} \n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_stats.sh")
            ]
        ).multi_run(1)
        return self

    def overall_coverage(self):
        """
        对二三代测序数据的整体覆盖度进行分析
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_coverage.py')} "
                f"-mcf {self.mgi_dir.joinpath('coverages', 'all.coverage.merged.tsv')} "
                f"-ncf {self.np_dir.joinpath('coverages', 'all.coverage.merged.tsv')} "
                f"-od {self.oa_dir.joinpath('coverage')} \n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_coverage.sh")
            ]
        ).multi_run(1)
        return self

    def overall_depth(self):
        """
        对二三代测序数据的深度进行分析
        :return:
        """
        Cmds(
            [
                f"{PYTHON3}  {PIPELINE.joinpath('overall_depth.py')} "
                f"-mdd {self.mgi_dir.joinpath('depths')} "
                f"-ndd {self.np_dir.joinpath('depths')} "
                f"-od {self.oa_dir.joinpath('depth')} "
                f"-c {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_depth.sh")
            ]
        ).multi_run(1)
        return self

    def overall_consensus(self):
        """
        对二三代基因组进行评估
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_consensus.py')} "
                f"-mnf {self.mgi_dir.joinpath('nextclade', 'nextclade.csv')} "
                f"-nnf {self.np_dir.joinpath('nextclade', 'nextclade.csv')} "
                f"-od {self.oa_dir.joinpath('consensus')} "
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_consensus.sh")
            ]
        ).multi_run(1)
        return self

    def overall_ct(self):
        """
        CT值与numreads和coverage的关联分析
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_ct.py')} "
                f"-mcf {self.mgi_dir.joinpath('coverages', 'all.coverage.merged.tsv')} "
                f"-ncf {self.np_dir.joinpath('coverages', 'all.coverage.merged.tsv')} "
                f"-csf {DATABASE.joinpath('ct', 'ct.score.txt')} "
                f"-od {self.oa_dir.joinpath('ct')}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_ct.sh")
            ]
        ).multi_run(1)
        return self

    def overall_mutation(self):
        """
        变异位点数量的统计分析
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_mutation.py')} "
                f"-mvd {self.mgi_dir.joinpath('mutations')} "
                f"-nvd {self.np_dir.joinpath('consensus_mutations')} "
                f"-od {self.oa_dir.joinpath('mutation')}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_mutation.sh")
            ]
        ).multi_run(1)
        return self

    def overall_tree(self):
        """
        使用二三代的基因组构建系统发育树
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_tree.py')} "
                f"-mmf {self.mgi_dir.joinpath('consensus', 'all.consensus.merged.fas')} "
                f"-nmf {self.np_dir.joinpath('consensus_mutations', 'all.consensus.merged.fas')} "
                f"-nnc {self.np_dir.joinpath('nextclade', 'nextclade.csv')} "
                f"-mnc {self.mgi_dir.joinpath('nextclade', 'nextclade.csv')} "
                f"-c {self.cpu} "
                f"-od {self.oa_dir.joinpath('tree')}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_tree.sh")
            ]
        ).multi_run(1)
        return self

    def genes_coverage(self):
        """
        绘制基因的覆盖率和深度图
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('genes_coverage.py')} "
                f"-mdd {self.mgi_dir.joinpath('depths')} "
                f"-ndd {self.np_dir.joinpath('depths')} "
                f"-od {self.oa_dir.joinpath('genes_coverage')}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_genes_coverage.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_coverage_x(self):
        """
        计算覆盖深度>=min_depth的覆盖率
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('coverage_x.py')} "
                f"-cd {self.mgi_dir.joinpath('coverages')} "
                f"-dd {self.mgi_dir.joinpath('depths')} "
                f"-od {self.mgi_dir.joinpath('coverages_x')} "
                f"-md {self.mgi_min_depth}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_coverages_x.sh")
            ]
        ).multi_run(1)
        return self

    def np_coverage_x(self):
        """
        计算覆盖深度>=30的覆盖率
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('coverage_x.py')} "
                f"-cd {self.np_dir.joinpath('coverages')} "
                f"-dd {self.np_dir.joinpath('depths')} "
                f"-od {self.np_dir.joinpath('coverages_x')} "
                f"-md {self.np_min_depth}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_coverages_x.sh")
            ]
        ).multi_run(1)
        return self

    def overall_stats_x(self):
        """
        对二三代测序数据进行一个整体的统计分析
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_stats.py')} "
                f"-mcf {self.mgi_dir.joinpath('coverages_x', 'all.coverage.merged.xls')} "
                f"-ncf {self.np_dir.joinpath('coverages_x', 'all.coverage.merged.xls')} "
                f"-msf {self.mgi_dir.joinpath('reads_stats', 'all.stats.merged.tsv')} "
                f"-nsf {self.np_dir.joinpath('reads_stats', 'all.stats.merged.tsv')} "
                f"-od {self.oa_dir.joinpath('stats_x')} \n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_stats_x.sh")
            ]
        ).multi_run(1)
        return self

    def overall_coverage_x(self):
        """
        对二三代测序数据的整体覆盖度进行分析
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_coverage.py')} "
                f"-mcf {self.mgi_dir.joinpath('coverages_x', 'all.coverage.merged.xls')} "
                f"-ncf {self.np_dir.joinpath('coverages_x', 'all.coverage.merged.xls')} "
                f"-od {self.oa_dir.joinpath('coverage_x')} \n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_coverage_x.sh")
            ]
        ).multi_run(1)
        return self

    def overall_ct_x(self):
        """
        CT值与numreads和coverage的关联分析
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_ct.py')} "
                f"-mcf {self.mgi_dir.joinpath('coverages_x', 'all.coverage.merged.xls')} "
                f"-ncf {self.np_dir.joinpath('coverages_x', 'all.coverage.merged.xls')} "
                f"-csf {DATABASE.joinpath('ct', 'ct.score.txt')} "
                f"-od {self.oa_dir.joinpath('ct_x')}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_ct_x.sh")
            ]
        ).multi_run(1)
        return self

    def genes_coverage_x(self):
        """
        绘制基因的覆盖率和深度图
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('genes_coverage.py')} "
                f"-mdd {self.mgi_dir.joinpath('depths')} "
                f"-ndd {self.np_dir.joinpath('depths')} "
                f"-od {self.oa_dir.joinpath('genes_coverage_x')} "
                f"-mmd {self.mgi_min_depth} "
                f"-nmd {self.np_min_depth} "
                f"\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_genes_coverage_x.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_delly(self):
        """
        使用delly检测二代数据的结构体变异
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('mgi_delly.py')} "
                f"--in_bam_dir {self.mgi_dir.joinpath('trimmed_bams')} "
                f"--out_vcf_dir {self.mgi_dir.joinpath('delly')} "
                f"--cpu {self.cpu} "
                f"--platform {PLATFORM2}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_delly.sh")
            ]
        ).multi_run(1)
        return self

    def np_sniffles(self):
        """
        使用sniffles检测三代数据的结构体变异
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('np_sniffles.py')} "
                f"--in_fqs_dir {self.np_dir.joinpath('clean_fqs')} "
                f"--out_dir {self.np_dir.joinpath('sniffles')} "
                f"--cpu {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_sniffles.sh")
            ]
        ).multi_run(1)
        return self

    def overall_sv(self):
        """
        使用sniffles检测三代数据的结构体变异
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_sv.py')} "
                f"--in_np_file {self.np_dir.joinpath('sniffles', 'all.merged.xls')} "
                f"--in_mgi_file {self.mgi_dir.joinpath('delly', 'all.merged.xls')} "
                f"--out_sv_file {self.oa_dir.joinpath('sv', 'sv.xls')}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_sniffles.sh")
            ]
        ).multi_run(1)
        return self

    def np_consensus_vcf(self):
        """
        使用 medaka longshot 获取 三代数据的 consensus vcf 文件
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('np_consensus_vcf.py')} "
                f"--in_dir {self.np_dir.joinpath('trimmed_bams')} "
                f"--out_dir {self.np_dir.joinpath('consensus_vcfs')} "
                f"--cpu {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_np_consensus_vcf.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_consensus_vcf(self):
        """
        使用 bcftools 获取二代数据的vcf文件
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('mgi_consensus_vcf.py')} "
                f"--in_dir {self.mgi_dir.joinpath('trimmed_bams')} "
                f"--out_dir {self.mgi_dir.joinpath('consensus_vcfs')} "
                f"--cpu {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_mgi_consensus_vcf.sh")
            ]
        ).multi_run(1)
        return self

    def mgi_all_coverages(self):
        """
        统计二代数据所有的覆盖度信息
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('all_coverages.py')} "
                f"--in_dir {self.mgi_dir.joinpath('depths')} "
                f"--out_dir {self.mgi_dir.joinpath('all_coverages')} "
                f"--platform {PLATFORM2} "
                f"--cutoff 30 "
                f"--threads {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_{PLATFORM2}_all_coverages.sh")
            ]
        ).multi_run(1)
        return self

    def np_all_coverages(self):
        """
        统计二代数据所有的覆盖度信息
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('all_coverages.py')} "
                f"--in_dir {self.np_dir.joinpath('depths')} "
                f"--out_dir {self.np_dir.joinpath('all_coverages')} "
                f"--platform {PLATFORM3} "
                f"--cutoff 30 "
                f"--threads {self.cpu}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_{PLATFORM3}_all_coverages.sh")
            ]
        ).multi_run(1)
        return self

    def overall_consensus_vcf(self):
        """
        比较二三代变异位点结果
        :return:
        """
        Cmds(
            [
                f"{PYTHON3} {PIPELINE.joinpath('overall_consensus_vcf.py')} "
                f"--mgi_coverage_file {self.mgi_dir.joinpath('all_coverages', 'MGI.coverages.merged.xls')} "
                f"--np_coverage_file {self.np_dir.joinpath('all_coverages', 'QNome.coverages.merged.xls')} "
                f"--mgi_vcf_file {self.mgi_dir.joinpath('consensus_vcfs', 'all.merged.consensus.vcf')} "
                f"--np_vcf_file {self.np_dir.joinpath('consensus_vcfs', 'all.merged.consensus.vcf')} "
                f"--out_file {self.oa_dir.joinpath('consensus_vcf', 'result.xls')}\n"
            ],
            [
                self.sh_dir.joinpath(f"step{self.counter.step()}_overall_consensus_vcf.sh")
            ]
        ).multi_run(1)
        return self


if __name__ == "__main__":
    (
        Pipeline()
        .mgi_rename()
        .mgi_fastp()
        .mgi_reads_stats()
        .mgi_bwa()
        .mgi_trim()
        .mgi_coverage()
        .mgi_depth()
        .mgi_ivar()
        .mgi_bcftools()
        .mgi_nextclade()
        .np_rename()
        .np_porechop()
        .np_reads_stats()
        .np_minimap2()
        .np_trim()
        .np_coverage()
        .np_depth()
        .np_medaka()
        .np_nextclade()
        .overall_stats()
        .overall_coverage()
        .overall_depth()
        .overall_consensus()
        .overall_ct()
        .overall_mutation()
        .overall_tree()
        .genes_coverage()
        .mgi_coverage_x()
        .np_coverage_x()
        .overall_stats_x()
        .overall_coverage_x()
        .overall_ct_x()
        .genes_coverage_x()
        .np_sniffles()
        .mgi_delly()
        .np_consensus_vcf()
        .mgi_consensus_vcf()
        .mgi_all_coverages()
        .np_all_coverages()
        .overall_consensus_vcf()
    )
