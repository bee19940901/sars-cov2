from pathlib import Path


SCRIPT = Path(__file__).absolute()
PIPELINE = SCRIPT.parent
PROJECT = PIPELINE.parent
DATABASE = PROJECT.joinpath("database")
SOFTWARE = PROJECT.joinpath("software")
WORKSPACE = PROJECT.joinpath("workspace")
PYTHON3 = SOFTWARE.joinpath('bin', 'python3')
TOOLS = PIPELINE.joinpath("tools.py")
RSCRIPT = SOFTWARE.joinpath('bin', 'Rscript')
PLATFORM2 = "MGI"
PLATFORM3 = "QNome"

NGMLR = SOFTWARE.joinpath('bin', 'ngmlr')
IVAR = SOFTWARE.joinpath('bin', 'ivar')
SAMTOOLS109 = SOFTWARE.joinpath('bin', 'samtools109')
SAMTOOLS = SOFTWARE.joinpath('bin', 'samtools')
SNIFFLES = SOFTWARE.joinpath('bin', 'sniffles')
DELLY = SOFTWARE.joinpath('bin', 'delly')
BCFTOOLS = SOFTWARE.joinpath('bin', 'bcftools')
MINIMAP2 = SOFTWARE.joinpath('bin', 'minimap2')
BWA = SOFTWARE.joinpath('bin', 'bwa')

HOMO_FA = DATABASE.joinpath("hg38", "ref.fa")
REF_FA = DATABASE.joinpath("cons_ref", "sars-cov2.ref.fa")
REF_LEN = 29903
PRIMERS_BED = DATABASE.joinpath("primers", "mgiseq.primers.bed")
BLACKLIST = DATABASE.joinpath("cons_ref", "LCS", "sars-cov2.Low-Complexity.txt")