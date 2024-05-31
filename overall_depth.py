"""
为每个深度文件绘制图片
"""
import re
from pathlib import Path
from collections import defaultdict

import pandas as pd

from configs import PIPELINE, SOFTWARE, PLATFORM3, PLATFORM2
from utils import Cmds, Results


class OverallDepth:

    def __init__(self, args):

        self.cpu = int(args.cpu)
        self.mgi_depth_dir = Path(args.mgi_depth_dir).absolute()
        self.np_depth_dir = Path(args.np_depth_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()

        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.mgi_plots_dir = self.out_dir.joinpath(f"{PLATFORM2}_plots")
        self.np_plots_dir = self.out_dir.joinpath(f"{PLATFORM3}_plots")
        for i in [self.sh_dir, self.mgi_plots_dir, self.np_plots_dir]:
            i.mkdir(parents=True, exist_ok=True)

        self.cv_file = self.out_dir.joinpath("cv.xls")
        self.cv_pdf = self.out_dir.joinpath("cv.pdf")
        self.cv_png = self.out_dir.joinpath("cv.png")

        self.mgi_depths = list(self.mgi_depth_dir.glob("*.depth"))
        self.mgi_samples = [re.sub(r"\.depth$", "", i.name) for i in self.mgi_depths]
        self.mgi_depth_pdfs = [self.mgi_plots_dir.joinpath(f"{s}.pdf") for s in self.mgi_samples]
        self.mgi_depth_pngs = [self.mgi_plots_dir.joinpath(f"{s}.png") for s in self.mgi_samples]

        self.np_depths = list(self.np_depth_dir.glob("*.depth"))
        self.np_samples = [re.sub(r"\.depth$", "", i.name) for i in self.np_depths]
        self.np_depth_pdfs = [self.np_plots_dir.joinpath(f"{s}.pdf") for s in self.np_samples]
        self.np_depth_pngs = [self.np_plots_dir.joinpath(f"{s}.png") for s in self.np_samples]

    def get_cv_file(self):
        cv_dict = defaultdict(list)
        for s, d in zip(self.mgi_samples, self.mgi_depths):
            cv_dict["sample"].append(s)
            cv_dict["platform"].append(PLATFORM2)
            df = pd.read_table(d, header=None)
            cv_dict["CV"].append(df.iloc[:, 2].std() / df.iloc[:, 2].mean())
        for s, d in zip(self.np_samples, self.np_depths):
            cv_dict["sample"].append(s)
            cv_dict["platform"].append(PLATFORM3)
            df = pd.read_table(d, header=None)
            cv_dict["CV"].append(df.iloc[:, 2].std() / df.iloc[:, 2].mean())
        pd.DataFrame(cv_dict).to_csv(self.cv_file, sep="\t", header=True, index=False)
        return self

    def cv_plot(self):
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'Rscript')} {PIPELINE.joinpath('overall_depth_cv_plot.R')} "
                f"-if {self.cv_file} "
                f"-of {self.cv_pdf} "
                f"-og {self.cv_png} "
            ],
            [
                self.sh_dir.joinpath("cv_plot.sh")
            ]
        ).multi_run(1)
        Results([self.cv_png, self.cv_pdf]).check_exists().check_empty()
        return self

    def mgi_depth_plot(self):
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'Rscript')} {PIPELINE.joinpath('overall_depth_single_plot.R')} "
                f"-if {depth_file} "
                f"-of {out_pdf} "
                f"-og {out_png} "
                for depth_file, out_pdf, out_png in zip(self.mgi_depths, self.mgi_depth_pdfs, self.mgi_depth_pngs)
            ],
            [
                self.sh_dir.joinpath(f"{s}_mgi_depth_plot.sh") for s in self.mgi_samples
            ]
        ).multi_run(self.cpu)
        Results(self.mgi_depth_pdfs + self.mgi_depth_pngs).check_exists().check_empty()
        return self

    def np_depth_plot(self):
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'Rscript')} {PIPELINE.joinpath('overall_depth_single_plot.R')} "
                f"-if {depth_file} "
                f"-of {out_pdf} "
                f"-og {out_png} "
                for depth_file, out_pdf, out_png in zip(self.np_depths, self.np_depth_pdfs, self.np_depth_pngs)
            ],
            [
                self.sh_dir.joinpath(f"{s}_np_depth_plot.sh") for s in self.np_samples
            ]
        ).multi_run(self.cpu)
        Results(self.np_depth_pdfs + self.np_depth_pngs).check_exists().check_empty()
        return self


def main(args):
    (
        OverallDepth(args)
        .get_cv_file()
        .cv_plot()
        .mgi_depth_plot()
        .np_depth_plot()
    )
