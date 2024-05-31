from argparse import ArgumentParser


if __name__ == "__main__":

    parser = ArgumentParser(prog="cov2_xiehe", description="协和新冠项目工具集")
    subparser = parser.add_subparsers(title="subcommands", dest="subcommand", description="协和新冠项目工具集子命令")

    get_sample_list_parser = subparser.add_parser("get_sample_list", help="获取样本名与下机数据对应关系")
    get_sample_list_parser.add_argument("-id", "--in_dir", type=str, required=True)
    get_sample_list_parser.add_argument("-of", "--out_file", type=str, required=True)

    rename_fq_parser = subparser.add_parser("rename_fq", help="根据样本与下机数据对应关系拷贝并改名")
    rename_fq_parser.add_argument("-if", "--in_file", type=str, required=True)
    rename_fq_parser.add_argument("-od", "--out_dir", type=str, required=True)
    rename_fq_parser.add_argument("-c", "--cpu", type=int, default=16)

    fastp_qc_parser = subparser.add_parser("fastp_qc", help="使用fastp对原始数据进行质控")
    fastp_qc_parser.add_argument("-id", "--in_dir", type=str, required=True, help="输入目录，包含所有的原始fq文件")
    fastp_qc_parser.add_argument("-od", "--out_dir", type=str, required=True, help="输出目录，包含所有的结果文件")
    fastp_qc_parser.add_argument("-c", "--cpu", type=int, default=8, help="使用的cpu核心数量")

    bwa_mem_parser = subparser.add_parser("bwa_mem", help="使用bwa将reads比对到参考基因组上")
    bwa_mem_parser.add_argument("-id", "--in_dir", type=str, required=True, help="输入目录（*.fq.gz）")
    bwa_mem_parser.add_argument("-od", "--out_dir", type=str, required=True, help="（*.sorted.bam, *.sorted.bam.bai）")
    bwa_mem_parser.add_argument("-bf", "--bed_file", type=str, default=None, help="mgiseq.primers.bed")
    bwa_mem_parser.add_argument("-c", "--cpu", type=int, default=8)

    coverage_parser = subparser.add_parser("coverage", help="计样本的覆盖率")
    coverage_parser.add_argument("-id", "--in_dir", type=str, required=True)
    coverage_parser.add_argument("-od", "--out_dir", type=str, required=True)
    coverage_parser.add_argument("-c", "--cpu", type=int, default=8)

    depth_parser = subparser.add_parser("depth", help="统计样本的测序深度")
    depth_parser.add_argument("-id", "--in_dir", type=str, required=True, help="bams")
    depth_parser.add_argument("-od", "--out_dir", type=str, required=True, help="depths")
    depth_parser.add_argument("-c", "--cpu", type=int, default=1)

    bcftools_call = subparser.add_parser("bcftools_call", help="获取变一下信息")
    bcftools_call.add_argument("-id", "--in_dir", type=str, required=True, help="输入目录（*.sorted.bam, *.sorted.bam.bai")
    bcftools_call.add_argument("-od", "--out_dir", type=str, required=True, help="输出目录（*.vcf）")
    bcftools_call.add_argument("-c", "--cpu", type=str, required=True, help="CPU数量")

    ivar_consensus_parser = subparser.add_parser("ivar_consensus", help="使用ivar对reads进行组装")
    ivar_consensus_parser.add_argument("-id", "--in_dir", type=str, required=True, help="输入目录，包含所有的fq文件")
    ivar_consensus_parser.add_argument("-od", "--out_dir", type=str, required=True, help="输出目录，把含所有的结果文件")
    ivar_consensus_parser.add_argument("-c", "--cpu", type=int, default=4, help="使用的cpu数量")

    select_reads_parser = subparser.add_parser("select_reads", help="从原始数据中获取新冠序列")
    select_reads_parser.add_argument("-id", "--in_dir", type=str, required=True, help="输入目录，包含所有的原始数据")
    select_reads_parser.add_argument("-od", "--out_dir", type=str, required=True, help="输出目录，包含所有的新冠序列文件")
    select_reads_parser.add_argument("-c", "--cpu", type=int,  default=8, help="使用的CPU的数量，默认为8")

    stats_reads_parser = subparser.add_parser("stats_reads", help="统计比对上新冠的reads信息")
    stats_reads_parser.add_argument("-id", "--in_dir", type=str, required=True, help="输入目录，包含所有的新冠reads")
    stats_reads_parser.add_argument("-od", "--out_dir", type=str, required=True, help="输出目录包含所有的结果文件")
    stats_reads_parser.add_argument("-c", "--cpu", type=int, default=4, help="使用的cpu核心数量")

    porechop_parser = subparser.add_parser("porechop", help="过滤三代reads")
    porechop_parser.add_argument("-id", "--in_dir", type=str, required=True)
    porechop_parser.add_argument("-od", "--out_dir", type=str, required=True)
    porechop_parser.add_argument("-c", "--cpu", type=int, default=1)

    minimap2_parser = subparser.add_parser("minimap2", help="将三代reads比对到参考基因组上")
    minimap2_parser.add_argument("-id", "--in_dir", type=str, required=True)
    minimap2_parser.add_argument("-od", "--out_dir", type=str, required=True)
    minimap2_parser.add_argument("-c", "--cpu", type=int, default=1)

    medaka_parser = subparser.add_parser("medaka", help="基因组组装和变异calling")
    medaka_parser.add_argument("-id", "--in_dir", type=str, required=True)
    medaka_parser.add_argument("-od", "--out_dir", type=str, required=True)
    medaka_parser.add_argument("-c", "--cpu", type=int, default=1)

    overall_stats_parser = subparser.add_parser("overall_stats", help="二三代数据整体统计")
    overall_stats_parser.add_argument("-mcf", "--mgi_coverage_file", type=str, required=True, help="二代覆盖率文件")
    overall_stats_parser.add_argument("-ncf", "--np_coverage_file", type=str, required=True, help="三代覆盖率文件")
    overall_stats_parser.add_argument("-msf", "--mgi_stats_file", type=str, required=True, help="二代统计文件")
    overall_stats_parser.add_argument("-nsf", "--np_stats_file", type=str, required=True, help="三代统计文件")
    overall_stats_parser.add_argument("-od", "--out_dir", type=str, required=True, help="输出目录")

    overall_coverage_analysis_parser = subparser.add_parser("overall_coverage_analysis", help="对二三代的覆盖度进行分析")
    overall_coverage_analysis_parser.add_argument("-mcf", "--mgi_coverage_file", required=True)
    overall_coverage_analysis_parser.add_argument("-ncf", "--np_coverage_file", required=True)
    overall_coverage_analysis_parser.add_argument("-od", "--out_dir", required=True)

    overall_depth_parser = subparser.add_parser("overall_depth", help="二三代数据整体测序深度的分析")
    overall_depth_parser.add_argument("-mdd", "--mgi_depth_dir", type=str, required=True)
    overall_depth_parser.add_argument("-ndd", "--np_depth_dir", type=str, required=True)
    overall_depth_parser.add_argument("-od", "--out_dir", type=str, required=True)
    overall_depth_parser.add_argument("-c", "--cpu", type=int, default=1)

    overall_consensus_parser = subparser.add_parser("overall_consensus", help="二三代基因组的比较")
    overall_consensus_parser.add_argument("-mnf", "--mgi_nextclade_file", type=str, required=True)
    overall_consensus_parser.add_argument("-nnf", "--np_nextclade_file", type=str, required=True)
    overall_consensus_parser.add_argument("-od", "--out_dir", type=str, required=True)

    overall_ct_parser = subparser.add_parser("overall_ct", help="绘制reads和coverages与ct值之间的对应关系")
    overall_ct_parser.add_argument("-mcf", "--mgi_coverage_file", type=str, required=True)
    overall_ct_parser.add_argument("-ncf", "-np_coverage_file", type=str, required=True)
    overall_ct_parser.add_argument("-od", "--out_dir", type=str, required=True)

    overall_mutation_parser = subparser.add_parser("overall_mutation", help="获取个样本的突变数量信息并绘图")
    overall_mutation_parser.add_argument("-mvd", "--mgi_vcf_dir", type=str, required=True)
    overall_mutation_parser.add_argument("-nvd", "--np_vcf_dir", type=str, required=True)
    overall_mutation_parser.add_argument("-od", "--out_dir", type=str, required=True)

    args = parser.parse_args()
    if args.subcommand == "select_reads":
        import select_reads
        select_reads.main(args)
    elif args.subcommand == "get_sample_list":
        import get_sample_list
        get_sample_list.main(args)
    elif args.subcommand == "rename_fq":
        import rename_fq
        rename_fq.main(args)
    elif args.subcommand == "fastp_qc":
        import fastp_qc
        fastp_qc.main(args)
    elif args.subcommand == "ivar_consensus":
        import ivar_consensus
        ivar_consensus.main(args)
    elif args.subcommand == "bwa_mem":
        import bwa_mem
        bwa_mem.main(args)
    elif args.subcommand == "coverage":
        import coverage
        coverage.main(args)
    elif args.subcommand == "bcftools_call":
        import bcftools_call
        bcftools_call.main(args)
    elif args.subcommand == "porechop":
        import porechop
        porechop.main(args)
    elif args.subcommand == "minimap2":
        import minimap2
        minimap2.main(args)
    elif args.subcommand == "medaka":
        import medaka
        medaka.main(args)
    elif args.subcommand == "stats_reads":
        import stats_reads
        stats_reads.main(args)
    elif args.subcommand == "overall_coverage_analysis":
        import overall_coverage_analysis
        overall_coverage_analysis.main(args)
    elif args.subcommand == "depth":
        import depth
        depth.main(args)
    elif args.subcommand == "overall_stats":
        import overall_stats
        overall_stats.main(args)
    elif args.subcommand == "overall_depth":
        import overall_depth
        overall_depth.main(args)
    elif args.subcommand == "overall_consensus":
        import overall_consensus
        overall_consensus.main(args)
    elif args.subcommand == "overall_ct":
        import overall_ct
        overall_ct.main(args)
    elif args.subcommand == "overall_mutation":
        import overall_mutation
        overall_mutation.main(args)
    else:
        parser.print_usage()
        raise Exception("Error:请输入正确的子命令！")
