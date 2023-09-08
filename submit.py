import argparse
from submit.submit_coffea_casa import run_coffea_casa
from submit.submit_lxplus import run_lxplus


def main(args):
    """
    Submit jobs at Coffea-Casa or at lxplus using HTCondor
    """
    if args.facility == "coffea-casa":
        run_coffea_casa(args)
    elif args.facility == "lxplus":
        run_lxplus(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--facility",
        dest="facility",
        type=str,
        default="coffea-casa",
        help="facility to run jobs {'coffea-casa', 'lxplus'} (default coffea-casa)",
    )
    parser.add_argument(
        "--redirector",
        dest="redirector",
        type=str,
        default="xcache",
        help="redirector to find CMS datasets {use 'xcache' at coffea-casa. use 'cmsxrootd.fnal.gov', 'xrootd-cms.infn.it' or 'cms-xrd-global.cern.ch' at lxplus} (default xcache)",
    )
    parser.add_argument(
        "--processor",
        dest="processor",
        type=str,
        default="ttbar",
        help="processor to be used {trigger, ttbar, candle, btag_eff} (default ttbar)",
    )
    parser.add_argument(
        "--executor",
        dest="executor",
        type=str,
        default="iterative",
        help="executor to be used {iterative, futures, dask} (default iterative)",
    )
    parser.add_argument(
        "--workers",
        dest="workers",
        type=int,
        default=4,
        help="number of workers to use with futures executor (default 4)",
    )
    parser.add_argument(
        "--year",
        dest="year",
        type=str,
        default="2017",
        help="year of the data {2016, 2017, 2018} (default 2017)",
    )
    parser.add_argument(
        "--yearmod",
        dest="yearmod",
        type=str,
        default="",
        help="year modifier {'', 'APV'} (default '')",
    )
    parser.add_argument(
        "--channel",
        dest="channel",
        type=str,
        default="2b1l",
        help="channel to be processed {'2b1l', '1b1e1mu'}",
    )
    parser.add_argument(
        "--lepton_flavor",
        dest="lepton_flavor",
        type=str,
        default="mu",
        help="lepton flavor to be processed {'mu', 'ele'}",
    )
    parser.add_argument(
        "--fileset",
        dest="fileset",
        type=str,
        default="UL",
        help="name of a json file at `wprime_plus_b/fileset` (default `wprime_plus_b/fileset/fileset_{year}_UL_NANO.json`)",
    )
    parser.add_argument(
        "--sample",
        dest="sample",
        type=str,
        default="all",
        help="sample key to be processed {'all', 'mc' or <sample_name>} (default all)",
    )
    parser.add_argument(
        "--nfiles",
        dest="nfiles",
        type=int,
        default=1,
        help="number of .root files to be processed by sample. To run all files use -1 (default 1)",
    )
    parser.add_argument(
        "--nsplit",
        dest="nsplit",
        type=int,
        default=1,
        help="number of subsets to divide the fileset into (default 1)",
    )
    parser.add_argument(
        "--tag",
        dest="tag",
        type=str,
        default="test",
        help="tag of the submitted jobs (default test)",
    )
    parser.add_argument(
        "--eos",
        dest="eos",
        type=bool,
        default=False,
        help="wheter to copy or not output files to EOS (default False)",
    )
    parser.add_argument(
        "--nsample",
        dest="nsample",
        type=list,
        default=[],
        help="nsample",
    )
    parser.add_argument(
        "--chunksize",
        dest="chunksize",
        type=int,
        default=50000,
        help="number of chunks to process",
    )
    parser.add_argument(
        "--output_type",
        dest="output_type",
        type=str,
        default="hist",
        help="type of output {hist, array}",
    )
    parser.add_argument(
        "--syst",
        dest="syst",
        type=str,
        default="nominal",
        help="systematic to apply {'nominal', 'jet', 'met', 'full'}",
    )
    main(parser.parse_args())