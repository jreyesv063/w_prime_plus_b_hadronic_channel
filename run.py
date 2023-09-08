import json
import time
import dask     
import pickle
import argparse
import numpy as np
import datetime
from pathlib import Path
from coffea import processor
from dask.distributed import Client
from humanfriendly import format_timespan
from distributed.diagnostics.plugin import UploadDirectory
from wprime_plus_b.processors.ttbar_analysis import TtbarAnalysis
from wprime_plus_b.processors.ttbar_hadronic import TtbarAnalysis_h


def main(args):
    np.seterr(divide="ignore", invalid="ignore")
    # load and process filesets
    fileset = {}
    with open(args.fileset, "r") as handle:
        data = json.load(handle)
    for sample, val in data.items():
        if args.nfiles != -1:
            val = val[: args.nfiles]
        fileset[sample] = [f"root://{args.redirector}/" + file for file in val]

    # define processors
    processors = {
        "ttbar": TtbarAnalysis,
        "ttbar_h": TtbarAnalysis_h
    }
    processor_kwargs = {
        "year": args.year,
        "yearmod": args.yearmod,
        "channel": args.channel,
        "lepton_flavor": args.lepton_flavor,
        "syst": args.syst,
        "output_type": args.output_type,
    }
    if args.processor in ["ztoll", "btag_eff", "qcd"]:
        del processor_kwargs["channel"]
        del processor_kwargs["syst"]
        del processor_kwargs["output_type"]
    if args.processor == "btag_eff":
        del processor_kwargs["lepton_flavor"]
    if args.processor == "ttbar_h":
        del processor_kwargs["channel"]
        del processor_kwargs["lepton_flavor"]

    # define executors
    executors = {
        "iterative": processor.iterative_executor,
        "futures": processor.futures_executor,
        "dask": processor.dask_executor,
    }
    executor_args = {
        "schema": processor.NanoAODSchema,
    }
    if args.executor == "futures":
        executor_args.update({"workers": args.workers})
    if args.executor == "dask":
        client = Client(args.client)
        print(f"client: {args.client}")
        executor_args.update({"client": client})
        # upload local directory to dask workers
        print(f"trying to upload {Path.cwd()} directory")
        try:
            client.register_worker_plugin(
                UploadDirectory(f"{Path.cwd()}", restart=True, update_path=True),
                nanny=True,
            )
            print(f"Uploaded {Path.cwd()} succesfully")
        except OSError:
            print("Failed to upload the directory")
            
    # run processor
    t0 = time.monotonic()
    out = processor.run_uproot_job(
        fileset,
        treename="Events",
        processor_instance=processors[args.processor](**processor_kwargs),
        executor=executors[args.executor],
        executor_args=executor_args,
    )
    exec_time = format_timespan(time.monotonic() - t0)
    
    # get metadata
    metadata = {"walltime": exec_time}
    metadata.update({'events_before': float(out["metadata"]['events_before'])})
    metadata.update({'events_after': float(out["metadata"]['events_after'])})
    if "sumw" in out["metadata"]:
        metadata.update({'sumw': float(out["metadata"]['sumw'])})
    metadata.update({"fileset": fileset[sample]})
    
    # save args to metadata
    args_dict = vars(args).copy()
    del args_dict["fileset"]
    if args.processor in ["ztoll", "btag_eff"]:
        del args_dict["channel"]
    if args.processor == "btag_eff":
        del args_dict["lepton_flavor"]
    metadata.update(args_dict)
    
    # drop metadata from output
    del out["metadata"]
    
    # define output and metadata paths
    date = datetime.datetime.today().strftime("%Y-%m-%d")
    base_path = Path(
        args.output_location
        + "/"
        + args.tag
        + "/"
        + args.processor
        + "/"
        + date
        + "/"
    )
    ttbar_output_path = Path(
        str(base_path)
        + "/"
        + args.channel 
        + "/"
        + args.year
        + "/"
        + args.lepton_flavor
    )
    other_output_path = Path(
        str(base_path)
        + "/"
        + args.year
        + "/"
        + args.lepton_flavor
    )
    ttbar_h_output_path = Path(
        str(base_path)
        + "/"
        + args.year
    )
    output_path = {
        "ttbar": ttbar_output_path,
        "ztoll": other_output_path,
        "qcd": other_output_path,
        "ttbar_h": ttbar_h_output_path
    }
    # save output
    if not output_path[args.processor].exists():
        output_path[args.processor].mkdir(parents=True)
    with open(f"{str(output_path[args.processor])}/{sample}.pkl", "wb") as handle:
        pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    # save metadata
    metadata_path = Path(f"{str(output_path[args.processor])}/metadata")
    if not metadata_path.exists():
        metadata_path.mkdir(parents=True)
    with open(f"{metadata_path}/{sample}_metadata.json", "w") as f:
        f.write(json.dumps(metadata))

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--processor",
        dest="processor",
        type=str,
        default="ttbar_cr1",
        help="processor to be used {trigger, ttbar_cr1, ttbar_cr2, candle, btag_eff}",
    )
    parser.add_argument(
        "--executor",
        dest="executor",
        type=str,
        default="iterative",
        help="executor to be used {iterative, futures, dask}",
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
    parser.add_argument("--year", dest="year", type=str, default="2017", help="year")
    parser.add_argument(
        "--yearmod",
        dest="yearmod",
        type=str,
        default="",
        help="year modifier {'', 'APV'}",
    )
    parser.add_argument(
        "--nfiles",
        dest="nfiles",
        type=int,
        default=1,
        help="number of .root files to be processed by sample (default 1. To run all files use -1)",
    )
    parser.add_argument(
        "--workers",
        dest="workers",
        type=int,
        default=4,
        help="number of workers to use with futures executor (default 4)",
    )
    parser.add_argument(
        "--redirector",
        dest="redirector",
        type=str,
        default="xcache",
        help="redirector to find CMS datasets {use 'xcache' at coffea-casa. Use 'cmsxrootd.fnal.gov', 'xrootd-cms.infn.it' or 'cms-xrd-global.cern.ch' at lxplus}",
    )
    parser.add_argument(
        "--output_location",
        dest="output_location",
        type=str,
        default="./outfiles/",
        help="output location (default ./outfiles)",
    )
    parser.add_argument(
        "--tag",
        dest="tag",
        type=str,
        default="test",
        help="tag of the submitted jobs",
    )
    parser.add_argument(
        "--fileset",
        dest="fileset",
        type=str,
        help="json fileset",
    )
    parser.add_argument(
        "--client",
        dest="client",
        type=str,
        help="dask client to use with dask executor on coffea-casa",
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
    args = parser.parse_args()
    main(args)