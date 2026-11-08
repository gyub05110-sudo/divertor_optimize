"""Record the CPU, GPU and dependency environment used for timed prospective runs."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import numpy
import pyarrow
import scipy
import shapely


def capture(study):
    runtime = Path(os.environ.get("DIVERTOR_GPU_RUNTIME", "D:/DivertorOptimization/runtime/cupy"))
    gpu = subprocess.check_output(["nvidia-smi", "--query-gpu=name,driver_version,memory.total,compute_cap",
                                   "--format=csv,noheader"], text=True).strip()
    if len(gpu.splitlines()) != 1:
        raise ValueError("this serial benchmark expects one visible GPU")
    modules = (numpy, scipy, shapely, pyarrow)
    with Path(sys.executable).open("rb") as stream:
        binary = hashlib.file_digest(stream, "sha256").hexdigest()
    return dict(utc=datetime.now(timezone.utc).isoformat(), python=sys.version, executable=sys.executable,
                executable_sha256=binary, platform=platform.platform(), cpu=platform.processor(),
                logical_processors=os.cpu_count(), gpu_name_driver_memory_compute=gpu, gpu_runtime=str(runtime),
                packages={m.__name__: dict(version=m.__version__, path=m.__file__) for m in modules},
                gpu_distributions={d.metadata["Name"]:d.version for d in importlib.metadata.distributions(path=[str(runtime)])},
                start="first timestamp inside the public Python CLI; outer process observation is separate",
                omp_num_threads=1, openblas_num_threads=1, mkl_num_threads=1,
                study=str(study.resolve()), runtime_cache=str((study/"runtime").resolve()),
                cache_policy="fresh process per measured search; no shared solver state; OS and compilation caches shared")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(capture(args.study), stream, indent=2)
