import Init
import pickle
import os
import shutil
import sys

VERSION = "1.0"


class Steps:
    """
    Class that holds the parameters and the different steps of an analysis.
    """

    def __init__(
        self,
        step: str,
        data_file: str = None,
        dAlTree_file: str = None,
        query_file: str = None,
        log_file: str = None,
        config: str = None,
    ) -> None:
        """
        Steps class constructor.

        Args:
            step (str): step to be run.
            data_file (str): name of the pickle file containing data object.
            dAlTree_file (str): name of the pickle file containing dAlTree object.
            query_file (str): name of the query file.
            log_file (str): name of the log file.
            config (dict): snakemake config object.
        """
        # Get step and debug command line arguments values
        self.step = step
        self.debug = config["debug"]

        # Read parameters from config.yaml
        self.parameters = Init.paramDef(
            self.step, config, inf="", queryName="", outdir=""
        )

        # Init data object if None
        if data_file is None:
            self.Data = Init.initData(self.parameters)
            self.Data.setGenAttr(self.step)
        else:
            with open(data_file, "rb") as f:
                self.Data = pickle.load(f)

        # Init dAlTree object if None
        if dAlTree_file is None:
            self.dAlTree = {}
        else:
            with open(dAlTree_file, "rb") as f:
                self.dAlTree = pickle.load(f)

        # Update queryFile
        if query_file is not None:
            self.Data.queryFile = query_file

        # Init logger
        Init.initLogger(self.Data, log_file, self.debug, VERSION)

        # Fakerun
        # TODO : Remove when not needed anymore
        if config["fakerun"]:
            all_copied = True
            for f in config["outputs"]:
                source_file = f.replace("results/", "results_ok/")
                if os.path.exists(source_file):
                    print(f"FAKERUN --- Copying {source_file} to {f}")
                    shutil.copy(source_file, f)
                else:
                    all_copied = False
            if all_copied:
                sys.exit(0)

        # Force basename for output files
        self.Data.baseName = "out"

    def get_params(self, params: list[str]) -> dict:
        """
        Return a subset of parameters as a dict.
        """
        return {p: self.parameters[p] for p in params}

    def serialize_data(self, filename: str) -> None:
        """
        Serialize self.Data to a file with pickle

        Args:
            filename (str): file name to dump data to.
        """
        with open(filename, "wb") as f:
            pickle.dump(self.Data, f)

    def serialize_dAlTree(self, filename: str) -> None:
        """
        Serialize self.dAlTree to a file with pickle

        Args:
            filename (str): file name to dump data to.
        """
        with open(filename, "wb") as f:
            pickle.dump(self.dAlTree, f)
