"""
Logging utility functions.
"""

import logging


def setup_logger(logger: logging.Logger, filename: str):
    """Setup a logger object by adding formatted File and Stream handlers.

    Args:
        logger (logging.Logger): logger object.
        filename (str): log file name.
    """
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(filename)
    ch = logging.StreamHandler()
    fh_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh.setFormatter(fh_formatter)
    ch_formatter = logging.Formatter(
        "%(asctime)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    ch.setFormatter(ch_formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
