import os,sys,time
import logging

logger = logging.getLogger(__name__)

def log_setting(prefix):
    """
    Open run log file 

    Args:
        prefix( the prefix of the log file )(str): input
    """
    logger.setLevel(level = logging.INFO)
    handler = logging.FileHandler(prefix+"_log.txt")
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(funcName)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
