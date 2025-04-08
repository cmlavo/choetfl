"""
Function definitions for general operations.
"""

from box import Box
import os
import yaml

def load_config(path: str ="default.yml") -> Box:
    """
    Load a config .yml file which specifies configurations for all scripts.
    Args:
        path (str): Path to the config .yml file. Default is "configs/default.yml".
    Return:
        config (box.Box): Box dictionary containing the configurations.
    """
    
    with open(os.path.join("configs/", path), "r") as f:
        config = yaml.safe_load(f)
    return Box(config)
