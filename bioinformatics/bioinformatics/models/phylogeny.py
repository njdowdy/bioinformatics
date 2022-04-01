from pydantic import BaseModel
from enum import Enum


class PhyloSoftware(str, Enum):
    iqtree2 = "iqtree2"
    raxmlng = "raxml-ng"
