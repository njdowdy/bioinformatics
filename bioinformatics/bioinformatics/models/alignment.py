from pydantic import BaseModel
from enum import Enum


class AlignmentSoftware(str, Enum):
    muscle5 = "muscle5"
    muscle3 = "muscle3"
    clustal_omega = "clustal_omega"


class AlignmentOutputFormat(str, Enum):
    fa = "fasta"
    clu = "clu"
    msf = "msf"
    phy = "phylip"
    selex = "selex"
    st = "stockholm"
    vie = "vienna"


class TrimSoftware(str, Enum):
    trimAl = "trimal"
    bmge = "bmge"
