"""Single-cell RNA-seq analysis pipeline for paired-timepoint CellRanger outputs."""

from .config import Config
from .pipeline import run

__all__ = ["Config", "run"]
