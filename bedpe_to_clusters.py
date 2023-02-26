from __future__ import annotations

import json

import numpy as np
import pandas as pd


def bedpe_to_clusters(
    bedpe_filename: str, matrix_resolution: int
) -> dict[str, list[set[tuple[int, int]]]]:
    """
    Loads loop information from a .bedpe file from disk into the hic3defdr
    sparse cluster format.

    Trans loops are ignored.

    Parameters
    ----------
    bedpe_filename
        Path to a .bedpe file on disk to load loop information from.
    matrix_resolution
        The matrix resolution of the contact matrix in base pairs.

    Returns
    -------
    dict[list[set[tuple[int, int]]]]
        The outer dict's keys are chromosome names as strings, and the values
        for each key are the list of clusters found on that chromosome. Each
        cluster is represented as a set. The elements in the set are pixels,
        represented as `(int, int)` tuples, where the integers refer to the
        0-based bin indices of the contact matrix that identify the pixel.
    """
    # parse .bedpe file
    bedpe_df = pd.read_csv(
        bedpe_filename,
        sep="\t",
        usecols=range(6),
        names=["chrom1", "start1", "end1", "chrom2", "start2", "end2"],
    )

    # identify set of all chromosomes
    chroms = set(bedpe_df.chrom1.unique()).union(set(bedpe_df.chrom1.unique()))

    return {
        # one key, value pair per chromosome
        # key is the chromosome name
        chrom: [
            # value is a list of clusters on that chromosome
            # each cluster comes from a row of the .bedpe table
            bedpe_row_to_cluster(row, matrix_resolution)
            # subset the .bedpe table to only cis loops on this chromosome
            for _, row in bedpe_df.query(
                "chrom1 == @chrom and chrom2 == @chrom"
            ).iterrows()
        ]
        for chrom in chroms
    }


def bedpe_row_to_cluster(
    bedpe_row: pd.Series, matrix_resolution: int
) -> set[tuple[int, int]]:
    """
    Convert a single row from a .bedpe file (represented as a pandas Series) to
    a hic3defdr-format sparse cluster.

    Parameters
    ----------
    bedpe_row
        Row from a .bedpe file, represented as a pandas Series. Must have at
        least the following columns: start1, end1, start2, end2.
    matrix_resolution
        The matrix resolution of the contact matrix in base pairs.

    Returns
    -------
    set[tuple[int, int]]
        Each tuple in the set represents a pixel in the cluster as an
        `(int, int)` pair where the integers refer to the 0-based bin indices of
        the contact matrix that identify the pixel.
    """
    row_min = np.floor(bedpe_row["start1"] / matrix_resolution).astype(int)
    col_min = np.floor(bedpe_row["start2"] / matrix_resolution).astype(int)
    row_max = np.ceil(bedpe_row["end1"] / matrix_resolution).astype(int) - 1
    col_max = np.ceil(bedpe_row["end2"] / matrix_resolution).astype(int) - 1
    return {
        (i, j) for i in range(row_min, row_max + 1) for j in range(col_min, col_max + 1)
    }


class NumpyEncoder(json.JSONEncoder):
    """
    Pass this to `json.dump()` to correctly serialize numpy values.

    Credit: https://stackoverflow.com/a/27050186

    Copied from `hic3defdr.util.clusters.NumpyEncoder`.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NumpyEncoder, self).default(obj)


def save_clusters(clusters: list[set[tuple[int, int]]], outfile: str) -> None:
    """
    Saves cluster information to disk in sparse JSON format.

    Paraphrased from `hic3defdr.util.clusters.save_clusters()`.

    Parameters
    ----------
    clusters
        The sets are clusters, the tuples are the indices of pixels in that
        cluster.
    outfile
        File to write JSON output to.
    """
    with open(outfile, "w") as handle:
        json.dump([[[i, j] for i, j in c] for c in clusters], handle, cls=NumpyEncoder)


if __name__ == "__main__":
    # get clusters on all chromosomes from input .bedpe file
    all_clusters = bedpe_to_clusters("loops.bedpe", 10000)

    # save each chromosome's clusters to a separate .json file
    for chrom, clusters in all_clusters.items():
        save_clusters(clusters, f"{chrom}_clusters.json")
