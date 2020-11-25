import pandas as pd
import os
from script_utils import show_output


def get_off_coverage(df, min_coverage=100, max_dist=20):
    """
    find the regions with target == 0
    reduce them to blocks and output the center coords IGVnav
    """

    new_cols = ["Start", "End", "meanCov", "maxCov"]
    # get the off_coverage df
    off_df = df.query("onTarget == 0 and Coverage >= @min_coverage ")
    # get the distance between adjacent positions (0 for chrom jumps)
    off_df.loc[:, "dist"] = (off_df["Pos"] - off_df.shift(1)["Pos"]) * (
        off_df["Chr"] == off_df.shift(1)["Chr"]
    ).astype(int)
    # find the gaps where adjacent distance is greater max_dist (or chrom jump)
    off_df.loc[:, "gap"] = ((off_df["dist"] == 0) | (off_df["dist"] > max_dist)).astype(
        int
    )
    # make blocks from the gaps
    off_df.loc[:, "block"] = off_df["gap"].cumsum()
    # get aggregate info for each block
    block_df = off_df.groupby("block").agg(
        Chr=("Chr", "first"),
        meanCov=("Coverage", "mean"),
        maxCov=("Coverage", "max"),
        Start=("Pos", "min"),
        End=("Pos", "max"),
    )
    for col in new_cols:
        block_df.loc[:, col] = block_df[col].astype(int)
    cols = ["Chr"] + new_cols
    return block_df[cols]


def make_IGV_nav(df):
    """
    converts a Chr, Start, End, dataFrame into IGVnav txt file
    """

    # create the columns
    df = df.copy().loc[:, ["Chr", "Start", "End"]]
    for col in ["Call", "Tags", "Notes"]:
        df[col] = ""
    # remove the "chr"
    df.loc[:, "Chr"] = df["Chr"].str.replace("chr", "")
    return df


def main(s):
    """
    wrapped into function lest module be called by each subprocess
    """

    c = s.config
    cc = c["offTarget"]
    w = s.wildcards
    p = s.params
    i = s.input
    o = s.output

    show_output(f"Reading coverage data for {w.sample} from {i.cov}", end="")
    off_target_df = pd.read_csv(i.cov, sep="\t")
    show_output(". Done.")

    # converting to block df
    block_df = get_off_coverage(
        off_target_df, min_coverage=cc["minCov"], max_dist=cc["maxDist"]
    )
    if block_df.empty:
        show_output(
            f"No offTarget coverage found for {w.sample}. Creating empty file",
            color="warning",
        )
        block_df.to_csv(o.off, sep="\t")
        pd.DataFrame(columns=["Chr", "Start", "End", ["Call", "Tags", "Notes"]]).to_csv(
            o.IGVnav, sep="\t", index=False
        )

    # write block_df to file
    block_df.to_csv(o.off, sep="\t")

    # make IGVnav file
    igv_df = make_IGV_nav(block_df)

    # write IGVnav_df to file
    igv_df.to_csv(o.IGVnav, sep="\t", index=False)
    show_output(f"IGVnav file for {w.sample} written to {o.IGVnav}")


if __name__ == "__main__":
    main(snakemake)
