"""Module containing helper functions for data manipulation"""

import pandas as pd


def owu_from_csv(csv_path: str) -> pd.MultiIndex:
    """
    Generates an observation-wise unfolded matrix from the csv file obtained from 'generate_data'

    Parameters
    ----------
    csv_path : str
        path to the .csv file containing the generated data

    Returns
    --------
    owu : pd.MultiIndex
        Observation-wise unfolded matrix that resulted from the.csv data



    """
    pd.options.mode.chained_assignment = None  # default='warn'
    # Import the generated data into a pandas DataFrame
    df = pd.read_csv(csv_path)

    # Identify the indexes where the runs start and create a column for runs
    run_start_ix = df[df["timestamps"] == 0].index
    n_runs = len(run_start_ix)
    df["run"] = -1

    # Place the correct indexes in the "runs" columns
    for run in range(n_runs - 2):
        df["run"].loc[run_start_ix[run] : run_start_ix[run + 1]] = run

    # Place the run number on the last indexes
    df["run"].loc[run_start_ix[n_runs - 1] :] = n_runs - 1

    # Turn the timestamps and the runs into multiindexes
    indexes = pd.MultiIndex.from_frame(df[["run", "timestamps"]])
    df = df.drop(columns=["run", "timestamps"])
    owu = pd.DataFrame(df.values, columns=df.columns, index=indexes)

    return owu

    from matplotlib import pyplot as plt

    plt.x_axis


