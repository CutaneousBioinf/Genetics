import numpy as np
import pandas as pd
import click


@click.command()
@click.option("--path1")
@click.option("--path2")
@click.option("--gap", default=500)
def compareLoci(path1, path2, gap):
    df1 = pd.read_csv(path1)
    df2 = pd.read_csv(path2)

    loci1 = df1.to_numpy()
    loci2 = df2.to_numpy()

    # delete extra first column that Pandas creates
    np.delete(loci1, 0, 1)
    np.delete(loci2, 0, 1)

    similarloci = []
    for i in range(0, len(loci1)):
        for j in range(0, len(loci2)):
            if abs(float(loci1[i][0]) - float(loci2[j][0])) <= gap and int(loci1[i][2]) == int(loci2[i][2]):
                similarloci.append((np.asarray(loci1[i]), np.asarray(loci2[j])))

    return similarloci


if __name__ == "__main__":
    print(compareLoci("loci.csv", "loci.csv"))
