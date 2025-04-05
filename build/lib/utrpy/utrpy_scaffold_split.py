from logging import info
from pandas  import DataFrame

def scaffold_split(gff: DataFrame) -> dict[str, DataFrame]:
    """
    Splits a GFF DataFrame by scaffold name and returns a dictionary where:
    - Keys are scaffold names.
    - Values are DataFrame subsets containing only the corresponding scaffold's entries.

    Parameters:
        gff (DataFrame): A pandas DataFrame representing a GFF file.

    Returns:
        dict[str, DataFrame]: A dictionary mapping scaffold names to their respective GFF subsets.
    """
    info("Scaffold splitting ...")
    return {scaffold: gff[gff[0] == scaffold].sort_values(by=gff.columns[3])
            .reset_index(drop=True) 
            for scaffold in gff[0].unique()}
