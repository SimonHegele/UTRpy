import pandas

from .utrpy_gff_utils import overlapping_features, attributes_dict

class Transcript():

    def get_features(self, gff):

        features = overlapping_features(gff, self.data)

        return features[
            features.apply(
                lambda exon: False if "transcript_id" not in attributes_dict(exon) else attributes_dict(exon)["transcript_id"]==self.id,
                axis=1
            )
        ]

    def __init__(self, row: pandas.Series, gff: pandas.DataFrame):

        assert row["type"] == "transcript" or "RNA" in row["type"]

        self.id       : str              = attributes_dict(row)["ID"]
        self.data     : pandas.Series    = row
        self.features : pandas.DataFrame = self.get_features(gff)
        self.exons    : pandas.DataFrame = self.features.loc[self.features["type"]=="exon"]
    
    def delete(self, gff: pandas.DataFrame) -> pandas.DataFrame:

        for i, feature in self.features.iterrows():

            gff = gff.drop(feature.name)

        return gff.drop(self.data.name)