import logging
import pandas

from .utrpy_gff_utils import overlapping_features, attributes_dict

class Transcript():

    def is_transcript_feature(self, feature: pandas.Series) -> bool:

        attrs = attributes_dict(feature)

        if "transcript_id" in attrs.keys() and attrs["transcript_id"]==self.id:
            return True
        if "Parent" in attrs.keys() and attrs["Parent"]==self.id:
            return True
        return False

    def get_features(self, gff):

        features = overlapping_features(gff, self.data)

        return features[
            features.apply(
                lambda feature: self.is_transcript_feature(feature),
                axis=1
            )
        ]

    def __init__(self, row: pandas.Series, gff: pandas.DataFrame):

        assert row["type"] == "transcript" or "RNA" in row["type"]

        self.id       : str              = attributes_dict(row)["ID"]
        self.data     : pandas.Series    = row
        self.features : pandas.DataFrame = self.get_features(gff)
        self.exons    : pandas.DataFrame = self.features.loc[self.features["type"]=="exon"]

        if len(self.features)==0:
            logging.error(f"No features found for transcript {self.id}")
        if len(self.exons)==0:
            logging.error(f"No exons found for transcript {self.id}")
    
    def delete(self, gff: pandas.DataFrame) -> pandas.DataFrame:

        for i, feature in self.features.iterrows():

            gff = gff.drop(feature.name)

        return gff.drop(self.data.name)