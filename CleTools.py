import pandas

class AssayQC:
    def __init__(self, file=None, df=None):
        """
        Initialize the CovQC class with a pandas DataFrame, and two lists: genes and targets.

        :param dataframe: pandas DataFrame to be used in the class.
        :param genes: List of genes.
        :param targets: List of targets.
        """
        self.covdf = dataframe if isinstance(dataframe, pd.DataFrame) else pd.DataFrame()
        self.genes = genes if isinstance(genes, list) else []
        self.targets = targets if isinstance(targets, list) else []


class CovQC:
    def __init__(self, file=None, df=None):
        """
        Initialize the CovQC class with a pandas DataFrame, and two lists: genes and targets.

        :param dataframe: pandas DataFrame to be used in the class.
        :param genes: List of genes.
        :param targets: List of targets.
        """
        self.covdf = dataframe if isinstance(dataframe, pd.DataFrame) else pd.DataFrame(dataframe)
        self.genes = genes if isinstance(genes, list) else []
        self.targets = targets if isinstance(targets, list) else []

class SmallVariants:
    def __init__(self, file=None, df=None):
        """
        Initialize the CovQC class with a pandas DataFrame, and two lists: genes and targets.

        :param dataframe: pandas DataFrame to be used in the class.
        :param genes: List of genes.
        :param targets: List of targets.
        """
        self.covdf = dataframe if isinstance(dataframe, pd.DataFrame) else pd.DataFrame(dataframe)
        self.genes = genes if isinstance(genes, list) else []
        self.targets = targets if isinstance(targets, list) else []


    class SV:
        def __init__(self, file=None, df=None):
            """
            Initialize the CovQC class with a pandas DataFrame, and two lists: genes and targets.

            :param dataframe: pandas DataFrame to be used in the class.
            :param genes: List of genes.
            :param targets: List of targets.
            """
            self.covdf = dataframe if isinstance(dataframe, pd.DataFrame) else pd.DataFrame(dataframe)
            self.genes = genes if isinstance(genes, list) else []
            self.targets = targets if isinstance(targets, list) else []