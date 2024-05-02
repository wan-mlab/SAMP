import subprocess

class FeatureBuilder:
    def __init__(self, left_prop, middle_prop):
        """
        Initializes the FeatureBuilder with split proportions.

        Parameters:
        ----------------
        left_prop (float): First split (left-side) fraction of the peptide sequence.
        middle_prop (float): Second split (middle) fraction of the peptide sequence.

        * left_prop and middle_prop must be less than 1.
        """
        self.left_prop = left_prop
        self.middle_prop = middle_prop

    def build(self, ampfile, nonampfile, out):
        """
        Builds features by calling an R script.

        Parameters:
        ----------------
        ampfile (str): Path to a fasta format AMP file, or a text file in fasta format.
        nonampfile (str): Path to a fasta format Non-AMP file, or a text file in fasta format.
        out (str): Name of the output file.

        Returns:
        ----------------
        Writes the feature matrix to the 'features' folder.
        """
        command = f"Rscript buildFeatures.R {ampfile} {nonampfile} {out} {self.left_prop} {self.middle_prop}"
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

        # Print output
        for line in process.stdout:
            print(line.decode("utf-8"), end='')

        # Wait for the command to finish
        process.wait()