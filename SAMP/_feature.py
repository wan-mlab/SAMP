import subprocess

class FeatureBuilder:
    def __init__(self, split1_prop, split2_prop):
        """
        Initializes the FeatureBuilder with split proportions.

        Parameters:
        ----------------
        split1_prop (float): First split (left-side) fraction of the peptide sequence.
        split2_prop (float): Second split (middle) fraction of the peptide sequence.

        * split1_prop and split2_prop must be less than 1.
        """
        self.split1_prop = split1_prop
        self.split2_prop = split2_prop

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
        command = f"Rscript buildFeatures.R {ampfile} {nonampfile} {out} {self.split1_prop} {self.split2_prop}"
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

        # Print output
        for line in process.stdout:
            print(line.decode("utf-8"), end='')

        # Wait for the command to finish
        process.wait()