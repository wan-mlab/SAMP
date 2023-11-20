from setuptools import setup, find_packages

setup(
    name="SAMP",
    version="1.0.0",
    description="SAMP: a Split amino acid composition based AntiMicrobial Prediction model ",
    url="https://github.com/wan-mlab/SAMP",
    author="Junxi Feng, Mengtao Sun, Cong Liu, Shibiao Wan",
    author_email="junxifeng@hsph.harvard.edu",
    license="MIT",
    packages=find_packages(where='./SAMP'),
    package_dir={
        '': 'SAMP'
    },
    include_package_data=True,
    install_requires=[
        "scikit-learn==1.3.0",
        "pandas==1.1.3",
        "numpy==1.23.5"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.6"
)