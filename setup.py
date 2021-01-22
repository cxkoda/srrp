import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="srrp-cxkoda", # Replace with your own username
    version="0.0.1",
    author="David Huber",
    author_email="dave@yomi.eu",
    description="An exact solver for special relativistic riemann problems",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cxkoda/srrp",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)