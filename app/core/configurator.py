"""
Constants for submission details
"""

DRAFT_MODE_DETAILS = """In draft mode, submissions may be incomplete."""

NAME_DETAILS = """Please provide an informal but informative name of the method used.
The name must not exceed 40 characters."""

COMPUTE_TIME_DETAILS = """Please provide the average compute time across all of the molecules.
For physical methods, report the GPU and/or CPU compute time in hours.
For empirical methods, report the query time in hours.
Create a new line for each processor type."""

COMPUTING_AND_HARDWARE_DETAILS = """Please provide details of the computing resources that were used to train models and make predictions.
Please specify compute time for training models and querying separately for empirical prediction methods.
Provide a detailed description of the hardware used to run the simulations."""

SOFTWARE_DETAILS = """List all major software packages used and their versions.
Create a new line for each software."""

CATEGORY_OPTIONS = ("Physical (MM)", "Physical (QM)", "Empirical", "Mixed")
CATEGORY_CHOICES = list(zip(CATEGORY_OPTIONS, CATEGORY_OPTIONS))

CATEGORY_DETAILS = """State which method category your prediction method is better described as:
`Physical (MM)`, `Physical (QM)`, `Empirical`, or `Mixed`.
Pick only one category label."""

METHOD_DETAILS = """Methodology and computational details.
Level of details should be roughly equivalent to that used in a publication.
Please include the values of key parameters with units.
Please explain how statistical uncertainties were estimated.

If you have evaluated additional microstates, please report their SMILES strings and populations of all the microstates in this section.
If you used a microstate other than the challenge provided microstate (`SMXX_micro000`), please list your chosen `Molecule ID` (in the form of `SMXX_extra001`) along with the SMILES string in your methods description.

Use as many lines of text as you need."""

NOTES_DETAILS = """Any extra notes about submissions may be placed here. 
This section is intended as a "notes-to-self" section to add annotations about your submission or submission runs. Notes in this section will not be used or referenced by challenge administrators."""


RANKED_DETAILS = """All submissions must either be ranked or non-ranked.
Only one ranked submission per participant is allowed.
Multiple ranked submissions from the same participant will not be judged.
Non-ranked submissions are accepted so we can verify that they were made before the deadline."""
