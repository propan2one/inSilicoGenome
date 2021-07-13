from Bio.Seq import Seq
import numpy as np

def random_dnasequence(size):
    """
    Generate random sequence of DNA

    Parameters
    ----------
    size : Size of DNA in bp
      An integer between 1 and 100 000 000.

    Returns
    -------
    pandas.core.arrays.categorical.Categorical
      The new concatenated pandas categorical.

    Examples
    --------
    >>> from insilicogenome import insilicogenome
    >>> insilicogenome.random_dnasequence(10)
    Seq('GTTCTTGAT')
    """
    if ((size>= 0) and (size<= 100000000)):
        pass  # empty file
    elif (size > 100000000):  # on a title line
        raise ValueError(
            "The size of the sequence to generate is greater than the maximum size of 100000000bp"
            f"The value of '{size}' is too big, it must be reduced"
        )
    else:
        assert size > 0, "Size value negatif ; this should be impossible!"
    np.random.seed(123)
    sequence = ''.join(np.random.choice(('C','G','T','A'), size ))
    return Seq(sequence)