:mod:`insilicogenome.cli`
=========================

.. py:module:: insilicogenome.cli


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   insilicogenome.cli.main



Attributes
~~~~~~~~~~

.. autoapisummary::

   insilicogenome.cli.help_message
   insilicogenome.cli.output


.. data:: help_message
   :annotation: = Multiline-String

    .. raw:: html

        <details><summary>Show Value</summary>

    .. code-block:: text
        :linenos:

        
        usage: insilicogenome -o|--output <string> [-s|--size] <integer>

        A Python utility to create artificial genomes, this tool is designed to be used in bioinformatics benchmarking programs.

        Arguments :
            -h, --help    show this help message and exit.
            -o, --output  Name of the genome and header. [string]
            -s, --size    Size in base paire of the genome < 100 000 000 bp. [int]

         Examples
            --------
            insilicogenome -o genome.fasta -s 100


    .. raw:: html

        </details>

   

.. data:: output
   

   

.. function:: main()


