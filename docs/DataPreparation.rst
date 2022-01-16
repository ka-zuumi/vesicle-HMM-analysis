Data Preparation
================

.. warning::

  Files holding microseconds of atomistic MD simulations can be quite large. Estimate space requirements before proceeding and skip creating files for intermediate steps if they take too much space.

Identify the leaflet(s) of interest, particularly for lipid bilayers. The outer and inner leaflets are treated separately in the analysis so care must be taken to differentiate the two. A lipid, often cholesterol, may flip between inner and outer leaflets. Other molecules of interest (e.g. a surfactant) may spend long periods of time on neither leaflets: either within the hydrophobic tail region or in the bulk solvent.

Produce input-friendly coordinates. MD simulations often make use of periodic boundary conditions. For visualization, the vesicle's "true" center must be tracked over time and the coordinates unwrapped so as to keep the center fixed. This kind of unwrapping can be done for an example pdb file ``example.pdb`` as follows:

.. code-block::

  python blahblah.py example.pdb > unwrapped-example.pdb

Condense data into coordinates of interest. Often there is too much information to process and so only select coordinates are follows. For an atomistic simulation of phospholipids, this may mean only printing the coordinates of the phosphorus heads. In the example below, the script ``blahblah.py`` prints only the coordinates of the centers of mass of each phospholipid tail or cholesterol.

.. code-block::

  python blahblah.py unwrapped-example.pdb > COM-example.pdb
