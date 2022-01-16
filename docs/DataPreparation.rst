Data Preparation
================

Identify the leaflet(s) of interest, particularly for lipid bilayers. The outer and inner leaflets are treated separately in the analysis so care must be taken to differentiate the two. A lipid, often cholesterol, may flip between inner and outer leaflets. Other molecules of interest (e.g. a surfactant) may spend long periods of time on neither leaflets: either within the hydrophobic tail region or in the bulk solvent.

Produce input-friendly coordinates. MD simulations often make use of periodic boundary conditions. For visualization, the vesicle's "true" center must be tracked over time and the coordinates unwrapped so as to keep the center fixed. This kind of unwrapping can be done as follows:

``
blah blah blah
``

Condense coordinates into coordinates of interest. Often there is too much information to process and so only select coordinates are follows. For an atomistic simulation of phospholipids, this may mean only printing the coordinates of the phosphorus heads.
