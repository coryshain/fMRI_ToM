# No Evidence of Theory of Mind Reasoning in the Human Language Network

This repository contains code to reproduce analyses from Shain, Paunov, Chen et al. (2022).
Scripts assume the following dependencies:

R: lme4
Python: pandas, matplotlib, seaborn, statsmodels, scipy

All scripts are contained in the `fmri_tom` directory.
Python scripts can be run as follows:

    python -m fmri_tom.<SCRIPT_NAME> <ARGS>

For example, to generate plots of the distribution of linguistic features among the items of the non-verbal ToM task, run (no arguments):

    python -m fmri_tom.analyze_items

All source data are available on [OSF](https://osf.io/bzwm8/).
Scripts assume that source data are located at `../../data/fMRI_ToM/` relative to this directory, and that the directory is structured as it is on OSF.
To point them to a different directory, replace the path above with the correct path in each script.

The following subsections describe the function of each script.

## analyze_items.py

Create plots of linguistic feature distributions within the verbal ToM task:

    ppython -m fmri_tom.analyze_items

## dice.py

Compute DICE scores of overlap between fROI vs whole-brain uncorrected activation maps:

    python -m fmri_tom.dice

## extract_texts.py

Generate word-by-word tables of stimulus items:

    python -m fmri_tom.extract_texts

## extract_tom2.py

Generate word-by-word tables of stimulus items from the second verbal ToM experiment of Deen et al 2015:

    python -m fmri_tom.extract_texts

## irc.py

Compute inter-region correlations and generate plots:

    python -m fmri_tom.irc

## regress_l2.R

Run group-level mixed effects analyses:

    ./regress_l2.R

## regress_l2_001.R

Run group-level mixed effects analyses using a top 1% voxel selection threshold:

    ./regress_l2_001.R

## regress_ling.py

Residualize linguistic regressors from the language and ToM networks' responses to the verbal ToM tasks:

    python -m fmri_tom.regress_ling

## signif_table.py

Save table of effects and FDR-corrected _p_ values:

    python -m fmri_tom.signif_table

## signif_table_001.py

Save table of effects and FDR-corrected _p_ values using a top 1% voxel selection threshold:

    python -m fmri_tom.signif_table_001

## spcorr.py

Compute and plot cross-run spatial correlation values:

    python -m fmri_tom.spcorr

## References

Shain, C., Paunov, A., Chen, X., Lipkin, B., Fedorenko, E. (accepted). No Evidence of Theory of Mind Reasoning in the Human Language Network. _Cerebral Cortex_.
