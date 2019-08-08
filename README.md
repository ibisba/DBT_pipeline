# DBT_pipeline
All tools associated with the Design Build Test pipeline for pathway discovery through RetroPath

## galaxy_RetroPath2.0
Contains the Dockerfile to build RetroPath2.0 docker image used by a Galaxy instance. Note: the rules are hardcoded in the
pipeline by adding the rules_rall.csv in the cache folder for it to be copied to the docker. Point your galaxy
tools list to the retroPath_galaxy.xml and retroPath_forward_galaxy.xml to add the tools after building the docker image.

## galaxy_rp2paths
Tools and Dockerfile to run rp2paths in galaxy

## rpRanker
This is the core code for the analysis of RetroPpath pathways. It contains the functions to run convert the output of RetroPath
to SBML, add cofactors to monocomponent reactions, analyse the pathways using FBA, thermodynamics and Selenzyme. The folder 
contains a Dockerfile that is the basis of all the other tools (excluding the above two described). Indeed all other tools
call the ibisba/rpbase docker image that contains rpRanker inside it.
