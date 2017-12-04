#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

inputs:
  colore_config: File
  fastcat_config: File
  namaster_config: File
  likelihood_config: File

outputs:
  Omega:
    type: float
    outputSource: LikelihoodAnalysis/Omega

  bias:
    type: float
    outputSource: LikelihoodAnalysis/bias


steps:
  CoLoRe:
    run: cwl/colore.cwl
    in:
      colore_config: colore_config
    out: [colore_catalogs]

  fastcat:
    run: cwl/fastcat.cwl
    in:
      fastcat_config: fastcat_config
      colore_config: colore_config
      colore_catalogs: [colore_catalogs]
    out: [fastcat_catalogs]

  NaMaster:
    run: tools/namaster.cwl
    in:
      namaster_config: namaster_config
      fastcat_catalogs: [fastcat_catalogs]
    out: [sacc]

  LikeAnalysis:
    run: tools/like.cwl
    in:
      likelihood_config: likelihood_config
      sacc: [sacc]
    out: [Omega] 
    

