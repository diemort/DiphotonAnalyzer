# DiphotonAnalyzer

Simple EDAnalyzer for flashgg output

### Prerequisites
- a valid installation of [`flashgg`](https://github.com/cms-analysis/flashgg) to define all useful data formats; follow the installation procedure recommended for your release.
- [PPS proton reconstruction](https://github.com/forthommel/cmssw/tree/proton_reco_9_4_X)
  > run `git cms-merge-topic forthommel:proton_reco_9_4_X` and recompile.

### (Simplified set of) instructions
- Clone this directory inside your $CMSSW_BASE/src directory containing the flashgg installation (with cmsenv already run)
- `scram b`
- Edit [your input configuration file](TreeProducer/test/treeProducer_data_cfg.py) accordingly (cuts, input files, any future feature...)
- `cmsRun your_input_configuration_file_cfg.py`
