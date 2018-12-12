# DiphotonAnalyzer

Simple EDAnalyzer for flashgg output

### Prerequisites
- [`flashgg`](https://github.com/cms-analysis/flashgg)
- [PPS proton reconstruction algorithm](https://github.com/CTPPS/cmssw/tree/ctpps_initial_proton_reconstruction_CMSSW_9_4_11)
> run `git cms-merge-topic forthommel:ctpps_initial_proton_reconstruction_CMSSW_9_4_9` before your initial `flashgg` build
- ...

### (Simplified set of) instructions
- Clone this directory inside your $CMSSW_BASE/src directory containing the flashgg installation (with cmsenv already run)
- `scram b`
- Edit <TreeProducer/test/TreeProducer_cfg.py> accordingly (cuts, input files, any future feature...)
- `cmsRun TreeProducer/test/TreeProducer_cfg.py`
