# Requirements

Before [Installation](https://github.com/msk-access/nucleo/tree/master/docs/installation-and-usage.md) of the pipeline, make sure your system supports these requirements
=======


Following are the requirements for running the workflow:

* A system with either [docker](https://www.docker.com/) or [singularity](https://sylabs.io/docs/) configured.
* Python 3.6 \(for running [cwltool](https://github.com/common-workflow-language/cwltool)and running [toil-cwl-runner](https://toil.readthedocs.io/en/latest/running/introduction.html)\)
  * Python Packages \(will be installed as part of pipeline installation\):
    * toil\[cwl\]==5.1.0
    * pytz==2021.1
    * typing==3.7.4.3
    * ruamel.yaml==0.16.5
    * pip==20.2.3
    * bumpversion==0.6.0
    * wheel==0.35.1
    * watchdog==0.10.3
    * flake8==3.8.4
    * tox==3.20.0
    * coverage==5.3
    * twine==3.2.0
    * pytest==6.1.1
    * pytest-runner==5.2
    * coloredlogs==10.0
    * pytest-travis-fold==1.3.0
  * Python Virtual Environment using [virtualenv](https://virtualenv.pypa.io/) or [conda](https://docs.conda.io/en/latest/).

