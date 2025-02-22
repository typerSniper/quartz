# Ref: https://github.com/facebookresearch/hydra/blob/1.2_branch/examples/tutorials/structured_configs/5.2_structured_config_schema_different_config_group/database_lib.py
from dataclasses import dataclass, field
from typing import Any, List

import hydra
from config.base_config import BaseConfig
from config.ibm2_config import *
from config.ibm_config import *
from config.ionq_config import *
from config.nam2_config import *
from config.nam_config import *
from config.rig_config import *
from config.tdg_config import *
from hydra.core.config_store import ConfigStore
from omegaconf import MISSING, OmegaConf  # Do not confuse with dataclass.MISSING

defaults = [
    # config group name c will load config named base
    {"c": "base"}
]


@dataclass
class Config:
    # this is unfortunately verbose due to @dataclass limitations
    defaults: List[Any] = field(default_factory=lambda: defaults)
    # Hydra will populate this field based on the defaults list
    c: BaseConfig = MISSING


cs = ConfigStore.instance()
cs.store(name="config", node=Config)
cs.store(group="c", name="base", node=BaseConfig)
cs.store(group="c", name="convert", node=ConvertConfig)

# other config groups for training
cs.store(group="c", name="nam", node=NamConfig)
cs.store(group="c", name="nam_ft", node=NamFTConfig)
cs.store(group="c", name="nam_pret", node=NamPretrainConfig)
cs.store(group="c", name="nam_mp", node=NamMPConfig)
cs.store(group="c", name="nam_rm_mp", node=NamRMMPConfig)
cs.store(group="c", name="nam_test", node=NamTestConfig)
cs.store(group="c", name="nam_alltest", node=NamAllTestConfig)
cs.store(group="c", name="nam_convert", node=NamConvertConfig)

cs.store(group="c", name="nam2", node=Nam2Config)
cs.store(group="c", name="nam2_ft", node=Nam2FTConfig)
cs.store(group="c", name="nam2_mp", node=Nam2MPConfig)
cs.store(group="c", name="nam2_rm_mp", node=Nam2RMMPConfig)
cs.store(group="c", name="nam2_test", node=Nam2TestConfig)

cs.store(group="c", name="ibm", node=IBMConfig)
cs.store(group="c", name="ibm_ft", node=IBMFTConfig)
cs.store(group="c", name="ibm_pret", node=IBMPretrainConfig)
cs.store(group="c", name="ibm_mp", node=IBMMPConfig)
cs.store(group="c", name="ibm_test", node=IBMTestConfig)

cs.store(group="c", name="ibm2", node=IBM2Config)
cs.store(group="c", name="ibm2_ft", node=IBM2FTConfig)
cs.store(group="c", name="ibm2_pret", node=IBM2PretrainConfig)
cs.store(group="c", name="ibm2_mp", node=IBM2MPConfig)
cs.store(group="c", name="ibm2_test", node=IBM2TestConfig)

cs.store(group="c", name="tdg", node=TdgConfig)
cs.store(group="c", name="tdg_ft", node=TdgFTConfig)
cs.store(group="c", name="tdg_mp", node=TdgMPConfig)
cs.store(group="c", name="tdg_rm_mp", node=TdgRMMPConfig)
cs.store(group="c", name="tdg_test", node=TdgTestConfig)

cs.store(group="c", name="rig", node=RigConfig)
cs.store(group="c", name="rig_ft", node=RigFTConfig)
cs.store(group="c", name="rig_mp", node=RigMPConfig)
cs.store(group="c", name="rig_rm_mp", node=RigRMMPConfig)

cs.store(group="c", name="ionq", node=IonQConfig)
cs.store(group="c", name="ionq_ft", node=IonQFTConfig)
cs.store(group="c", name="ionq_mp", node=IonQMPConfig)

# cfg groups for test
# cs.store(group="c", name="nam_test", node=NamMultiPretrainConfig)
