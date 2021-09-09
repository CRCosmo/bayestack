#!/usr/bin/env python

import lumfuncUtils
import bayestack
import os,sys
import importlib

__name_cached=__name__
if __name__=='__main__':
    param_file=sys.argv[-1]
    settingsf=param_file.split('.')[-2]
    set_module=importlib.import_module(settingsf)
    globals().update(set_module.__dict__)
__name__=__name_cached


bayestack.run()

