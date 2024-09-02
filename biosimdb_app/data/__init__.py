#!/usr/bin/env python
from flask import Blueprint

#bp = Blueprint('data', __name__, url_prefix='/data')
bp = Blueprint('data', __name__)

from . import dataViewer
from . import dataSearcher
from . import showprovenance