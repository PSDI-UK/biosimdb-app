#!/usr/bin/env python
from flask import Blueprint

#bp = Blueprint('main', __name__, url_prefix='/main')
bp = Blueprint('main', __name__)

from . import head
