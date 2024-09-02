#!/usr/bin/env python
from . import bp
from flask import render_template, current_app
import os

@bp.route('/')
def index():
    return render_template('content/home.html', page_name="/")

@bp.route('/contact')
def contact():
    return render_template('content/contact.html', page_name="/contact")

@bp.route('/help')
def help():
    return render_template('content/help.html', page_name="/help")

@bp.route('/molstar')
def molstar():
    # https://github.com/molstar/molstar/blob/b1b19726842f9a266ffe31b3cd66a4f7536496d5/src/apps/viewer/app.ts#L467
    # top_name = os.path.join(current_app.instance_path, "simulations", "1AKI_minimised.gro")
    # traj_name = os.path.join(current_app.instance_path, "simulations", "1AKI_minimised.trr")
    # pdb_name = os.path.join(current_app.instance_path, "simulations", "7lcj.pdb")

    # with open(top_name, 'r') as file:
    #     top_data = file.read()
    # with open(traj_name, 'rb') as file:
    #     traj_data = file.read()
    # with open(pdb_name, 'r') as file:
    #     pdb_data = file.read()

    return render_template('data/molstar.html', page_name="/molstar")
