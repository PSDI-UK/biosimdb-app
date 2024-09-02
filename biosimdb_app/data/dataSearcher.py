#!/usr/bin/env python
from . import bp
from flask import render_template, request, redirect, url_for
# from extensions import mysql
import os
from . import dataViewer
# from biosimdb_app.utils.db import get_db


@bp.route('/search_results', methods=['GET', 'POST'])
def get_search_text():
    """
    Return search results
    """
    if request.method == "POST":
        search_text = request.form["search_text"]
        if search_text:
            # search_text = request.args.get('search_text', search_text, type=str)
            return redirect((url_for('data.search_data', search_text=search_text)))
        else:
            return redirect((url_for('data.display_data')))


    

@bp.route('/search_results/<search_text>')
def search_data(search_text):
    query = f"""SELECT title, abstract, datetime_added, 
            simulation_ID, topology_file 
            FROM
            `simulation`
            WHERE
            (`title` AND `software` AND `forcefields` AND 
            `abstract` AND `citations`) LIKE '{search_text}%'
            """

    records, page, total_pages, sort_column, sort_direction = dataViewer.get_data(query)

    return render_template('data/search_results.html', search_text=search_text, 
            records=records, page_name=f"/search_results/{search_text}", page=page, 
            total_pages=total_pages, sort_column=sort_column, 
            sort_direction=sort_direction)