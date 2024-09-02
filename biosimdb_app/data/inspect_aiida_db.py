#!/usr/bin/env python
"""
Inspect AiiDA databases.
"""

from aiida import orm, profile_context
from aiida.storage.sqlite_zip.backend import SqliteZipBackend
from aiida.tools.visualization import Graph
from pyvis.network import Network
import tempfile
import os


def aiida_attributes(archive_filename, start, end):
    """
    Returns the aiida attributes for each process node.
    """
    #  create a profile instance from the archive path
    archive_profile = SqliteZipBackend.create_profile(archive_filename)


    # As per a standard profile, we can now use the QueryBuilder, 
    # to find and query for data
    attributes_list = []
    with profile_context(archive_profile):
        processes = orm.QueryBuilder().append(orm.ProcessNode).all(flat=True)
        # processes = orm.QueryBuilder().append(orm.ProcessNode).add_filter(orm.ProcessNode, {'id': {'>': 41}}).all(flat=True)
        n_processes = len(processes)
        print("******", processes)

        graph, html_content, options = convert_aiida_graph(processes, start, end)
        for process in processes:
            attributes_list.append(process.attributes)
            # input_nodes = process.get_incoming().all_nodes()
            # output_nodes = process.get_outgoing().all_nodes()
        
        '''graph = Graph(graph_attr={"size": "8!,8!", "rankdir": "TB"})
        for process in processes:
            print(f"process node: {process}")
            graph.add_incoming(process, annotate_links="both")
            graph.add_outgoing(process, annotate_links="both")
            # node = load_node(process.uuid)
            #filename = node.filename
            # print(node)
            print("Attributes:")
            for key, value in process.attributes.items():
                print(f"\t{key}: {value}")
            print("Extras:")
            for key, value in process.extras.items():
                print(f"\t{key}: {value}")
            last_job_info = process.get_last_job_info()
            print("Last Job Information:")
            print(last_job_info)
            # for key, value in last_job_info.items():
            #     print(f"\t{key}: {value}")
            print("incoming nodes")
            input_nodes = process.get_incoming().all_nodes()
            for input_node in input_nodes:
                print(f"\tinput: {input_node, input_node.node_type}")
                print()
                if isinstance(input_node, List):
                    print(input_node.get_list())
                if isinstance(input_node, Str):
                    print(input_node.value)
            print("outgoing nodes")
            output_nodes = process.get_outgoing().all_nodes()
            for output_node in output_nodes:
                print(f"\toutput: {output_node, output_node.node_type}")
                print()
                if output_node in [orm.List, orm.Str]:
                    print(output_node.value)'''
    return attributes_list, graph, html_content, options, n_processes


def convert_aiida_graph(processes, start, end):
    graph = Graph(graph_attr={"size": "8!,8!", "rankdir": "TB"},
                  #global_node_style={"shape": "circle"},
                  )
    for process in processes[start:end]:
        graph.add_incoming(process) #, annotate_links="both") # label, type or both
        graph.add_outgoing(process) #, annotate_links="both")

    # graph.add_incoming(processes[start]) #, annotate_links="both") # label, type or both
    # graph.add_outgoing(processes[start]) #, annotate_links="both")
    # Create a graph with graphviz from aiida graph
    graph = graph.graphviz
    graph = graph.source # string of dot file contents
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file.write(graph)
        temp_file_path = temp_file.name
    # Create interactive pyvis graph from graphviz graph
    # network = Network(height='500px', width='100%')
    network = Network(directed=True)
    # Test example of a graph
    # l = 'testhsjdvfhdafviufheoafelriafbirafsjdflsfvdsfvksfvusygfureilfgreifhreiu'
    # network.add_node(0, group=0, label=l, shape='ellipse')
    # network.add_node(1, group=0, label=l, shape='ellipse')
    # network.add_node(2, group=1, label=l, shape='ellipse')
    # network.add_node(3, group=1, label=l, shape='ellipse')
    # network.add_node(4, group=2, label=l, shape='ellipse')
    # network.add_node(5, group=2, label=l, shape='ellipse')
    # network.add_edges([(0, 1), (2,1), (3,4), (3,2), (4,5), (1,4), (1,5)])

    # network.show_buttons(filter_=['physics']) # only works foe example graph
    # network.show_buttons()
    network.from_DOT(temp_file_path)
    options = network.options
    options_code = '''
        var options = {
            "physics": {
                //"enabled": false,
                "barnesHut": {
                //"springConstant": 0.05,
                //"avoidOverlap": 0.0
                },
                //"minVelocity": 0.75
            },
            "height": "100%",
            "width": "100%",
            "interaction": {
                "hover": true,
                //"navigationButtons": true,
                "zoomView": false
            }
        };
    '''

    os.remove(temp_file_path)
    html_content = network.generate_html()
    # if using from_DOT, need to manually update options as set_options doesn't work
    html_content = html_content.replace("var options = parsedData.options;", 
                    options_code)
    # html_content = html_content.replace("ellipse", "circle")
    return graph, html_content, options
