import cProfile
import time


def draw_call_graph(func_to_profile, path=None):
    """Use PyCallGraph [1] and GraphViz [2] to draw the call graph
    for the function specified in func_to_profile.
    [1]: https://github.com/gak/pycallgraph, https://pycallgraph.readthedocs.io/en/develop/index.html
    [2]: http://www.graphviz.org/about/

    Saves the call graph to path."""

    from pycallgraph import PyCallGraph, Color, Config
    from pycallgraph.output import GraphvizOutput

    def orange_green(node):
        """Make a higher total time have an orange colour
        and a higher number of calls have a green colour.

        """

        return Color(0.2 + node.time.fraction * 0.8,
                     0.2 + node.calls.fraction * 0.4 + node.time.fraction * 0.4,
                     0.2,)

    graphviz = GraphvizOutput()
    pycallgraph = PyCallGraph(output=graphviz, config=Config(include_stdlib=False))
    graphviz.edge_color_func = lambda e: Color(0, 0, 0)
    graphviz.node_color_func = orange_green
    graphviz.output_file = path if path else 'pycallgraph.png'
    graphviz.done()

    with PyCallGraph(output=graphviz):
        func_to_profile()


def run_cprofiler(str_to_run, sort_type='time'):
    """Use cProfile to profile function specified as a string
    in func_to_profile. Sort outcomes by sort_type
    (see cProfile for viable sort types).

    Outputs results to console."""

    print('Profiling starting...')
    cProfile.run(str_to_run, sort=sort_type)
    print('Profiling complete.')


def run_timed(func_to_run):
    """Run function and print execution time to console.

    Outputs results to console."""

    print('Starting...')
    start_time = time.time()
    func_to_run()
    print('Complete. Execution time: %s seconds' %
          (time.time() - start_time))
