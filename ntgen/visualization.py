import matplotlib.pyplot as plt
import networkx as nx


def draw(
    structure: str,
    sequence: str | None = None,
    colors: list | None = None,
    node_size: int = 10,
    ax: plt.Axes | None = None,
) -> None:
    """
    Draw a structure
    """

    compl = str.maketrans('ATUCG ', 'UAAGC ')

    if sequence is None:
        sequence = [''] * len(structure)

    if colors is None:
        colors = ['skyblue'] * len(structure)

    assert len(structure) == len(sequence)

    hide_axis = False
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
        hide_axis = True

    # Place nucleotides in a graph with pairings obtained from the structure
    G = nx.Graph()
    stack = []
    for i, (char, nt, color) in enumerate(zip(structure, sequence, colors, strict=True)):
        G.add_node(i, label=nt, color=color)
        if char == '(':
            stack.append(i)
        elif char == ')':
            j = stack.pop()
            G.add_edge(j, i)

        if i > 0:
            G.add_edge(i - 1, i)

    # This layout seems to work well for gRNAs
    pos = nx.kamada_kawai_layout(G)
    # but we need to rotate the structures
    flipped_pos = {node: (1 - x, y) for node, (x, y) in pos.items()}
    rotated_pos = {node: (-y, x) for node, (x, y) in flipped_pos.items()}

    norm_seq = [s.translate(compl) for s in sequence]
    edge_colors = []
    for n1, n2 in G.edges():
        if abs(n1 - n2) == 1:  # subsequent bases
            edge_colors.append('lightgray')
        else:
            edge_colors.append('red' if norm_seq[n1] != norm_seq[n2].translate(compl) else 'lightgray')

    nx.draw(
        G,
        rotated_pos,
        with_labels=False,
        node_size=node_size,
        node_color=[G.nodes[n]['color'] for n in G.nodes],
        edge_color=edge_colors,
        hide_ticks=False,
        ax=ax)
    nx.draw_networkx_labels(
        G,
        rotated_pos,
        labels={n: G.nodes[n]['label'] for n in G.nodes},
        font_size=node_size / 15,
        hide_ticks=False,
        ax=ax
    )

    if hide_axis:
        ax.set_axis_off()
