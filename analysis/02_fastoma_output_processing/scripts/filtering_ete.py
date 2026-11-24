"""ETE 4 operations and visual checks for Filtering_orthoxml.ipynb."""
import hashlib
import importlib.util

from ete4 import Tree
from lxml import etree
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

import sister_nhx as w


def read(path):
    return Tree(path.read_text(), parser=1)


def convert_and_check(r):
    """Feed each complete group to the same converter used by Irene's to-nhx command."""
    path = r.tools / "orthoxml/converters/to_nhx.py"
    spec = importlib.util.spec_from_file_location("irene_to_nhx", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    r.raw.mkdir()
    seen = set()
    registry_hash = hashlib.sha256()
    tips = 0
    with tqdm(total=len(r.selected), desc="Convert and verify HOGs") as bar:
        for node in w.records(r.xml):
            if node.tag == w.SPECIES:
                species = node.attrib["name"]
                for gene in node.iter(w.GENE):
                    gid = gene.attrib["id"]
                    w.require(gid not in r.registry, f"Duplicate registry ID {gid}")
                    r.registry[gid] = (gene.attrib["protId"], species)
                    w.registry_digest(registry_hash, gene, species)
                continue
            hog = node.attrib["id"]
            w.require(hog not in seen and w.tree_hash(node) == r.selected[hog], f"Incomplete HOG {hog}")
            seen.add(hog)
            converter = module.OrthoxmlToNewick(xref_tag="id", encode_levels_as_nhx=True)
            for event, element in etree.iterwalk(node, events=("start", "end")):
                if event == "start":
                    if element.tag == w.REF:
                        gid = element.attrib["id"]
                        converter.gene2xref[gid] = (gid, r.registry[gid][1])
                        tips += 1
                    converter.start(element.tag, element.attrib)
                else:
                    converter.end(element.tag)
            dest = w.tree_path(r.raw, hog)
            dest.write_text(converter.close()[hog])
            verify(read(dest), node, r.registry, hog)
            bar.update()
    w.require(seen == r.selected.keys() and tips == len(r.registry) == r.manifest["subset_registry_size"],
              "Subset roots or registry differ from source")
    w.require(registry_hash.hexdigest() == r.manifest["source_registry_sha256"], "Protein registry differs")
    r.checkpoint("Converted complete HOGs and verified all written NHX trees against XML")


def verify(tree, xml, registry, hog):
    expected = {"Ev": "0>1" if xml.tag in (w.ORTHO, w.REF) else "1>0"}
    if xml.tag == w.REF:
        gid = xml.attrib["id"]
        w.require(tree.is_leaf and tree.name == gid, f"{hog}: leaf identity differs {gid}")
        expected["S"] = registry[gid][1]
    else:
        for prop in xml.iterchildren(w.PROP):
            key, value = prop.attrib["name"], prop.attrib["value"]
            if key == "TaxRange":
                value = value.replace("'", "").replace("T=", "")
                if value:
                    expected["S"] = value
            elif key.lower() in {"taxid", "taxonid", "taxon_id", "ncbi_taxon_id"} and value:
                expected["T"] = value
        children = [child for child in xml if child.tag in w.TREE]
        for child, original in zip(tree.children, children, strict=True):
            verify(child, original, registry, hog)
    w.require({k: v for k, v in tree.props.items() if k in ("Ev", "S", "T")} == expected,
              f"{hog}: annotation differs {tree.props} vs {expected}")


def groups(tree, xml):
    # XML distinguishes fragments, which Irene's NHX labels as duplication-like.
    internal = (n for n in tree.traverse("preorder") if not n.is_leaf)
    for n, original in zip(internal, xml.iter(*w.GROUP), strict=True):
        event = "speciation" if n.props["Ev"] == "0>1" else w.xml_event(original)
        yield event, [c.name for c in n.children if c.is_leaf]


def prune(tree, removed, stats, expected):
    """Contract from the leaves upwards, retaining original daughter order."""
    def visit(node):
        if node.is_leaf:
            return None if node.name in removed else node
        children = list(node.children)
        node.remove_children()
        for child in children:
            retained = visit(child)
            if retained is not None:
                node.add_child(retained)
        if not node.children:
            stats["empty"] += 1
            return None
        if len(node.children) == 1:
            stats["unary"] += 1
            only = node.children[0].detach()
            if node is tree or node.dist is None or only.dist is None:
                only.props.pop("dist", None)
            else:
                only.dist += node.dist
            return only
        return node

    root = visit(tree)
    w.require(root is not None, "Training HOG became empty")
    remaining = list(root.leaf_names())
    w.require(len(remaining) == len(expected) and set(remaining) == expected, "Incorrect retained leaf set")
    w.require(all(len(n.children) != 1 for n in root.traverse()), "Unary node survived")
    return root, remaining


def signature(node):
    return node.props, tuple(signature(c) for c in node.children)


def shape(node):
    if node.is_leaf:
        return node.name
    return ["S" if node.props["Ev"] == "0>1" else "D", *(shape(c) for c in node.children)]


def write_checked(root, dest):
    dest.write_text(root.write(parser=1, props=["Ev", "S", "T"], format_root_node=True))
    w.require(signature(root) == signature(read(dest)), f"{dest.name}: annotated round trip failed")


def save_figure(r, fig, name):
    fig.tight_layout()
    fig.savefig(r.out / name, dpi=150, bbox_inches="tight")
    plt.show()
    plt.close(fig)


def plot_cohort(r):
    census = pd.read_csv(r.out / "root_census.tsv", sep="\t")
    subset = census.loc[census.train > 0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    bars = axes[0].bar(["With training", "Without training"], [len(subset), len(census)-len(subset)])
    axes[0].bar_label(bars, fmt="{:,.0f}")
    axes[0].set_ylabel("Root HOGs")
    hist = subset.train.value_counts().sort_index()
    axes[1].scatter(hist.index, hist.values, s=14)
    axes[1].set(xscale="log", yscale="log", xlabel="Training tips per selected HOG", ylabel="Root HOGs")
    save_figure(r, fig, "01_training_cohort.png")


def plot_results(r):
    counts = pd.read_csv(r.out / "hog_counts.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    kinds = ["train", "oma", "test"]
    x = np.arange(3)
    for i, prefix in enumerate(["before", "after"]):
        vals = [counts[f"{prefix}_{k}"].sum() for k in kinds]
        axes[0].bar(x + (i-.5)*.4, vals, width=.4, label=prefix.title())
    axes[0].set_xticks(x, ["Training", "OMA", "Test"])
    axes[0].set(ylabel="Tips (including singleton HOGs)", yscale="symlog")
    axes[0].legend()
    fates = counts.status.value_counts().reindex(["multi_tip","singleton"], fill_value=0)
    bars = axes[1].bar(["Multi-tip", "Singleton"], fates)
    axes[1].bar_label(bars, fmt="{:,.0f}")
    axes[1].set_ylabel("Root HOGs")
    axes[2].scatter(counts.before_tips, counts.after_tips, alpha=.3, s=8)
    upper = counts.before_tips.max()
    axes[2].plot([1,upper],[1,upper], color="grey", ls="--")
    axes[2].set(xscale="log", yscale="log", xlabel="Tips before pruning", ylabel="Tips after pruning")
    save_figure(r, fig, "03_pruning_outcomes.png")
    return counts


def plot_example(r, counts):
    small = counts.query("status == 'multi_tip' and before_tips <= 24 and removed_oma + removed_test > 0")
    if small.empty:
        print("No changed multi-tip HOG with at most 24 leaves.")
        return
    row = small.sort_values(["removed_oma","removed_test","unary_collapsed"], ascending=False).iloc[0]
    hog = row.root_hog
    trees = [read(w.tree_path(folder, hog)) for folder in (r.raw, r.final)]
    kept = set(trees[1].leaf_names())
    fig, axes = plt.subplots(1, 2, figsize=(13, max(4,row.before_tips*.33)))
    for tree, ax, title in zip(trees, axes, ["Original", "Pruned"]):
        ys = {n: i for i,n in enumerate(reversed(list(tree.leaves())))}
        depth = {tree: 0}
        for n in tree.traverse("preorder"):
            for c in n.children:
                depth[c] = depth[n]+1
        for n in tree.traverse("postorder"):
            if n.is_leaf:
                kind = w.source(r, n.name)
                label = f"{r.registry[n.name][0]}  [{kind}]"
                ax.text(depth[n]+.08, ys[n], label, va="center", fontsize=8,
                        color="#222222" if n.name in kept else "#a23426")
            else:
                cy = [ys[c] for c in n.children]
                ys[n] = (min(cy)+max(cy))/2
                ax.vlines(depth[n], min(cy), max(cy), color="0.45", lw=1)
                for c in n.children:
                    ax.hlines(ys[c], depth[n], depth[c], color="0.45", lw=1)
                ax.text(depth[n],ys[n],"S" if n.props["Ev"]=="0>1" else "D",
                        fontsize=8, ha="center", va="center", backgroundcolor="white")
        ax.set_xlim(-.2,max(depth.values())+4)
        ax.set_ylim(-1,len(ys) if tree.is_leaf else len(list(tree.leaves())))
        ax.set_title(title)
        ax.axis("off")
    fig.suptitle(f"{hog}   •   red = removed tip   •   S / D = speciation-like / duplication-like")
    save_figure(r,fig,"04_example_tree.png")
    return hog


def plot_conversion(r):
    census = pd.read_csv(r.out / "root_census.tsv", sep="\t")
    subset = census[census.selected].copy()
    sizes = subset.oma + subset.train + subset.test
    fig, ax = plt.subplots(figsize=(7, 3.5))
    bins = np.unique(np.geomspace(1, max(2, sizes.max()+1), 45).astype(int))
    ax.hist(sizes, bins=bins, color="#377b88")
    ax.set(xscale="log", yscale="log", xlabel="Tips per verified unpruned tree", ylabel="Root HOGs")
    save_figure(r, fig, "02_verified_tree_sizes.png")


def provenance(r):
    from importlib.metadata import version
    import json
    import sys
    from pathlib import Path
    for path in [Path(__file__), Path(w.__file__), r.extractor.parent/"parser.py", r.extractor.parent/"utils.py"]:
        r.manifest["code_sha256"][r.rel(path)] = w.sha(path)
    r.manifest["packages"].update({p: version(p) for p in ("numpy","pandas","matplotlib","tqdm")})
    r.manifest["converter"] = "Irene OrthoxmlToNewick(xref_tag='id', encode_levels_as_nhx=True), per root"
    r.manifest["pruning_backend"] = "ETE 4"
    r.manifest["python_executable"] = sys.executable
    cells = json.loads((r.root/r.manifest["notebook"]).read_text())["cells"]
    r.manifest["notebook_cells_sha256"] = hashlib.sha256(json.dumps(
        [(c["cell_type"], "".join(c["source"])) for c in cells]).encode()).hexdigest()
    r.checkpoint("ETE workflow inputs recorded")
