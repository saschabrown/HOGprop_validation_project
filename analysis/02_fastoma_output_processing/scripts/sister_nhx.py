"""Stage 02 operations for the complete-HOG sister NHX notebook."""
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import platform
import re
import subprocess
import sys
from collections import Counter, defaultdict
from contextlib import ExitStack
from datetime import datetime, timezone
from fractions import Fraction
from importlib.metadata import version
from pathlib import Path
from uuid import uuid4

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from lxml import etree

NS = "{http://orthoXML.org/2011/}"
ORTHO, PARA, REF, GENE, SPECIES, PROP = (
    NS + t for t in ("orthologGroup", "paralogGroup", "geneRef", "gene", "species", "property")
)
GROUP = (ORTHO, PARA)
TREE = (*GROUP, REF)
# Source fingerprint includes Irene's modified event encoding, which the release number alone cannot identify.
IRENE_SOURCE = "1c918112c561f61bd758e1af5256b78c3f54d11597e3dfe253b337efe7cb239f"
NOTEBOOK = "cafa_train_hog_sisters_nhx.ipynb"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def code_hash(source):
    h = hashlib.sha256()
    for file in sorted(source.rglob("*.py")):
        h.update(file.relative_to(source).as_posix().encode() + b"\0" + bytes.fromhex(sha(file)))
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def table(stack, path, columns):
    handle = stack.enter_context(path.open("w", newline=""))
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(columns.split())
    return writer.writerow


class Run:
    def __init__(self, root, mode, notebook_name=NOTEBOOK, ete=False):
        require(root.is_absolute(), "WP1_ROOT must be absolute")
        require(mode in ("stage01", "fixture"), "Mode must be stage01 or fixture")
        self.mode = mode
        self.root = root
        repo = root / "HOGprop_validation_project"
        self.stage = repo / "analysis/02_fastoma_output_processing"
        upstream = repo / "analysis/01_fastoma/results/20261105_101624_7305523"
        inputs = self.stage / "fixtures/sister_nhx" if mode == "fixture" else upstream
        self.source = inputs / "FastOMA_HOGs.orthoxml"
        self.fastas = inputs / "RootHOGsFasta"
        self.mapping = (inputs if mode == "fixture" else repo) / "irene_cafa_oma_fasta_pairs.tsv"
        self.extractor = self.stage / ("OMA_orthoxml/orthoxml_subset.py" if ete else "scripts/orthoxml_subset.py")
        self.tools = next((p for p in sorted((root / "local_config").glob(
            "irene_tree_reference.*/orthoxml_tools_src")) if code_hash(p) == IRENE_SOURCE), None)
        require(self.tools is not None, "Expected the pinned Irene orthoxml_tools_src in local_config/irene_tree_reference.*/")
        for path in (self.mapping, self.source, self.extractor):
            require(path.is_file(), f"Missing input: {path}")
        require(self.fastas.is_dir(), f"Missing Stage 01 sequence directory: {self.fastas}")
        run_id = os.getenv("WP1_STAGE02_RUN_ID") or f"sister_nhx_{mode}_{uuid4().hex[:12]}"
        require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id), "Invalid run ID")
        self.out = self.stage / "results" / run_id
        self.out.mkdir(parents=True, exist_ok=False)
        self.xml = self.out / "complete_training_hogs.orthoxml"
        self.raw = self.out / "unpruned"
        self.final = self.out / "pruned"
        self.final.mkdir()
        self.train = {}
        self.test = set()
        self.selected = {}
        self.registry = {}
        self.assigned = {}
        self.coverage = {}
        self.taxon = {}
        with self.mapping.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                code, tax = row["oma_code"], row["ncbi_taxid"]
                require(code and code not in self.taxon and tax.isdigit(), f"Invalid mapping: {code}, {tax}")
                self.taxon[code] = tax
        notebook = self.stage / "notebooks" / notebook_name
        cells = json.loads(notebook.read_text())["cells"]
        sources = [(c["cell_type"], "".join(c["source"])) for c in cells]
        self.manifest = dict(
            status="incomplete", started=datetime.now(timezone.utc).isoformat(), mode=mode,
            run_directory=self.rel(self.out), platform=platform.platform(), python=sys.version,
            packages={p: version(p) for p in ("lxml", "biopython", ("ete4" if ete else "dendropy"), "ipykernel")},
            irene_source=self.rel(self.tools), irene_python_source_sha256=IRENE_SOURCE,
            inputs_sha256={self.rel(p): sha(p) for p in (self.source, self.mapping)},
            code_sha256={self.rel(p): sha(p) for p in (self.extractor, Path(__file__))},
            notebook=self.rel(notebook), notebook_cells_sha256=hashlib.sha256(
                json.dumps(sources).encode()).hexdigest(),
            sequence_inputs={}, random_seed=None,
            parameters=dict(focal="every registered training tip in a selected root HOG",
                            eligible="direct OMA sibling, mapped taxon equal, orthologGroup parent",
                            matrix="BLOSUM62", mode="global", open_gap=-11, extend_gap=-1,
                            identity="first optimal alignment; exact matches / alignment columns",
                            tie="lowest numeric XML gene ID", missing_sequence="no removal",
                            test="remove all", unary="collapse all", minimum_final_tips=2),
        )
        if (repo / ".git").exists():
            for key, args in (("git_commit", ["rev-parse", "HEAD"]),
                              ("git_status", ["status", "--short"])):
                self.manifest[key] = subprocess.check_output(["git", "-C", repo, *args], text=True).strip()
        self.checkpoint("inputs recorded")

    def rel(self, path):
        return str(path.relative_to(self.root))

    def checkpoint(self, stage):
        self.manifest["last_completed_step"] = stage
        write_json(self.out / "provenance.json", self.manifest)


def records(path):
    # Retain one complete species registry or root HOG, clearing preceding siblings.
    for _, node in etree.iterparse(str(path), events=("end",), tag=(SPECIES, *GROUP),
                                   remove_blank_text=True, resolve_entities=False, no_network=True):
        parent = node.getparent()
        if parent.tag not in (NS + "orthoXML", NS + "groups"):
            continue
        yield node
        node.clear()
        while node.getprevious() is not None:
            del parent[0]


def registry_digest(digest, gene, species):
    for value in (gene.attrib["id"], gene.attrib["protId"], species):
        digest.update(value.encode() + b"\0")


def tree_hash(node):
    return hashlib.sha256(etree.tostring(node, method="c14n")).hexdigest()


def select_and_extract(r):
    """Scan original XML and feed complete root IDs and complete reference sets to Irene's filter."""
    refs = set()
    source_roots = set()
    with ExitStack() as stack:
        census = table(stack, r.out / "root_census.tsv", "root_hog oma train test selected")
        for node in records(r.source):
            if node.tag == SPECIES:
                species = node.attrib["name"]
                if species.isdigit():
                    for gene in node.iter(GENE):
                        gid, protein = gene.attrib["id"], gene.attrib["protId"]
                        _, taxon, split = protein.rsplit("_", 2)
                        require(taxon == species and split in ("train", "test"), f"Invalid CAFA label {protein}")
                        require(gid not in r.train and gid not in r.test, f"Duplicate CAFA gene {gid}")
                        if split == "train":
                            r.train[gid] = (protein, species)
                        else:
                            r.test.add(gid)
                continue
            hog = node.attrib["id"]
            require(hog not in source_roots, f"Duplicate root HOG {hog}")
            source_roots.add(hog)
            ids = [ref.attrib["id"] for ref in node.iter(REF)]
            training = [i for i in ids if i in r.train]
            tests = sum(i in r.test for i in ids)
            census([hog, len(ids)-len(training)-tests, len(training), tests, bool(training)])
            if training:
                require(node.tag == ORTHO and re.fullmatch(r"[\w.:-]+", hog), f"Unsupported root {hog}")
                require(len(set(ids)) == len(ids) and refs.isdisjoint(ids), f"Repeated selected gene in {hog}")
                r.selected[hog] = tree_hash(node)
                refs.update(ids)
                r.assigned.update(dict.fromkeys(training, hog))
    (r.out / "root_hog_ids.txt").write_text("".join(h+"\n" for h in r.selected))
    # Same two filter inputs as extract_subset_from_orthoxml. Gathering refs during the
    # census avoids its extra pass retaining the full original XML gene registry in lxml.
    sys.path.insert(0, str(r.tools))
    sys.path.insert(0, str(r.extractor.parent))
    import orthoxml
    require(Path(orthoxml.__file__).parent.parent == r.tools, "Restart kernel to load pinned orthoxml source")
    spec = importlib.util.spec_from_file_location("irene_subset", r.extractor)
    subset = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(subset)
    registry_hash = hashlib.sha256()

    class AuditedFilter(subset.FilterHOGSubsetWithContainedGenes):
        def process_species(self, elem):
            species = elem.attrib["name"]
            for gene in elem.iter(GENE):
                if gene.attrib["id"] in self.gene_ids:
                    registry_digest(registry_hash, gene, species)
            return super().process_species(elem)

    subset.process_stream_orthoxml(
        str(r.source), str(r.xml), parser_cls=AuditedFilter,
        parser_kwargs={"roothog_ids": set(r.selected), "gene_ids": refs})
    r.manifest.update(source_roots=len(source_roots), selected_roots=len(r.selected),
                      subset_registry_size=len(refs), source_registry_sha256=registry_hash.hexdigest(),
                      subset_sha256=sha(r.xml),
                      extractor_entry="orthoxml_subset.FilterHOGSubsetWithContainedGenes via process_stream_orthoxml")
    r.checkpoint("complete-HOG extraction")


def convert(r):
    # Internal IDs are unique even when a protein label is shared. The leaf audit maps every tip.
    command = [str(r.tools / "cli.py"), "to-nhx", "--infile", str(r.xml), "--outdir", str(r.raw),
               "--encode-levels", "--xref-tag", "id"]
    subprocess.run([sys.executable, "-B", *command], check=True)
    r.manifest["converter_arguments"] = [r.rel(Path(a)) if a.startswith("/") else a for a in command]
    r.checkpoint("NHX conversion")


def read_tree(path):
    import dendropy
    return dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True,
                             suppress_leaf_node_taxa=True, extract_comment_metadata=False)


def tree_path(folder, hog):
    return folder / f"tree_{hog}.nwk"


def xml_event(node):
    if node.tag == ORTHO:
        return "speciation"
    if any(prop.get("name") == "Type" for prop in node.iterchildren(PROP)):
        return "fragment"
    return "duplication"


def verify_nodes(node, xml, registry, hog):
    """Positional comparison preserves multiplicity and distinguishes same-tip unary clades."""
    if xml.tag == REF:
        gid = xml.attrib["id"]
        protein, species = registry[gid]
        require(node.is_leaf() and node.label == gid and protein, f"{hog}: wrong tip for {gid}")
        expected = {"Ev": "0>1", "S": species}
    else:
        children = [c for c in xml if c.tag in TREE]
        require(children and len(children) == node.num_child_nodes(), f"{hog}: child count differs")
        expected = {"Ev": "0>1" if xml.tag == ORTHO else "1>0"}
        for p in xml.iterchildren(PROP):
            name, val = p.attrib["name"], p.attrib["value"]
            if name == "TaxRange":
                level = val.replace("'", "").replace("T=", "")
                if level:
                    expected["S"] = level
            elif name.lower() in ("taxid", "taxonid", "taxon_id", "ncbi_taxon_id") and val:
                expected["T"] = val
        for child, x in zip(node.child_node_iter(), children, strict=True):
            verify_nodes(child, x, registry, hog)
    expected_comment = "&&NHX" + "".join(f":{k}={expected[k]}" for k in ("Ev", "S", "T") if k in expected)
    require(node.comments == [expected_comment], f"{hog}: {node.comments} differs from {expected_comment}")


def validate_unpruned(r):
    expected_files = {tree_path(r.raw, h).name for h in r.selected}
    require({p.name for p in r.raw.iterdir()} == expected_files, "NHX filenames differ from exact root HOG IDs")
    seen = set()
    referenced = 0
    registry_hash = hashlib.sha256()
    for node in records(r.xml):
        if node.tag == SPECIES:
            for gene in node.iter(GENE):
                registry_digest(registry_hash, gene, node.attrib["name"])
                gid = gene.attrib["id"]
                require(gid not in r.registry, f"Duplicate subset registry ID {gid}")
                r.registry[gid] = (gene.attrib["protId"], node.attrib["name"])
            continue
        hog = node.attrib["id"]
        require(hog not in seen and tree_hash(node) == r.selected[hog], f"Incomplete or repeated HOG {hog}")
        seen.add(hog)
        verify_nodes(read_tree(tree_path(r.raw, hog)).seed_node, node, r.registry, hog)
        referenced += sum(1 for _ in node.iter(REF))
    require(seen == r.selected.keys() and len(r.registry) == referenced == r.manifest["subset_registry_size"],
            "Selected roots, reference count or registry size differ")
    require(registry_hash.hexdigest() == r.manifest["source_registry_sha256"],
            "Subset protein IDs or species differ from the Stage 01 gene registry")
    r.checkpoint("all complete roots, NHX tips, topology, events and reported levels verified")


def make_aligner():
    aligner = Align.PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def identity(a, b, aligner):
    if a == b:
        return Fraction(1)
    alphabet = aligner.substitution_matrix.alphabet
    a, b = ("".join(x if x in alphabet else "X" for x in seq) for seq in (a, b))
    alignment = aligner.align(a, b)[0]
    return Fraction(alignment.counts().identities, alignment.length)


def sequences(r, hog, wanted):
    path = r.fastas / (hog.split("_", 1)[0].replace(":", "") + ".fa.gz")
    if not path.is_file():
        r.manifest["sequence_inputs"][r.rel(path)] = "missing"
        return {}
    r.manifest["sequence_inputs"][r.rel(path)] = sha(path)
    seqs = {}
    with gzip.open(path, "rt") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            _, meta, *_ = rec.description.split(maxsplit=2)
            gid = meta.rsplit("||", 1)[-1]
            if gid in wanted:
                require(gid not in seqs and rec.id == r.registry[gid][0], f"{hog}: FASTA identity {gid}")
                seqs[gid] = str(rec.seq).upper()
    return seqs


def choose(focal, options, seqs, aligner):
    if not options:
        return None, {}, "no_eligible_sister"
    if len(options) == 1:
        return options[0], {}, "single"
    if any(not seqs.get(g) for g in (focal, *options)):
        return None, {}, "missing_sequence"
    scores = {g: identity(seqs[focal], seqs[g], aligner) for g in options}
    best = max(scores.values())
    tied = [g for g in options if scores[g] == best]
    return min(tied, key=int), scores, "identity_tie" if len(tied) > 1 else "highest_identity"


def source(r, gid):
    if gid in r.train:
        return "train"
    if gid in r.test:
        return "test"
    return "oma"


def decisions(r, hog, xml, write, aligner, groups=None):
    if groups is None:
        groups = ((xml_event(p), [ref.attrib["id"] for ref in p.iterchildren(REF)])
                  for p in xml.iter(*GROUP))
    parent_info, candidates = {}, {}
    for i, (event, direct_leaves) in enumerate(groups, 1):
        siblings = defaultdict(list)
        for gid in direct_leaves:
            parent_info[gid] = (i, event)
            species = r.registry[gid][1]
            if source(r, gid) == "oma" and species in r.taxon and event == "speciation":
                siblings[r.taxon[species]].append(gid)
        for gid in direct_leaves:
            if gid in r.train:
                candidates[gid] = sorted(siblings[r.train[gid][1]], key=int)
    wanted = {g for focal, opts in candidates.items() if len(opts) > 1 for g in (focal, *opts)}
    seqs = sequences(r, hog, wanted) if wanted else {}
    selected = defaultdict(list)
    for focal, options in candidates.items():
        gid, scores, reason = choose(focal, options, seqs, aligner)
        if gid is not None:
            selected[gid].append(focal)
        write([hog, *parent_info[focal], r.train[focal][1], focal, ",".join(options), gid, reason,
               json.dumps({g: str(s) for g, s in scores.items()})])
    return parent_info, selected


def prune_node(node, removed, stats):
    if node.is_leaf():
        return None if node.label in removed else node
    children = [kept for child in node.child_nodes() if (kept := prune_node(child, removed, stats)) is not None]
    if not children:
        stats["empty"] += 1
        return None
    if len(children) == 1:
        stats["unary"] += 1
        child = children[0]
        # The Stage 01 hierarchy is branchless. Sum only fully measured contracted paths.
        if node.parent_node is None or node.edge_length is None or child.edge_length is None:
            child.edge_length = None
        else:
            child.edge_length += node.edge_length
        return child
    node.set_child_nodes(children)
    return node


def signature(node):
    return (node.label, tuple(node.comments), node.edge_length, tuple(signature(c) for c in node.child_node_iter()))


def filename_set(folder):
    return {p.name for p in folder.iterdir()}


def prune_and_audit(r, ete_backend=None):
    aligner = make_aligner()
    totals = Counter()
    multi = set()
    with ExitStack() as stack:
        assignment = table(stack, r.out / "pruning_assignments.tsv",
                           "root_hog parent_node event taxon cafa_gene_id candidate_oma_gene_ids selected_oma_gene_id decision identity_by_oma_gene")
        audit = table(stack, r.out / "leaf_audit.tsv", "root_hog gene_id protId species mapped_taxon source parent_node parent_event fate selected_by")
        counts = table(stack, r.out / "hog_counts.tsv",
                       "root_hog before_tips before_train before_oma before_test removed_oma removed_test empty_groups_removed unary_collapsed after_tips after_train after_oma after_test final_multitip_tips status unpruned_sha256 pruned_sha256")
        singleton = table(stack, r.out / "singletons.tsv", "root_hog gene_id protId taxon")
        for xml in records(r.xml):
            if xml.tag == SPECIES:
                continue
            hog = xml.attrib["id"]
            path = tree_path(r.raw, hog)
            tree = ete_backend.read(path) if ete_backend else read_tree(path)
            # Read the immutable XML direct-parent relationships whose NHX correspondence was verified.
            groups = ete_backend.groups(tree, xml) if ete_backend else None
            parents, chosen = decisions(r, hog, xml, assignment, aligner, groups)
            tests = parents.keys() & r.test
            removed = tests | chosen.keys()
            stats = Counter()
            if ete_backend:
                root, remaining = ete_backend.prune(tree, removed, stats, parents.keys() - removed)
            else:
                root = prune_node(tree.seed_node, removed, stats)
                require(root is not None, f"Training root became empty: {hog}")
                root.parent_node = None
                tree.seed_node = root
                remaining = [n.label for n in tree.leaf_node_iter()]
                require(len(remaining) == len(set(remaining)) and set(remaining) == parents.keys() - removed,
                        f"{hog}: retained membership differs")
                require(not r.test.intersection(remaining) and all(n.num_child_nodes() != 1 for n in tree),
                        f"{hog}: test tip or unary node survived")
            status = "singleton" if len(remaining) == 1 else "multi_tip"
            final_hash = ""
            if status == "multi_tip":
                dest = tree_path(r.final, hog)
                if ete_backend:
                    ete_backend.write_checked(root, dest)
                else:
                    tree.write(path=str(dest), schema="newick", suppress_rooting=True, suppress_item_comments=False,
                               suppress_leaf_node_labels=False, unquoted_underscores=True)
                    require(signature(read_tree(dest).seed_node) == signature(root), f"{hog}: NHX round trip failed")
                final_hash = sha(dest)
                multi.add(dest.name)
            else:
                gid = remaining[0]
                singleton([hog, gid, *r.train[gid]])
            before = Counter(source(r, gid) for gid in parents)
            after = Counter(source(r, gid) for gid in remaining)
            require(before["train"] == after["train"] and after["oma"] == before["oma"]-len(chosen),
                    f"{hog}: counts inconsistent")
            for gid, parent in parents.items():
                protein, species = r.registry[gid]
                kind = source(r, gid)
                taxon = r.taxon.get(species) if kind == "oma" else species
                fate = ("remove_test" if gid in tests else "remove_sister" if gid in chosen
                        else "exclude_singleton" if status == "singleton" else "retain")
                audit([hog, gid, protein, species, taxon, kind, *parent, fate, ",".join(chosen.get(gid, ()))])
                if kind == "train":
                    require(gid not in r.coverage and gid not in removed, f"Training gene {gid} lost or repeated")
                    r.coverage[gid] = (hog, status)
            counts([hog, len(parents), before["train"], before["oma"], before["test"], len(chosen), len(tests),
                    stats["empty"], stats["unary"], len(remaining), after["train"], after["oma"], after["test"],
                    len(remaining) if status == "multi_tip" else 0, status, sha(path), final_hash])
            totals.update(removed_oma=len(chosen), removed_test=len(tests), empty_groups=stats["empty"],
                          collapsed_unary=stats["unary"], final_tips=len(remaining) if status == "multi_tip" else 0)
        require(filename_set(r.final) == multi, "Final file set differs from multi-tip HOGs")
        require({gid: h for gid, (h, _) in r.coverage.items()} == r.assigned, "Training coverage differs from Stage 01")
        coverage = table(stack, r.out / "training_coverage.tsv", "gene_id protId taxon root_hog status")
        for gid, info in r.train.items():
            hog, status = r.coverage.get(gid, ("", "unassigned"))
            coverage([gid, *info, hog, status])
            totals["training_" + status] += 1
    r.manifest["counts"] = dict(totals, training_registry=len(r.train), final_multi_tip_hogs=len(multi),
                                singleton_hogs=len(r.selected)-len(multi))
    r.checkpoint("pruning, training coverage and final annotated tree round trips")


def shape(node):
    if node.is_leaf():
        return node.label
    event = "S" if node.comments[0].startswith("&&NHX:Ev=0>1") else "D"
    return [event, *(shape(child) for child in node.child_node_iter())]


def verify_fixture(r, ete_backend=None):
    """Hand-specified expected trees and decisions, independent of the pruning implementation."""
    expected_path = r.source.parent / "expected.json"
    oracle = json.loads(expected_path.read_text())
    def rows(name):
        with (r.out / name).open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    require(set(r.selected) == set(oracle["selected_roots"]), "Fixture root cohort")
    require(filename_set(r.final) == {tree_path(r.final, h).name for h in oracle["final_shapes"]},
            "Fixture final root set")
    for h, topology in oracle["final_shapes"].items():
        path = tree_path(r.final, h)
        observed = ete_backend.shape(ete_backend.read(path)) if ete_backend else shape(read_tree(path).seed_node)
        require(observed == topology, f"Fixture topology: {h}")
    require({row["cafa_gene_id"]: [row["decision"], row["selected_oma_gene_id"]]
             for row in rows("pruning_assignments.tsv")} == oracle["decisions"], "Fixture assignments")
    audit = rows("leaf_audit.tsv")
    for action, key in (("remove_test", "removed_test"), ("remove_sister", "removed_oma")):
        require(sorted(x["gene_id"] for x in audit if x["fate"] == action) == sorted(oracle[key]),
                f"Fixture {action}")
    require(sorted(x["gene_id"] for x in rows("singletons.tsv")) == oracle["singletons"], "Fixture singleton IDs")
    require(sorted(x["gene_id"] for x in rows("training_coverage.tsv") if x["status"] == "unassigned")
            == oracle["unassigned"], "Fixture unassigned IDs")
    r.manifest["fixture_expected_sha256"] = sha(expected_path)


def example(r):
    """Show a small changed HOG as two annotated topology expressions with protein names."""
    with (r.out / "hog_counts.tsv").open() as handle:
        candidates = (row for row in csv.DictReader(handle, delimiter="\t")
                      if row["status"] == "multi_tip" and int(row["before_tips"]) <= 40
                      and int(row["before_tips"]) > int(row["after_tips"]))
        row = next(candidates, None)
    if row is None:
        print("No changed multi-tip HOG with at most 40 tips in this run.")
        return
    def labelled(node):
        if node.is_leaf():
            return r.registry[node.label][0]
        event = "S" if node.comments[0].startswith("&&NHX:Ev=0>1") else "D"
        return event + "(" + ", ".join(labelled(c) for c in node.child_node_iter()) + ")"
    hog = row["root_hog"]
    print(hog)
    for label, folder in (("Before", r.raw), ("After", r.final)):
        print(label, labelled(read_tree(tree_path(folder, hog)).seed_node))
    print("OMA sisters removed:", row["removed_oma"], "| Test tips removed:", row["removed_test"],
          "| Unary nodes collapsed:", row["unary_collapsed"])


def finish(r, ete_backend=None):
    if r.mode == "fixture":
        verify_fixture(r, ete_backend)
    r.manifest["table_sha256"] = {p.name: sha(p) for p in sorted(r.out.glob("*.tsv"))}
    r.manifest.update(status="complete", completed=datetime.now(timezone.utc).isoformat())
    r.checkpoint("all applicable checks passed")
    print(json.dumps(r.manifest["counts"], indent=2))
    print("Run directory:", r.out)
