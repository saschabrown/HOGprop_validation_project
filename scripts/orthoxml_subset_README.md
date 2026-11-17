# OrthoXML subset script

Checkpoint recorded 2026-11-17 before restructuring.

`orthoxml_subset.py` preserves Irene's supplied script and author header. It selects complete root HOGs by a supplied ID list and writes their associated gene registry. No scientific behavior was changed for this checkpoint.

The script requires `lxml` and the external `parser` module exposing `StreamOrthoXMLParser`, `OrthoXMLStreamWriter`, and `process_stream_orthoxml`. That parser implementation is not bundled here. Dependency packaging and notebook integration remain pending; this is not yet a standalone reproducible installation.

The function docstring mentions a `prot_ids` argument that the current signature does not accept. Use root HOG IDs through the existing `-t` option. Syntax was checked for this checkpoint; no extraction was run.

The same checkpoint preserves existing deletions of an empty `notebooks/Untitled.ipynb` and the generated `config/fastoma_slurm_large_resolved.config` snapshot. Both remain available in prior Git history. Parent-workspace analysis notebooks, data and agent support are outside this repository and are not included.
