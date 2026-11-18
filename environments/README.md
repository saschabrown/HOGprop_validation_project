# Environments

- FastOMA source and its existing Conda recipe are under [software/fastoma](../software/fastoma). The historical launch activated `fastoma_omamer212`. The recipe still needs an exact environment lock and a verified setup on both machines.
- Stage 02 source uses lxml, NumPy, pandas, Matplotlib, IPython and the external OrthoXML parser. Verify and pin the exact parser implementation before a fresh notebook run.
- Stage 03 software versions are pending the alignment and fitting pilot.
- Stage 04 has an [upstream HOGprop revision](../software/hogprop/README.md). Its full environment lock is pending.
